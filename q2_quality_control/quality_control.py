# ----------------------------------------------------------------------------
# Copyright (c) 2017-2022, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import biom
from q2_types.feature_data import DNAFASTAFormat
import pandas as pd
import os
import tempfile
import subprocess
import qiime2.util

from qiime2.plugin.util import transform
from ._blast import _search_seqs
from ._utilities import (
    _evaluate_composition, _visualize, _pointplot_multiple_y)
from ._evaluate_seqs import _evaluate_seqs
from ._evaluate_taxonomy import _evaluate_taxonomy
from ._utilities import _run_command
from ._stats import DecontamScoreFormat

left_justify_supported_methods = {'vsearch'}


def exclude_seqs(query_sequences: DNAFASTAFormat,
                 reference_sequences: DNAFASTAFormat, method: str = 'blast',
                 perc_identity: float = 0.97, evalue: float = None,
                 perc_query_aligned: float = 0.97, threads: str = 1,
                 left_justify: bool = False,
                 ) -> (pd.Series, pd.Series):

    if left_justify and (method not in left_justify_supported_methods):
        raise ValueError("Enabling left_justify is not compatible with "
                         "method=%r, check the documentation for valid "
                         "combinations" % (method, ))
    # BLAST query seqs vs. ref db of contaminants (or targets)
    hit_ids = _search_seqs(
        query_sequences, reference_sequences, evalue=evalue,
        perc_identity=perc_identity, threads=threads,
        perc_query_aligned=perc_query_aligned, method=method,
        left_justify=left_justify)

    # convert query_sequences to series for filtering
    query_series = query_sequences.view(pd.Series)

    # if no hits are in hit_ids, return empty hits and query_series as misses
    if len(hit_ids) < 1:
        hits_seqs = pd.Series(dtype='string')
        return hits_seqs, query_series
    # if all query seqs are hits, return query_series as hits and empty misses
    elif len(hit_ids) == len(query_series):
        misses_seqs = pd.Series(dtype='string')
        return query_series, misses_seqs
    # otherwise filter seqs from seq file
    else:
        hits_seqs = {}
        misses_seqs = {}
        for seq_id, seq in query_series.items():
            seq = str(seq)
            if seq_id in hit_ids:
                hits_seqs[seq_id] = seq
            else:
                misses_seqs[seq_id] = seq
        return (pd.Series(hits_seqs, dtype='string'),
                pd.Series(misses_seqs, dtype='string'))


def evaluate_composition(
        output_dir: str, expected_features: pd.DataFrame,
        observed_features: pd.DataFrame, depth: int = 7, palette: str = 'Set1',
        plot_tar: bool = True, plot_tdr: bool = True,
        plot_r_value: bool = False, plot_r_squared: bool = True,
        plot_bray_curtis: bool = False, plot_jaccard: bool = False,
        plot_observed_features: bool = False,
        plot_observed_features_ratio: bool = True,
        metadata: qiime2.CategoricalMetadataColumn = None) -> None:

    # results, fn_features, misclassifications, underclassifications,
    # composition_regression, score_plot, mismatch_histogram
    results = _evaluate_composition(
        expected_features, observed_features, depth=depth, palette=palette,
        metadata=metadata, plot_tar=plot_tar, plot_tdr=plot_tdr,
        plot_r_value=plot_r_value, plot_r_squared=plot_r_squared,
        plot_bray_curtis=plot_bray_curtis, plot_jaccard=plot_jaccard,
        plot_observed_features=plot_observed_features,
        plot_observed_features_ratio=plot_observed_features_ratio)

    _visualize(output_dir, 'Feature evaluation results',
               'evaluate_composition', *results)


def evaluate_seqs(output_dir: str, query_sequences: DNAFASTAFormat,
                  reference_sequences: DNAFASTAFormat,
                  show_alignments: bool = False) -> None:

    results, alignments, mismatch_histogram = _evaluate_seqs(
        query_sequences, reference_sequences, show_alignments)

    _visualize(output_dir, 'Sequence evaluation results',
               'evaluate_seqs', results=results,
               false_negative_features=None,
               misclassifications=None, underclassifications=None,
               composition_regression=None, score_plot=None,
               mismatch_histogram=mismatch_histogram, alignments=alignments)


def evaluate_taxonomy(output_dir: str, expected_taxa: pd.DataFrame,
                      observed_taxa: pd.DataFrame, depth: int,
                      palette: str = 'Set1', require_exp_ids: bool = True,
                      require_obs_ids: bool = True,
                      feature_table: biom.Table = None, sample_id: str = None
                      ) -> None:
    prf = _evaluate_taxonomy(expected_taxa, observed_taxa, require_exp_ids,
                             require_obs_ids, feature_table, sample_id,
                             level_range=range(0, depth))

    score_plot = _pointplot_multiple_y(
        prf, xval='level', yvals=['Precision', 'Recall', 'F-measure'],
        palette=palette)

    _visualize(output_dir, 'Taxonomic accuracy results',
               'evaluate_taxonomy', results=prf,
               false_negative_features=None,
               misclassifications=None, underclassifications=None,
               composition_regression=None, score_plot=score_plot,
               mismatch_histogram=None, alignments=None)


# decontam added code
_WHOLE_NUM = (lambda x: x >= 0, 'non-negative')
_PER_NUM = (lambda x: 1 >= x >= 0, 'between 0 and 1')
_COL_STR = (lambda x: x in {'column_name', 'column_number'},
            'sample_name or column_name or column_number')
_DECON_METHOD_STR = (lambda x: x in {'frequency', 'prevalence', 'combined'},
                     'freqeuncy, prevalence, combined')
_BOOLEAN = (lambda x: type(x) is bool, 'True or False')
# Better to choose to skip, than to implicitly ignore things that KeyError
_SKIP = (lambda x: True, '')
_valid_inputs = {
    'table': _SKIP,
    'metadata': _SKIP,
    'threshold': _PER_NUM,
    'method': _DECON_METHOD_STR,
    'freq_concentration_column': _SKIP,
    'prev_control_column': _SKIP,
    'prev_control_indicator': _SKIP,
}


def _check_inputs(**kwargs):
    for param, arg in kwargs.items():
        check_is_valid, explanation = _valid_inputs[param]
        if not check_is_valid(arg):
            raise ValueError('Argument to %r was %r, should be %s.'
                             % (param, arg, explanation))


def _decontam_identify_helper(track_fp, method):
    df = pd.read_csv(track_fp, sep='\t', index_col=0)
    df.index.name = '#OTU ID'
    # removes last column containing true/false information from the dataframe
    df = df.drop(df.columns[[(len(df.columns)-1)]], axis=1)
    if method == 'combined':
        df = df.fillna(0)

    # removes all columns that are completely empty
    temp_transposed_table = df.transpose()
    temp_transposed_table = temp_transposed_table.dropna()
    df = temp_transposed_table.transpose()
    metadata = transform(df, from_type=pd.DataFrame,
                         to_type=DecontamScoreFormat)

    return metadata


def decontam_identify(table: pd.DataFrame,
                      metadata: qiime2.Metadata,
                      method: str = 'prevalence',
                      freq_concentration_column: str = 'NULL',
                      prev_control_column: str = 'NULL',
                      prev_control_indicator: str = 'NULL'
                      ) -> (DecontamScoreFormat):
    _check_inputs(**locals())
    with tempfile.TemporaryDirectory() as temp_dir_name:
        track_fp = os.path.join(temp_dir_name, 'track.tsv')
        ASV_dest = os.path.join(temp_dir_name, 'temp_ASV_table.csv')
        transposed_table = table.transpose()
        transposed_table.to_csv(os.path.join(ASV_dest))
        metadata = metadata.to_dataframe()
        meta_dest = os.path.join(temp_dir_name, 'temp_metadata.csv')
        metadata.to_csv(os.path.join(meta_dest))

        cmd = ['run_decontam.R',
               '--asv_table_path', str(ASV_dest),
               '--threshold', str(0.1),
               '--decon_method', method,
               '--output_track', track_fp,
               '--meta_table_path', str(meta_dest),
               '--freq_con_column', str(freq_concentration_column),
               '--prev_control_or_exp_sample_column',
               str(prev_control_column),
               '--prev_control_sample_indicator',
               str(prev_control_indicator)]
        try:
            _run_command(cmd)
        except subprocess.CalledProcessError as e:
            if e.returncode == 2:
                raise ValueError("There was an issue running "
                                 "run_decontam.R please check your inputs")
            else:
                raise Exception("An error was encountered "
                                "while running Decontam in R "
                                "(return code %d), please inspect stdout"
                                " and stderr to learn more." % e.returncode)
        return _decontam_identify_helper(track_fp, method)


def decontam_remove(decontam_scores: qiime2.Metadata,
                    table: pd.DataFrame,
                    threshold: float = 0.1
                    ) -> (biom.Table):
    with tempfile.TemporaryDirectory() as temp_dir_name:
        df = decontam_scores.to_dataframe()
        df.loc[(df['p'].astype(float) <= threshold),
               'contaminant_seq'] = 'True'
        df.loc[(df['p'].astype(float) > threshold),
               'contaminant_seq'] = 'False'
        df = df[df.contaminant_seq == 'True']
        remove_these = df.index
        for bad_seq in list(remove_these):
            table = (table[table.index != bad_seq])
        output = os.path.join(temp_dir_name, 'temp.tsv.biom')
        temp_transposed_table = table.transpose()
        temp_transposed_table.to_csv(output, sep="\t")
        with open(output) as fh:
            no_contam_table = biom.Table.from_tsv(fh, None, None, None)
        return no_contam_table
