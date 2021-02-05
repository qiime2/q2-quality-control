# ----------------------------------------------------------------------------
# Copyright (c) 2017-2021, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import qiime2
import biom
from q2_types.feature_data import DNAFASTAFormat
from q2_types.feature_data._transformer import _dnafastaformats_to_series
import pandas as pd

from ._blast import _search_seqs
from ._utilities import (
    _evaluate_composition, _visualize, _pointplot_multiple_y)
from ._evaluate_seqs import _evaluate_seqs
from ._evaluate_taxonomy import _evaluate_taxonomy

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
    query_series = _dnafastaformats_to_series(query_sequences)

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
