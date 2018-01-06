# ----------------------------------------------------------------------------
# Copyright (c) 2017-2018, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import qiime2
from q2_types.feature_data import DNAFASTAFormat
from q2_types.feature_data._transformer import _dnafastaformats_to_series
import pandas as pd
from ._blast import _search_seqs
from ._utilities import _evaluate_composition, _visualize
from ._evaluate_seqs import _evaluate_seqs


def exclude_seqs(query_sequences: DNAFASTAFormat,
                 reference_sequences: DNAFASTAFormat, method: str='blast',
                 perc_identity: float=0.97, evalue: float=None,
                 perc_query_aligned: float=0.97, threads: str=1
                 ) -> (pd.Series, pd.Series):

    # BLAST query seqs vs. ref db of contaminants (or targets)
    hit_ids = _search_seqs(
        query_sequences, reference_sequences, evalue=evalue,
        perc_identity=perc_identity, threads=threads,
        perc_query_aligned=perc_query_aligned, method=method)

    # convert query_sequences to series for filtering
    query_series = _dnafastaformats_to_series(query_sequences)

    # if no hits are in hit_ids, return empty hits and query_series as misses
    if len(hit_ids) < 1:
        hits_seqs = pd.Series()
        return hits_seqs, query_series
    # if all query seqs are hits, return query_series as hits and empty misses
    elif len(hit_ids) == len(query_series):
        misses_seqs = pd.Series()
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
        return pd.Series(hits_seqs), pd.Series(misses_seqs)


def evaluate_composition(
        output_dir: str, expected_features: pd.DataFrame,
        observed_features: pd.DataFrame, depth: int=7, palette: str='Set1',
        plot_tar: bool=True, plot_tdr: bool=True, plot_r_value: bool=False,
        plot_r_squared: bool=True, plot_observed_features: bool=False,
        plot_observed_features_ratio: bool=True,
        metadata: qiime2.CategoricalMetadataColumn=None) -> None:

    # results, fn_features, misclassifications, underclassifications,
    # composition_regression, score_plot, mismatch_histogram
    results = _evaluate_composition(
        expected_features, observed_features, depth=depth, palette=palette,
        metadata=metadata, plot_tar=plot_tar, plot_tdr=plot_tdr,
        plot_r_value=plot_r_value, plot_r_squared=plot_r_squared,
        plot_observed_features=plot_observed_features,
        plot_observed_features_ratio=plot_observed_features_ratio)

    _visualize(output_dir, 'Feature evaluation results',
               'evaluate_composition', *results)


def evaluate_seqs(output_dir: str, query_sequences: DNAFASTAFormat,
                  reference_sequences: DNAFASTAFormat,
                  show_alignments: bool=False) -> None:

    results, alignments, mismatch_histogram = _evaluate_seqs(
        query_sequences, reference_sequences, show_alignments)

    _visualize(output_dir, 'Sequence evaluation results',
               'evaluate_seqs', results=results,
               false_negative_features=None,
               misclassifications=None, underclassifications=None,
               composition_regression=None, score_plot=None,
               mismatch_histogram=mismatch_histogram, alignments=alignments)
