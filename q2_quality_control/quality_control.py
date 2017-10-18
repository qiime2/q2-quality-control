# ----------------------------------------------------------------------------
# Copyright (c) 2017, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import qiime2
# import q2_feature_classifier as qfc
from q2_types.feature_data import DNAFASTAFormat
from q2_types.feature_data._transformer import _dnafastaformats_to_series
from q2_feature_table import filter_seqs
import pandas as pd
from ._utilities import _evaluate_taxonomic_composition, _visualize
from ._blast import _search_seqs


def exclude_seqs(feature_sequences: DNAFASTAFormat,
                 reference_sequences: DNAFASTAFormat, method='blast',
                 perc_identity: float=0.97, evalue: float=None,
                 perc_query_aligned: float=0.97, threads: str=1
                 ) -> (pd.Series, pd.Series):

    # BLAST query seqs vs. ref db of contaminants (or targets)
    res = _search_seqs(
        feature_sequences, reference_sequences, evalue=evalue,
        perc_identity=perc_identity, threads=threads,
        perc_query_aligned=perc_query_aligned, method=method)

    # convert feature_sequences to series for filtering
    query_series = _dnafastaformats_to_series(feature_sequences)

    # if no hits are in res, early return empty hits and query_series as misses
    if len(res) < 1:
        hits_seqs = pd.Series()
        return hits_seqs, query_series
    # if all query seqs are hits, return query_series as hits and empty misses
    elif len(res) == len(query_series):
        misses_seqs = pd.Series()
        return query_series, misses_seqs

    # filter seqs from seq file
    res_md = qiime2.Metadata(res)
    hits_seqs = filter_seqs(query_series, res_md, exclude_ids=False)
    misses_seqs = filter_seqs(query_series, res_md, exclude_ids=True)

    # output hits/rejects
    return hits_seqs, misses_seqs


def evaluate_taxonomic_composition(
        output_dir: str, expected_features: pd.DataFrame,
        observed_features: pd.DataFrame, depth: int=7, palette: str='Set1',
        yvals: str='TAR,TDR,R,Observed / Expected Taxa',
        ) -> None:

    # results, fn_features, misclassifications, underclassifications,
    # composition_regression, score_plot, mismatch_histogram
    results = _evaluate_taxonomic_composition(
        expected_features, observed_features, depth=depth, palette=palette,
        yvals=yvals)

    _visualize(output_dir, *results)
