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
from ._blast import _blast_seqs


def exclude_seqs(feature_sequences: DNAFASTAFormat,
                 reference_sequences: DNAFASTAFormat, method='blast',
                 perc_identity: float=0.97, evalue: float=0.001,
                 threads: str=1) -> (pd.Series, pd.Series):

    # BLAST query seqs vs. ref db of contaminants (or targets)
    res = _blast_seqs(
        feature_sequences, reference_sequences, evalue=evalue,
        perc_identity=perc_identity, threads=threads, method=method)

    # convert feature_sequences to series for filtering
    query_series = _dnafastaformats_to_series(feature_sequences)

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
