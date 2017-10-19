# ----------------------------------------------------------------------------
# Copyright (c) 2017, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import qiime2
from q2_types.feature_data import DNAFASTAFormat
from q2_types.feature_data._transformer import _dnafastaformats_to_series
import pandas as pd
from ._utilities import _evaluate_taxonomic_composition, _visualize
from ._blast import _search_seqs


def exclude_seqs(feature_sequences: DNAFASTAFormat,
                 reference_sequences: DNAFASTAFormat, method='blast',
                 perc_identity: float=0.97, evalue: float=None,
                 perc_query_aligned: float=0.97, threads: str=1
                 ) -> (pd.Series, pd.Series):

    # BLAST query seqs vs. ref db of contaminants (or targets)
    hit_ids = _search_seqs(
        feature_sequences, reference_sequences, evalue=evalue,
        perc_identity=perc_identity, threads=threads,
        perc_query_aligned=perc_query_aligned, method=method)

    # convert feature_sequences to series for filtering
    query_series = _dnafastaformats_to_series(feature_sequences)

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


def evaluate_taxonomic_composition(
        output_dir: str, expected_features: pd.DataFrame,
        observed_features: pd.DataFrame, depth: int=7, palette: str='Set1',
        yvals: str='TAR,TDR,R,Observed / Expected Taxa',
        metadata: qiime2.MetadataCategory=None) -> None:

    # results, fn_features, misclassifications, underclassifications,
    # composition_regression, score_plot, mismatch_histogram
    results = _evaluate_taxonomic_composition(
        expected_features, observed_features, depth=depth, palette=palette,
        yvals=yvals, metadata=metadata)

    _visualize(output_dir, *results)
