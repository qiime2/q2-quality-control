# ----------------------------------------------------------------------------
# Copyright (c) 2017, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import qiime2
import q2_feature_classifier as qfc
from q2_types.feature_data import DNAFASTAFormat
from q2_types.feature_data._transformer import _dnafastaformats_to_series
from q2_feature_table import filter_seqs
import pandas as pd
from ._utilities import _evaluate_taxonomic_composition, _visualize


def exclude_seqs(feature_sequences: DNAFASTAFormat,
                 reference_sequences: DNAFASTAFormat, method='blast',
                 perc_identity: float=0.97, threads: str=1,
                 ) -> (pd.Series, pd.Series):

    # q2-feature-classifier requires a reference taxonomy but we don't really
    # care about assigning taxonomy here — just identifying hits/misses.
    # So first we generate a fake taxonomy file to avoid requiring a real
    # taxonomy.
    reference_taxonomy = _dnafastaformats_to_series(reference_sequences).index
    reference_taxonomy = pd.Series(
        reference_taxonomy, index=reference_taxonomy, name='Taxon')

    # BLAST query seqs vs. ref db of contaminants (or targets)
    if method == 'blast':
        res = qfc._blast.classify_consensus_blast(
            feature_sequences, reference_sequences, reference_taxonomy,
            maxaccepts=1, perc_identity=perc_identity)
    elif method == 'vsearch':
        res = qfc._vsearch.classify_consensus_vsearch(
            feature_sequences, reference_sequences, reference_taxonomy,
            maxaccepts=1, perc_identity=perc_identity, threads=threads)

    # convert feature_sequences to series for filtering
    query_series = _dnafastaformats_to_series(feature_sequences)

    # filter seqs from seq file
    res_md = qiime2.Metadata(res)
    hits_seqs = filter_seqs(
        query_series, res_md, where="Taxon='Unassigned'", exclude_ids=True)
    misses_seqs = filter_seqs(
        query_series, res_md, where="Taxon='Unassigned'", exclude_ids=False)

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
