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
from q2_feature_table import filter_seqs
import pandas as pd


def exclude_seqs(feature_sequences: DNAFASTAFormat,
                 reference_sequences: DNAFASTAFormat, method='blast',
                 perc_identity: float=0.97, threads: str=1):

    # q2-feature-classifier requires a reference taxonomy but we don't really
    # care about assigning taxonomy here — just identifying hits/misses.
    # So first we generate a fake taxonomy file to avoid requiring a real
    # taxonomy.
    reference_taxonomy = qiime2.Artifact.import_data(
        "FeatureData[Sequence]", reference_sequences).view(pd.Series).index
    reference_taxonomy = pd.Series(
        reference_taxonomy, index=reference_taxonomy, name='Taxon')

    # BLAST query seqs vs. ref db of contaminants (or targets)
    if method == 'blast':
        # blast uses perc_identity as a percentage, vsearch as a decimal
        perc_identity = perc_identity * 100
        res = qfc._blast.classify_consensus_blast(
            reference_sequences, feature_sequences, reference_taxonomy,
            maxaccepts=1, perc_identity=perc_identity)
    elif method == 'vsearch':
        res = qfc._vsearch.classify_consensus_vsearch(
            reference_sequences, feature_sequences, reference_taxonomy,
            maxaccepts=1, perc_identity=perc_identity, threads=threads)

    # convert feature_sequences to series for filtering
    query_series = qiime2.Artifact.import_data(
        "FeatureData[Sequence]", feature_sequences).view(pd.Series)

    # filter seqs from seq file
    hits_seqs = filter_seqs(query_series, res, exclude_ids=True)
    misses_seqs = filter_seqs(query_series, res, exclude_ids=False)

    # output hits/rejects
    return hits_seqs, misses_seqs
