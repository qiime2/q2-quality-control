# ----------------------------------------------------------------------------
# Copyright (c) 2017, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import q2_quality_control
from qiime2.plugin import Str, Plugin, Choices, Range, Float, Int
from q2_types.feature_data import FeatureData, Sequence
from .quality_control import exclude_seqs


plugin = Plugin(
    name='quality-control',
    version=q2_quality_control.__version__,
    website="https://github.com/qiime2/q2-quality-control",
    package='q2_quality_control',
    description=(
        'This QIIME 2 plugin supports methods for assessing and controlling '
        'the quality of feature and sequence data.'),
    short_description=(
        'Plugin for quality control of feature and sequence data.')
)


plugin.methods.register_function(
    function=exclude_seqs,
    inputs={'feature_sequences': FeatureData[Sequence],
            'reference_sequences': FeatureData[Sequence]},
    parameters={'method': Str % Choices(["blast", "vsearch", "blastn-short"]),
                'perc_identity': Float % Range(0.0, 1.0, inclusive_end=True),
                'evalue': Float,
                'perc_query_aligned': Float,
                'threads': Int % Range(1, None)},
    outputs=[('sequence_hits', FeatureData[Sequence]),
             ('sequence_misses', FeatureData[Sequence])],
    input_descriptions={
        'feature_sequences': 'Sequences to test for exclusion',
        'reference_sequences': ('Reference sequences to align against feature '
                                'sequences')},
    parameter_descriptions={
        'method': ('Alignment method to use for matching feature sequences '
                   'against reference sequences'),
        'perc_identity': ('Reject match if percent identity to reference is '
                          'lower. Must be in range [0.0, 1.0]'),
        'evalue': ('BLAST expectation (E) value threshold for saving hits. '
                   'Reject if E value is higher than threshold. This '
                   'threshold is disabled by default.'),
        'perc_query_aligned': (
            'Percent of query sequence that must align to reference in order '
            'to be accepted as a hit.'),
        'threads': (
            'Number of jobs to execute. Only applies to vsearch method.'),
    },
    output_descriptions={
        'sequence_hits': (
            'Subset of feature sequences that align to reference sequences'),
        'sequence_misses': (
            'Subset of feature sequences that do not align to reference '
            'sequences')
    },
    name='Exclude sequences by alignment',
    description=(
        'This method aligns feature sequences to a set of reference sequences '
        'to identify sequences that hit/miss the reference within a specified '
        'perc_identity, evalue, and perc_query_aligned. This method could '
        'be used to define a positive filter, e.g., extract only feature '
        'sequences that align to a certain clade of bacteria; or to define a '
        'negative filter, e.g., identify sequences that align to contaminant '
        'or human DNA sequences that should be excluded from subsequent '
        'analyses. Note that filtering is performed based on the '
        'perc_identity, perc_query_aligned, and evalue thresholds (the '
        'latter only if method==BLAST and an evalue is set). Set '
        'perc_identity==0 and/or perc_query_aligned==0 to disable these '
        'filtering thresholds as necessary.')
)
