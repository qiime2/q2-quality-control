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
from q2_types.feature_table import FeatureTable, RelativeFrequency
from .quality_control import (
    exclude_seqs, evaluate_taxonomic_composition)
from ._utilities import _results_columns


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
    parameters={'method': Str % Choices(["blast", "vsearch"]),
                'perc_identity': Float % Range(0.0, 1.0, inclusive_end=True),
                'evalue': Float,
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
                          'lower. (or if E value is higher than E value '
                          'threshold.)  Must be in range [0.0, 1.0]'),
        'evalue': ('BLAST expectation (E) value threshold for saving hits. '
                   'Reject if E value is higher than threshold (or '
                   'perc_identity is lower than identity threshold.)'),
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
        'percent identity. This method could be used to define a positive '
        'filter, e.g., extract only feature sequences that align to a certain '
        'clade of bacteria; or to define a negative filter, e.g., identify '
        'sequences that align to contaminant or human DNA sequences that '
        'should be excluded from subsequent analyses. Note that filtering is '
        'performed based on both the perc_identity threshold and E value (the '
        'latter only if method==BLAST). Set E value to a higher value or '
        'perc_identity==0 to disable one or the other filtering method.')
)

plugin.visualizers.register_function(
    function=evaluate_taxonomic_composition,
    inputs={'expected_features': FeatureTable[RelativeFrequency],
            'observed_features': FeatureTable[RelativeFrequency]},
    parameters={'depth': Int,
                'palette': Str % Choices([
                    'Set1', 'Set2', 'Set3', 'Pastel1', 'Pastel2', 'Paired',
                    'Accent', 'Dark2', 'tab10', 'tab20', 'tab20b', 'tab20c',
                    'viridis', 'plasma', 'inferno', 'magma', 'terrain',
                    'rainbow']),
                'yvals': Str},
    input_descriptions={
        'expected_features': 'Expected feature compositions',
        'observed_features': 'Observed feature compositions'},
    parameter_descriptions={
        'depth': ("Maximum depth of semicolon-delimited taxonomic ranks to "
                  "test (e.g., 1 = root, 7 = species for the greengenes "
                  "reference sequence database)."),
        'palette': "Color palette to utilize for plotting.",
        'yvals': ("Comma-separated list of accuracy metrics to plot on score "
                  "plots. Enclose in quotes and do not type spaces between "
                  "values. Must be one more more of [{0}].".format(
                    ','.join(_results_columns())))},
    name='Evaluate expected vs. observed taxonomic composition of samples',
    description=(
        "This visualizer compares the feature composition of pairs of "
        "observed and expected samples containing the same sample ID in two "
        "separate feature tables. Typically, this feature composition will be "
        "taxonomy classifications or other semicolon-delimited feature "
        "annotations. Taxon accuracy rate, taxon detection rate, and linear "
        "regression scores between expected and observed observations are "
        "calculated at each semicolon-delimited rank, and plots of per-level "
        "accuracy and observation correlations are plotted. A histogram of "
        "distance between false positive observations and the nearest "
        "expected feature is also generated, where distance equals the number "
        "of rank differences between the observed feature and the nearest "
        "common lineage in the expected feature. This visualizer is most "
        "suitable for testing per-run data quality on sequencing runs that "
        "contain mock communities or other samples with known composition. "
        "Also suitable for sanity checks of bioinformatic pipeline "
        "performance.")
)
