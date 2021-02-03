# ----------------------------------------------------------------------------
# Copyright (c) 2017-2021, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import q2_quality_control
from qiime2.plugin import (Str, Plugin, Choices, Range, Float, Int, Bool,
                           MetadataColumn, Categorical, Citations, TypeMap,
                           Visualization, TypeMatch)
from q2_types.feature_data import FeatureData, Sequence, Taxonomy
from q2_types.sample_data import SampleData
from q2_types.per_sample_sequences import (
    SequencesWithQuality, PairedEndSequencesWithQuality)
from q2_types.feature_table import FeatureTable, RelativeFrequency
from q2_types.bowtie2 import Bowtie2Index

from .quality_control import (exclude_seqs, evaluate_composition,
                              evaluate_seqs, evaluate_taxonomy)
from ._filter import bowtie2_build, filter_reads


citations = Citations.load('citations.bib', package='q2_quality_control')

plugin = Plugin(
    name='quality-control',
    version=q2_quality_control.__version__,
    website='https://github.com/qiime2/q2-quality-control',
    package='q2_quality_control',
    description=(
        'This QIIME 2 plugin supports methods for assessing and controlling '
        'the quality of feature and sequence data.'),
    short_description=(
        'Plugin for quality control of feature and sequence data.')
)


seq_inputs = {'query_sequences': FeatureData[Sequence],
              'reference_sequences': FeatureData[Sequence]}

seq_inputs_descriptions = {
    'query_sequences': 'Sequences to test for exclusion',
    'reference_sequences': ('Reference sequences to align against feature '
                            'sequences')}

taxa_inputs = {'depth': Int,
               'palette': Str % Choices([
                    'Set1', 'Set2', 'Set3', 'Pastel1', 'Pastel2', 'Paired',
                    'Accent', 'Dark2', 'tab10', 'tab20', 'tab20b', 'tab20c',
                    'viridis', 'plasma', 'inferno', 'magma', 'terrain',
                    'rainbow'])}

taxa_inputs_descriptions = {
    'depth': 'Maximum depth of semicolon-delimited taxonomic ranks to '
             'test (e.g., 1 = root, 7 = species for the greengenes '
             'reference sequence database).',
    'palette': 'Color palette to utilize for plotting.'}

filter_input = {'demultiplexed_sequences': 'The sequences to be trimmed.',
                'database': 'Bowtie2 indexed database.'}
filter_output = {'filtered_sequences': 'The resulting filtered sequences.'}

filter_parameters = {
    'n_threads': Int % Range(1, None),
    'mode': Str % Choices(['local', 'global']),
    'sensitivity': Str % Choices([
        'very-fast', 'fast', 'sensitive', 'very-sensitive']),
    'ref_gap_open_penalty': Int % Range(1, None),
    'ref_gap_ext_penalty': Int % Range(1, None),
    'exclude_seqs': Bool,
}

filter_parameter_descriptions = {
    'n_threads': 'Number of alignment threads to launch.',
    'mode': 'Bowtie2 alignment settings. See bowtie2 manual for more details.',
    'sensitivity': 'Bowtie2 alignment sensitivity. See bowtie2 manual for '
                   'details.',
    'ref_gap_open_penalty': 'Reference gap open penalty.',
    'ref_gap_ext_penalty': 'Reference gap extend penalty.',
    'exclude_seqs': 'Exclude sequences that align to reference. Set this '
                    'option to False to exclude sequences that do not align '
                    'to the reference database.'
}

P_method, P_left_justify, _ = TypeMap({
    (Str % Choices("blast", "blastn-short"), Bool % Choices(False)):
        Visualization,
    (Str % Choices("vsearch"), Bool): Visualization,
})

plugin.methods.register_function(
    function=exclude_seqs,
    inputs=seq_inputs,
    parameters={'method': P_method,
                'perc_identity': Float % Range(0.0, 1.0, inclusive_end=True),
                'evalue': Float,
                'perc_query_aligned': Float,
                'threads': Int % Range(1, None),
                'left_justify': P_left_justify,
                },
    outputs=[('sequence_hits', FeatureData[Sequence]),
             ('sequence_misses', FeatureData[Sequence])],
    input_descriptions=seq_inputs_descriptions,
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
        'left_justify': ('Reject match if the pairwise alignment begins with '
                         'gaps'),
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
        'filtering thresholds as necessary.'),
    citations=[citations['camacho2009blast+']]
)

plugin.visualizers.register_function(
    function=evaluate_composition,
    inputs={'expected_features': FeatureTable[RelativeFrequency],
            'observed_features': FeatureTable[RelativeFrequency]},
    parameters={**taxa_inputs,
                'plot_tar': Bool,
                'plot_tdr': Bool,
                'plot_r_value': Bool,
                'plot_r_squared': Bool,
                'plot_bray_curtis': Bool,
                'plot_jaccard': Bool,
                'plot_observed_features': Bool,
                'plot_observed_features_ratio': Bool,
                'metadata': MetadataColumn[Categorical]},
    input_descriptions={
        'expected_features': 'Expected feature compositions',
        'observed_features': 'Observed feature compositions'},
    parameter_descriptions={
        **taxa_inputs_descriptions,
        'plot_tar': 'Plot taxon accuracy rate (TAR) on score plot. TAR is '
                    'the number of true positive features divided by the '
                    'total number of observed features (TAR = true positives '
                    '/ (true positives + false positives)).',
        'plot_tdr': 'Plot taxon detection rate (TDR) on score plot. TDR is '
                    'the number of true positive features divided by the '
                    'total number of expected features (TDR = true positives '
                    '/ (true positives + false negatives)).',
        'plot_r_value': 'Plot expected vs. observed linear regression r '
                        'value on score plot.',
        'plot_r_squared': 'Plot expected vs. observed linear regression r-'
                          'squared value on score plot.',
        'plot_bray_curtis': 'Plot expected vs. observed Bray-Curtis '
                            'dissimilarity scores on score plot.',
        'plot_jaccard': 'Plot expected vs. observed Jaccard distances scores '
                        'on score plot.',
        'plot_observed_features':
            'Plot observed features count on score plot.',
        'plot_observed_features_ratio':
            'Plot ratio of observed:expected features on score plot.',
        'metadata': 'Optional sample metadata that maps observed_features '
                    'sample IDs to expected_features sample IDs.'},
    name='Evaluate expected vs. observed taxonomic composition of samples',
    description='This visualizer compares the feature composition of pairs of '
        'observed and expected samples containing the same sample ID in two '
        'separate feature tables. Typically, feature composition will consist '
        'of taxonomy classifications or other semicolon-delimited feature '
        'annotations. Taxon accuracy rate, taxon detection rate, and linear '
        'regression scores between expected and observed observations are '
        'calculated at each semicolon-delimited rank, and plots of per-level '
        'accuracy and observation correlations are plotted. A histogram of '
        'distance between false positive observations and the nearest '
        'expected feature is also generated, where distance equals the number '
        'of rank differences between the observed feature and the nearest '
        'common lineage in the expected feature. This visualizer is most '
        'suitable for testing per-run data quality on sequencing runs that '
        'contain mock communities or other samples with known composition. '
        'Also suitable for sanity checks of bioinformatics pipeline '
        'performance.',
    citations=[citations['bokulich2018optimizing']]
)

plugin.visualizers.register_function(
    function=evaluate_seqs,
    inputs=seq_inputs,
    parameters={'show_alignments': Bool},
    input_descriptions=seq_inputs_descriptions,
    parameter_descriptions={
        'show_alignments': 'Option to plot pairwise alignments of query '
                           'sequences and their top hits.'},
    name='Compare query (observed) vs. reference (expected) sequences.',
    description='This action aligns a set of query (e.g., observed) sequences '
        'against a set of reference (e.g., expected) sequences to evaluate '
        'the quality of alignment. The intended use is to align observed '
        'sequences against expected sequences (e.g., from a mock community) '
        'to determine the frequency of mismatches between observed sequences '
        'and the most similar expected sequences, e.g., as a measure of '
        'sequencing/method error. However, any sequences may be provided as '
        'input to generate a report on pairwise alignment quality against '
        'a set of reference sequences.',
    citations=[citations['camacho2009blast+']]
)

plugin.visualizers.register_function(
    function=evaluate_taxonomy,
    inputs={'expected_taxa': FeatureData[Taxonomy],
            'observed_taxa': FeatureData[Taxonomy],
            'feature_table': FeatureTable[RelativeFrequency]},
    parameters={**taxa_inputs,
                'require_exp_ids': Bool,
                'require_obs_ids': Bool,
                'sample_id': Str},
    input_descriptions={
        'expected_taxa': 'Expected taxonomic assignments',
        'observed_taxa': 'Observed taxonomic assignments',
        'feature_table': 'Optional feature table containing relative '
                         'frequency of each feature, used to weight accuracy '
                         'scores by frequency. Must contain all features '
                         'found in expected and/or observed taxa. Features '
                         'found in the table but not the expected/observed '
                         'taxa will be dropped prior to analysis.'},
    parameter_descriptions={
        **taxa_inputs_descriptions,
        'require_obs_ids': 'Require that all features found in expected taxa '
                           'must be found in observed taxa or raise error.',
        'require_exp_ids': 'Require that all features found in observed taxa '
                           'must be found in expected taxa or raise error.',
        'sample_id': 'Optional sample ID to use for extracting frequency data '
                     'from feature table, and for labeling accuracy results. '
                     'If no sample_id is provided, feature frequencies are '
                     'derived from the sum of all samples present in the '
                     'feature table.'},
    name='Evaluate expected vs. observed taxonomic assignments',
    description='This visualizer compares a pair of observed and expected '
        'taxonomic assignments to calculate precision, recall, and F-measure '
        'at each taxonomic level, up to maximum level specified by the depth '
        'parameter. These metrics are calculated at each semicolon-delimited '
        'rank. This action is useful for comparing the accuracy of taxonomic '
        'assignment, e.g., between different taxonomy classifiers or other '
        'bioinformatics methods. Expected taxonomies should be derived from '
        'simulated or mock community sequences that have known taxonomic '
        'affiliations.',
    citations=[citations['bokulich2018optimizing']]
)

T = TypeMatch([SequencesWithQuality, PairedEndSequencesWithQuality])
plugin.methods.register_function(
    function=filter_reads,
    inputs={'demultiplexed_sequences': SampleData[T],
            'database': Bowtie2Index},
    parameters=filter_parameters,
    outputs=[('filtered_sequences', SampleData[T])],
    input_descriptions=filter_input,
    parameter_descriptions=filter_parameter_descriptions,
    output_descriptions=filter_output,
    name='Filter demultiplexed sequences by alignment to reference database.',
    description=(
        'Filter out (or keep) demultiplexed single- or paired-end sequences '
        'that align to a reference database, using bowtie2 and samtools. This '
        'method can be used to filter out human DNA sequences and other '
        'contaminants in any FASTQ sequence data (e.g., shotgun genome or '
        'amplicon sequence data), or alternatively (when exclude_seqs is '
        'False) to only keep sequences that do align to the reference.'),
    citations=[citations['langmead2012fast'], citations['heng2009samtools']]
)

plugin.methods.register_function(
    function=bowtie2_build,
    inputs={'sequences': FeatureData[Sequence]},
    parameters={'n_threads': Int % Range(1, None)},
    outputs=[('database', Bowtie2Index)],
    input_descriptions={
        'sequences': 'Reference sequences used to build bowtie2 index.'},
    parameter_descriptions={'n_threads': 'Number of threads to launch'},
    output_descriptions={'database': 'Bowtie2 index.'},
    name='Build bowtie2 index from reference sequences.',
    description='Build bowtie2 index from reference sequences.',
    citations=[citations['langmead2012fast']]
)
