# ----------------------------------------------------------------------------
# Copyright (c) 2017-2022, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import gzip
import itertools
import unittest

from q2_types.per_sample_sequences import (
    SingleLanePerSampleSingleEndFastqDirFmt,
    SingleLanePerSamplePairedEndFastqDirFmt,
    FastqGzFormat,
)
from qiime2 import Artifact
from .test_quality_control import QualityControlTestsBase

import numpy as np
import numpy.testing as npt
import pandas as pd
import qiime2
import biom
import os
import tempfile
from qiime2.plugin.testing import TestPluginBase
from qiime2.plugin.util import transform
from q2_quality_control.quality_control import (decontam_identify,
                                                decontam_remove)
from q2_quality_control._stats import DecontamScoreFormat


class TestBowtie2Build(QualityControlTestsBase):

    # This test just makes sure this runs without error, which will include
    # validating the contents.
    def test_build(self):
        genome = Artifact.load(self.get_data_path('sars2-genome.qza'))
        self.plugin.methods['bowtie2_build'](genome)


seq_ids_that_map = ['SARS2:6:73:941:1973#', 'SARS2:6:73:231:3321#',
                    'SARS2:6:73:233:3421#', 'SARS2:6:73:552:2457#',
                    'SARS2:6:73:567:7631#']
seq_id_that_does_not_map = 'SARS2:6:73:356:9806#'


class TestFilterSingle(QualityControlTestsBase):

    def setUp(self):
        super().setUp()

        self.demuxed_art = Artifact.load(self.get_data_path('single-end.qza'))
        self.indexed_genome = Artifact.load(
            self.get_data_path('sars2-indexed.qza'))

    def test_filter_single_exclude_seqs(self):
        obs_art, = self.plugin.methods['filter_reads'](
            self.demuxed_art, self.indexed_genome, exclude_seqs=True)
        obs = obs_art.view(SingleLanePerSampleSingleEndFastqDirFmt)
        obs_seqs = obs.sequences.iter_views(FastqGzFormat)
        for _, obs_fp in obs_seqs:
            with gzip.open(str(obs_fp), 'rt') as obs_fh:
                self.assertNotEqual(len(obs_fh.readlines()), 0)
                # Iterate over expected and observed reads, side-by-side
                for records in itertools.zip_longest(*[obs_fh] * 4):
                    (obs_seq_h, obs_seq, _, obs_qual) = records
                    # Make sure seqs that map to genome were removed
                    obs_id = obs_seq_h.strip('@/012\n')
                    self.assertTrue(obs_id not in seq_ids_that_map)
                    self.assertTrue(obs_id in seq_id_that_does_not_map)

    def test_filter_single_keep_seqs(self):
        obs_art, = self.plugin.methods['filter_reads'](
            self.demuxed_art, self.indexed_genome, exclude_seqs=False)
        obs = obs_art.view(SingleLanePerSampleSingleEndFastqDirFmt)
        obs_seqs = obs.sequences.iter_views(FastqGzFormat)
        for _, obs_fp in obs_seqs:
            with gzip.open(str(obs_fp), 'rt') as obs_fh:
                self.assertNotEqual(len(obs_fh.readlines()), 0)
                # Iterate over expected and observed reads, side-by-side
                for records in itertools.zip_longest(*[obs_fh] * 4):
                    (obs_seq_h, obs_seq, _, obs_qual) = records
                    # Make sure seqs that do not map to genome were removed
                    obs_id = obs_seq_h.strip('@/012\n')
                    self.assertTrue(obs_id in seq_ids_that_map)
                    self.assertTrue(obs_id not in seq_id_that_does_not_map)


class TestFilterPaired(QualityControlTestsBase):

    def setUp(self):
        super().setUp()

        self.demuxed_art = Artifact.load(self.get_data_path('paired-end.qza'))
        self.indexed_genome = Artifact.load(
            self.get_data_path('sars2-indexed.qza'))

    def test_filter_paired_exclude_seqs(self):
        obs_art, = self.plugin.methods['filter_reads'](
            self.demuxed_art, self.indexed_genome, exclude_seqs=True)
        obs = obs_art.view(SingleLanePerSamplePairedEndFastqDirFmt)
        obs_seqs = obs.sequences.iter_views(FastqGzFormat)
        for _, obs_fp in obs_seqs:
            with gzip.open(str(obs_fp), 'rt') as obs_fh:
                self.assertNotEqual(len(obs_fh.readlines()), 0)
                # Iterate over expected and observed reads, side-by-side
                for records in itertools.zip_longest(*[obs_fh] * 4):
                    (obs_seq_h, obs_seq, _, obs_qual) = records
                    # Make sure seqs that map to genome were removed
                    obs_id = obs_seq_h.strip('@/012\n')
                    self.assertTrue(obs_id not in seq_ids_that_map)
                    self.assertTrue(obs_id in seq_id_that_does_not_map)

    def test_filter_paired_keep_seqs(self):
        obs_art, = self.plugin.methods['filter_reads'](
            self.demuxed_art, self.indexed_genome, exclude_seqs=False)
        obs = obs_art.view(SingleLanePerSamplePairedEndFastqDirFmt)
        obs_seqs = obs.sequences.iter_views(FastqGzFormat)
        for _, obs_fp in obs_seqs:
            with gzip.open(str(obs_fp), 'rt') as obs_fh:
                self.assertNotEqual(len(obs_fh.readlines()), 0)
                # Iterate over expected and observed reads, side-by-side
                for records in itertools.zip_longest(*[obs_fh] * 4):
                    (obs_seq_h, obs_seq, _, obs_qual) = records
                    # Make sure seqs that do not map to genome were removed
                    obs_id = obs_seq_h.strip('@/012\n')
                    self.assertTrue(obs_id in seq_ids_that_map)
                    self.assertTrue(obs_id not in seq_id_that_does_not_map)


# Decontam tests


class TestIdentify(TestPluginBase):
    package = 'q2_quality_control.tests'

    def setUp(self):
        super().setUp()
        table = qiime2.Artifact.load(
            self.get_data_path('expected/decon_default_ASV_table.qza'))
        self.asv_table = table.view(qiime2.Metadata).to_dataframe()
        self.metadata_input = qiime2.Metadata.load(
            self.get_data_path('expected/test_metadata.tsv'))

    def test_prevalence(self):
        exp_table = pd.read_csv(
            self.get_data_path('expected/prevalence-score-table.tsv'),
            sep='\t', index_col=0)
        temp_transposed_table = exp_table.transpose()
        temp_transposed_table = temp_transposed_table.dropna()
        exp_table = temp_transposed_table.transpose()
        output_feature_table = decontam_identify(
            asv_or_otu_table=self.asv_table,
            meta_data=self.metadata_input,
            decon_method='prevalence',
            prev_control_or_exp_sample_column='Sample_or_ConTrol',
            prev_control_sample_indicator='Control')
        df_output_feature_table = transform(
            output_feature_table,
            from_type=DecontamScoreFormat,
            to_type=pd.DataFrame)
        df_output_feature_table = df_output_feature_table.round(decimals=6)
        exp_table = exp_table.round(decimals=6)
        with tempfile.TemporaryDirectory() as temp_dir_name:
            test_biom_fp = os.path.join(temp_dir_name, 'test_output.tsv')
            expected_biom_fp = os.path.join(temp_dir_name,
                                            'expected_output.tsv')
            df_output_feature_table.to_csv(test_biom_fp, sep="\t")
            exp_table.to_csv(expected_biom_fp, sep="\t")
            with open(test_biom_fp) as fh:
                test_table = biom.Table.from_tsv(fh, None, None, None)
            with open(expected_biom_fp) as th:
                expecter_table = biom.Table.from_tsv(th, None, None, None)

            self.assertEqual(test_table, expecter_table)

    def test_frequency(self):
        exp_table = pd.read_csv(
            self.get_data_path('expected/frequency-score-table.tsv'),
            sep='\t', index_col=0)
        temp_transposed_table = exp_table.transpose()
        temp_transposed_table = temp_transposed_table.dropna()
        exp_table = temp_transposed_table.transpose()
        output_feature_table = decontam_identify(
            asv_or_otu_table=self.asv_table,
            meta_data=self.metadata_input,
            decon_method='frequency',
            freq_concentration_column='quant_reading')
        df_output_feature_table = transform(
            output_feature_table,
            from_type=DecontamScoreFormat,
            to_type=pd.DataFrame)
        df_output_feature_table = df_output_feature_table.round(decimals=6)
        exp_table = exp_table.round(decimals=6)
        with tempfile.TemporaryDirectory() as temp_dir_name:
            test_biom_fp = os.path.join(temp_dir_name, 'test_output.tsv')
            expected_biom_fp = os.path.join(temp_dir_name,
                                            'expected_output.tsv')
            df_output_feature_table.to_csv(test_biom_fp, sep="\t")
            exp_table.to_csv(expected_biom_fp, sep="\t")
            with open(test_biom_fp) as fh:
                test_table = biom.Table.from_tsv(fh, None, None, None)
            with open(expected_biom_fp) as th:
                expecter_table = biom.Table.from_tsv(th, None, None, None)
            self.assertEqual(test_table, expecter_table)

    def test_combined(self):
        exp_table = pd.read_csv(
            self.get_data_path('expected/combined-score-table.tsv'),
            sep='\t', index_col=0)
        temp_transposed_table = exp_table.transpose()
        temp_transposed_table = temp_transposed_table.dropna()
        exp_table = temp_transposed_table.transpose()
        output_feature_table = decontam_identify(
            asv_or_otu_table=self.asv_table,
            meta_data=self.metadata_input,
            decon_method='combined',
            prev_control_or_exp_sample_column='Sample_or_ConTrol',
            prev_control_sample_indicator='Control',
            freq_concentration_column='quant_reading')
        df_output_feature_table = transform(
            output_feature_table,
            from_type=DecontamScoreFormat,
            to_type=pd.DataFrame)
        df_output_feature_table = df_output_feature_table.round(decimals=6)
        exp_table = exp_table.round(decimals=6)
        with tempfile.TemporaryDirectory() as temp_dir_name:
            test_biom_fp = os.path.join(temp_dir_name, 'test_output.tsv')
            expected_biom_fp = os.path.join(temp_dir_name,
                                            'expected_output.tsv')
            df_output_feature_table.to_csv(test_biom_fp, sep="\t")
            exp_table.to_csv(expected_biom_fp, sep="\t")
            with open(test_biom_fp) as fh:
                test_table = biom.Table.from_tsv(fh, None, None, None)
            with open(expected_biom_fp) as th:
                expecter_table = biom.Table.from_tsv(th, None, None, None)
            self.assertEqual(test_table, expecter_table)


class TestRemove(TestPluginBase):
    package = 'q2_quality_control.tests'

    def setUp(self):
        super().setUp()
        table = qiime2.Artifact.load(
            self.get_data_path('expected/decon_default_ASV_table.qza'))
        self.asv_table = table.view(qiime2.Metadata).to_dataframe()
        id_table = qiime2.Artifact.load(
            self.get_data_path('expected/decon_default_score_table.qza'))
        self.identify_table = id_table.view(qiime2.Metadata)

    def test_remove(self):
        exp_table = pd.read_csv(
            self.get_data_path('expected/no-contaminant-asv-table.tsv'),
            sep='\t', index_col=0)
        output_asv_table = decontam_remove(
            asv_or_otu_table=self.asv_table,
            decon_identify_table=self.identify_table,
            threshold=0.1)
        temp_table = output_asv_table.to_dataframe()
        with tempfile.TemporaryDirectory() as temp_dir_name:
            test_biom_fp = os.path.join(temp_dir_name, 'test_output.tsv')
            expected_biom_fp = os.path.join(temp_dir_name,
                                            'expected_output.tsv')
            temp_table.to_csv(test_biom_fp, sep="\t")
            exp_table.to_csv(expected_biom_fp, sep="\t")
            with open(test_biom_fp) as fh:
                test_table = biom.Table.from_tsv(fh, None, None, None)
            with open(expected_biom_fp) as th:
                expecter_table = biom.Table.from_tsv(th, None, None, None)
            self.assertEqual(test_table, expecter_table)


if __name__ == '__main__':
    unittest.main()
