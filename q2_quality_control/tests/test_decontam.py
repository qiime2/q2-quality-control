# ----------------------------------------------------------------------------
# Copyright (c) 2017-2022, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import unittest
import pandas as pd
import qiime2
import biom
from qiime2.plugin.testing import TestPluginBase
from q2_quality_control.decontam import (decontam_identify,
                                         decontam_remove)
import os
import tempfile


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
        df_output_feature_table = decontam_identify(
            table=self.asv_table,
            metadata=self.metadata_input,
            method='prevalence',
            prev_control_column='Sample_or_Control',
            prev_control_indicator='Control Sample')
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
        df_output_feature_table = decontam_identify(
            table=self.asv_table,
            metadata=self.metadata_input,
            method='frequency',
            freq_concentration_column='quant_reading')
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
        df_output_feature_table = decontam_identify(
            table=self.asv_table,
            metadata=self.metadata_input,
            method='combined',
            prev_control_column='Sample_or_Control',
            prev_control_indicator='Control Sample',
            freq_concentration_column='quant_reading')
        df_output_feature_table = df_output_feature_table.fillna(0)
        df_output_feature_table = df_output_feature_table.round(decimals=6)
        exp_table = exp_table.fillna(0)
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
            table=self.asv_table,
            decontam_scores=self.identify_table,
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


class TestIdentify_mixed_names(TestPluginBase):
    package = 'q2_quality_control.tests'

    def setUp(self):
        super().setUp()
        table = qiime2.Artifact.load(
            self.get_data_path('expected/decon_default_ASV_table.qza'))
        self.asv_table = table.view(qiime2.Metadata).to_dataframe()
        self.metadata_input = qiime2.Metadata.load(
            self.get_data_path(
                'expected/test_metadata_different_ordered_names.tsv'))

    def test_prevalence(self):
        exp_table = pd.read_csv(
            self.get_data_path('expected/prevalence-score-table.tsv'),
            sep='\t', index_col=0)
        temp_transposed_table = exp_table.transpose()
        temp_transposed_table = temp_transposed_table.dropna()
        exp_table = temp_transposed_table.transpose()
        df_output_feature_table = decontam_identify(
            table=self.asv_table,
            metadata=self.metadata_input,
            method='prevalence',
            prev_control_column='Sample_or_Control',
            prev_control_indicator='Control Sample')
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
        df_output_feature_table = decontam_identify(
            table=self.asv_table,
            metadata=self.metadata_input,
            method='frequency',
            freq_concentration_column='quant_reading')
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
        df_output_feature_table = decontam_identify(
            table=self.asv_table,
            metadata=self.metadata_input,
            method='combined',
            prev_control_column='Sample_or_Control',
            prev_control_indicator='Control Sample',
            freq_concentration_column='quant_reading')
        df_output_feature_table = df_output_feature_table.fillna(0)
        df_output_feature_table = df_output_feature_table.round(decimals=6)
        exp_table = exp_table.fillna(0)
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


class TestIdentify_more_names(TestPluginBase):
    package = 'q2_quality_control.tests'

    def setUp(self):
        super().setUp()
        table = qiime2.Artifact.load(
            self.get_data_path('expected/decon_default_ASV_table.qza'))
        self.asv_table = table.view(qiime2.Metadata).to_dataframe()
        self.metadata_input = qiime2.Metadata.load(
            self.get_data_path(
                'expected/test_metadata_more_sample_names.tsv'))

    def test_prevalence(self):
        exp_table = pd.read_csv(
            self.get_data_path('expected/prevalence-score-table.tsv'),
            sep='\t', index_col=0)
        temp_transposed_table = exp_table.transpose()
        temp_transposed_table = temp_transposed_table.dropna()
        exp_table = temp_transposed_table.transpose()
        df_output_feature_table = decontam_identify(
            table=self.asv_table,
            metadata=self.metadata_input,
            method='prevalence',
            prev_control_column='Sample_or_Control',
            prev_control_indicator='Control Sample')
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
        df_output_feature_table = decontam_identify(
            table=self.asv_table,
            metadata=self.metadata_input,
            method='frequency',
            freq_concentration_column='quant_reading')
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
        df_output_feature_table = decontam_identify(
            table=self.asv_table,
            metadata=self.metadata_input,
            method='combined',
            prev_control_column='Sample_or_Control',
            prev_control_indicator='Control Sample',
            freq_concentration_column='quant_reading')
        df_output_feature_table = df_output_feature_table.fillna(0)
        df_output_feature_table = df_output_feature_table.round(decimals=6)
        exp_table = exp_table.fillna(0)
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


if __name__ == '__main__':
    unittest.main()
