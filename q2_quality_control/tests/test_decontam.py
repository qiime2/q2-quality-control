# ----------------------------------------------------------------------------
# Copyright (c) 2017-2022, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import unittest
import pandas as pd
import pandas.testing as pdt
import qiime2
import qiime2.plugin.util
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


def _sort_seqs(seqs):
    return sorted(list(seqs), key=lambda x: x.metadata['id'])


class TestRemove(TestPluginBase):
    package = 'q2_quality_control.tests'

    def setUp(self):
        super().setUp()
        self.input_table = pd.DataFrame(
            [[1, 2, 3, 4, 5], [9, 10, 11, 12, 13]],
            columns=['abc', 'def', 'jkl', 'mno', 'pqr'],
            index=['sample-1', 'sample-2'])
        self.input_seqs = pd.Series(
            ['ACGT', 'TTTT', 'AAAA', 'CCCC', 'GGG'],
            index=['abc', 'def', 'jkl', 'mno', 'pqr'])
        self.input_scores = pd.DataFrame(
            [[13.0, 0.969179],
             [16.0, 0.566067],
             [25.0, 0.019475],
             [10.0, 0.383949],
             [13.0, 0.969179]],
            index=['abc', 'def', 'jkl', 'mno', 'pqr'],
            columns=['prev', 'p'])

    def test_remove(self):
        exp_table = pd.DataFrame(
            [[1, 2, 4, 5], [9, 10, 12, 13]],
            columns=['abc', 'def', 'mno', 'pqr'],
            index=['sample-1', 'sample-2'])
        exp_seqs = pd.Series(['ACGT', 'TTTT', 'CCCC', 'GGG'],
                             index=['abc', 'def', 'mno', 'pqr'])

        obs_table, obs_seqs = decontam_remove(
            table=self.input_table,
            decontam_scores=self.input_scores,
            threshold=0.1,
            rep_seqs=self.input_seqs)

        pdt.assert_series_equal(obs_seqs, exp_seqs)
        pdt.assert_frame_equal(obs_table, exp_table)

    def test_remove_alt_threshold(self):
        exp_table = pd.DataFrame(
            [[1, 2, 5], [9, 10, 13]],
            columns=['abc', 'def', 'pqr'],
            index=['sample-1', 'sample-2'])
        exp_seqs = pd.Series(['ACGT', 'TTTT', 'GGG'],
                             index=['abc', 'def', 'pqr'])

        obs_table, obs_seqs = decontam_remove(
            table=self.input_table,
            decontam_scores=self.input_scores,
            threshold=0.383949,
            rep_seqs=self.input_seqs)

        pdt.assert_series_equal(obs_seqs, exp_seqs)
        pdt.assert_frame_equal(obs_table, exp_table)

    def test_remove_all(self):
        exp_table = pd.DataFrame(
            [], columns=[], index=['sample-1', 'sample-2'])
        exp_seqs = pd.Series([], index=[])

        obs_table, obs_seqs = decontam_remove(
            table=self.input_table,
            decontam_scores=self.input_scores,
            threshold=1.0,
            rep_seqs=self.input_seqs)

        pdt.assert_series_equal(obs_seqs, exp_seqs, check_dtype=False)
        pdt.assert_frame_equal(obs_table, exp_table, check_dtype=False)

    def test_remove_none(self):
        exp_table = pd.DataFrame(
            [[1, 2, 3, 4, 5], [9, 10, 11, 12, 13]],
            columns=['abc', 'def', 'jkl', 'mno', 'pqr'],
            index=['sample-1', 'sample-2'])
        exp_seqs = pd.Series(['ACGT', 'TTTT', 'AAAA', 'CCCC', 'GGG'],
                             index=['abc', 'def', 'jkl', 'mno', 'pqr'])

        obs_table, obs_seqs = decontam_remove(
            table=self.input_table,
            decontam_scores=self.input_scores,
            threshold=0.0,
            rep_seqs=self.input_seqs)

        pdt.assert_series_equal(obs_seqs, exp_seqs)
        pdt.assert_frame_equal(obs_table, exp_table)


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
