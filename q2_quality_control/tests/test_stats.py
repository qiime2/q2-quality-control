# ----------------------------------------------------------------------------
# Copyright (c) 2016-2022, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import os
import pandas as pd
import qiime2
import unittest
from qiime2.plugin.testing import TestPluginBase
from q2_quality_control._stats import DecontamScoreFormat
from qiime2.plugin.util import transform


class TestStatsBoilerplate(TestPluginBase):
    package = 'q2_quality_control.tests'

    def test_decontam_table_format_validate_positive(self):
        filenames = ['score-table-format.tsv', 'prevalence-score-table.tsv']
        filepaths = [self.get_data_path(os.path.join('expected', filename))
                     for filename in filenames]

        for filepath in filepaths:
            format = DecontamScoreFormat(filepath, mode='r')
            # Should pass without error
            format.validate()
            self.assertTrue(True)

    def test_decontam_table_format_to_metadata(self):
        _, obs = self.transform_format(
            DecontamScoreFormat, qiime2.Metadata,
            os.path.join('expected', 'score-table-format.tsv'))

        index = pd.Index(['Seq1', 'Seq2', 'Seq3', 'Seq4', 'Seq5'],
                         name='#OTU ID', dtype=object)
        cols = ['freq', 'prev', 'p.freq', 'p.prev', 'p']
        exp_df = pd.DataFrame(
            [[0.323531342, 549, 0.1, 1, 1],
             [0.098873007, 538, 0.2, 1, 1],
             [0.003674889, 160, 0.3, 0.328517354, 0.328517354],
             [0.067621014, 519, 0.4, 1, 1],
             [0.045234472, 354, 0.5, 1, 1]],
            index=index, columns=cols, dtype=float)
        exp = qiime2.Metadata(exp_df)
        self.assertEqual(exp, obs)

    def test_metadata_to_decontam_table_format(self):
        _, obs_df = self.transform_format(
            DecontamScoreFormat, qiime2.Metadata,
            os.path.join('expected', 'score-table-format.tsv'))
        index = pd.Index(['Seq1', 'Seq2', 'Seq3', 'Seq4', 'Seq5'],
                         name='#OTU ID', dtype=object)
        cols = ['freq', 'prev', 'p.freq', 'p.prev', 'p']
        exp_df = pd.DataFrame(
            [[0.323531342, 549, 0.1, 1, 1],
             [0.098873007, 538, 0.2, 1, 1],
             [0.003674889, 160, 0.3, 0.328517354, 0.328517354],
             [0.067621014, 519, 0.4, 1, 1],
             [0.045234472, 354, 0.5, 1, 1]],
            index=index, columns=cols, dtype=float)
        exp = qiime2.Metadata(exp_df)
        decontam_table = transform(
            exp, from_type=qiime2.Metadata, to_type=DecontamScoreFormat)
        obs_table = transform(
            obs_df, from_type=qiime2.Metadata, to_type=DecontamScoreFormat)
        decontam_table = transform(
            decontam_table, from_type=DecontamScoreFormat,
            to_type=qiime2.Metadata)
        obs_table = transform(
            obs_table, from_type=DecontamScoreFormat, to_type=qiime2.Metadata)
        self.assertEqual(decontam_table, obs_table)

    def test_decontam_table_format_to_df(self):
        _, obs_df = self.transform_format(
            DecontamScoreFormat, pd.DataFrame,
            os.path.join('expected', 'score-table-format.tsv'))

        index = pd.Index(['Seq1', 'Seq2', 'Seq3', 'Seq4', 'Seq5'],
                         name='#OTU ID', dtype=object)
        cols = ['freq', 'prev', 'p.freq', 'p.prev', 'p']
        exp_df = pd.DataFrame(
            [[0.323531342, 549, 0.1, 1, 1],
             [0.098873007, 538, 0.2, 1, 1],
             [0.003674889, 160, 0.3, 0.328517354, 0.328517354],
             [0.067621014, 519, 0.4, 1, 1],
             [0.045234472, 354, 0.5, 1, 1]],
            index=index, columns=cols, dtype=float)
        exp = qiime2.Metadata(exp_df)
        obs = qiime2.Metadata(obs_df)

        self.assertEqual(exp, obs)

    def test_df_to_decontam_table_format(self):
        _, obs_df = self.transform_format(
            DecontamScoreFormat, pd.DataFrame,
            os.path.join('expected', 'score-table-format.tsv'))
        index = pd.Index(['Seq1', 'Seq2', 'Seq3', 'Seq4', 'Seq5'],
                         name='#OTU ID', dtype=object)
        cols = ['freq', 'prev', 'p.freq', 'p.prev', 'p']
        exp_df = pd.DataFrame(
            [[0.323531342, 549, 0.1, 1, 1],
             [0.098873007, 538, 0.2, 1, 1],
             [0.003674889, 160, 0.3, 0.328517354, 0.328517354],
             [0.067621014, 519, 0.4, 1, 1],
             [0.045234472, 354, 0.5, 1, 1]],
            index=index, columns=cols, dtype=float)
        exp_table = transform(
            exp_df, from_type=pd.DataFrame, to_type=DecontamScoreFormat)
        obs_table = transform(
            obs_df, from_type=pd.DataFrame, to_type=DecontamScoreFormat)
        decontam_table = transform(
            exp_table, from_type=DecontamScoreFormat, to_type=pd.DataFrame)
        obs_table = transform(
            obs_table, from_type=DecontamScoreFormat, to_type=pd.DataFrame)
        exp_test = qiime2.Metadata(decontam_table)
        obs = qiime2.Metadata(obs_table)

        self.assertEqual(exp_test, obs)


if __name__ == '__main__':
    unittest.main()
    