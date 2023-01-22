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
from qiime2.plugin.testing import TestPluginBase
from q2_quality_control._stats import DecontamScoreFormat


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
        _, obs = self.transform_format(DecontamScoreFormat, qiime2.Metadata,
                                       os.path.join('expected',
                                                    'score-table-format.tsv'))

        index = pd.Index(['Seq1', 'Seq2', 'Seq3', 'Seq4', 'Seq5'],
                         name='#OTU ID', dtype=object)
        cols = ['freq', 'prev', 'p.freq', 'p.prev', 'p']
        exp_df = pd.DataFrame(
            [[0.323531341965556, 549, 0.1, 1, 1],
             [0.0988730071632247, 538, 0.2, 1, 1],
             [0.0036748885600188, 160, 0.3, 0.328517354005742,
              0.328517354005742],
             [0.067621013578574, 519, 0.4, 1, 1],
             [0.045234471701338, 354, 0.5, 0.99999999979982,
              0.99999999979982]],
            index=index, columns=cols, dtype=float)
        exp = qiime2.Metadata(exp_df)
        self.assertEqual(exp, obs)

    def test_metadata_to_decontam_table_format(self):
        transformer = self.get_transformer(qiime2.Metadata, DecontamScoreFormat)
        index = pd.Index(['Seq1', 'Seq2', 'Seq3', 'Seq4', 'Seq5'],
                         name='#OTU ID', dtype=object)
        cols = ['freq', 'prev', 'p.freq', 'p.prev', 'p']
        md = pd.DataFrame(
            [[0.323531341965556, 549, 0.1, 1, 1],
             [0.0988730071632247, 538, 0.2, 1, 1],
             [0.0036748885600188, 160, 0.3, 0.328517354005742,
              0.328517354005742],
             [0.067621013578574, 519, 0.4, 1, 1],
             [0.045234471701338, 354, 0.5, 0.99999999979982,
              0.99999999979982]],
            index=index, columns=cols, dtype=float)
        transformer(md)
        self.assertTrue(True)

    def test_decontam_table_format_to_df(self):
        _, obs = self.transform_format(DecontamScoreFormat, pd.DataFrame,
                                       os.path.join('expected',
                                                    'score-table-format.tsv'))

        index = pd.Index(['Seq1', 'Seq2', 'Seq3', 'Seq4', 'Seq5'],
                         name='#OTU ID', dtype=object)
        cols = ['freq', 'prev', 'p.freq', 'p.prev', 'p']
        exp_df = pd.DataFrame(
            [[0.323531341965556, 549, 0.1, 1, 1],
             [0.0988730071632247, 538, 0.2, 1, 1],
             [0.0036748885600188, 160, 0.3, 0.328517354005742,
              0.328517354005742],
             [0.067621013578574, 519, 0.4, 1, 1],
             [0.045234471701338, 354, 0.5, 0.99999999979982,
              0.99999999979982]],
            index=index, columns=cols, dtype=float)
        exp = exp_df
        self.assertEqual(exp, obs)

    def test_df_to_decontam_table_format(self):
        transformer = self.get_transformer(pd.DataFrame, DecontamScoreFormat)
        index = pd.Index(['Seq1', 'Seq2', 'Seq3', 'Seq4', 'Seq5'],
                         name='#OTU ID', dtype=object)
        cols = ['freq', 'prev', 'p.freq', 'p.prev', 'p']
        df = pd.DataFrame(
            [[0.323531341965556, 549, 0.1, 1, 1],
             [0.0988730071632247, 538, 0.2, 1, 1],
             [0.0036748885600188, 160, 0.3, 0.328517354005742,
              0.328517354005742],
             [0.067621013578574, 519, 0.4, 1, 1],
             [0.045234471701338, 354, 0.5, 0.99999999979982,
              0.99999999979982]], index=index, columns=cols, dtype=float)
        # It shouldn't error
        transformer(df)
        self.assertTrue(True)
