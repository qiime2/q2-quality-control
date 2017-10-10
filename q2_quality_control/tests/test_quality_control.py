# ----------------------------------------------------------------------------
# Copyright (c) 2017, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from qiime2.plugin.testing import TestPluginBase
from q2_quality_control.quality_control import exclude_seqs
from q2_types.feature_data import DNAFASTAFormat


class QualityControlTests(TestPluginBase):
    package = 'q2_quality_control.tests'

    def setUp(self):
        super().setUp()

        def _load_DNAFASTAFormat(reads_fn):
            reads_fp = self.get_data_path(reads_fn)
            return DNAFASTAFormat(reads_fp, mode='r')

        self.query_seqs = _load_DNAFASTAFormat('query-sequences.fasta')
        self.bacterial_ref = _load_DNAFASTAFormat(
            'bacterial-ref-sequences.fasta')
        self.bacterial_exp = _load_DNAFASTAFormat(
            'bacterial-query-sequences.fasta')
        self.fungal_ref = _load_DNAFASTAFormat('fungal-ref-sequences.fasta')
        self.fungal_exp = _load_DNAFASTAFormat('fungal-query-sequences.fasta')


    def test_exclude_seqs_blast(self):
        obs, missed = exclude_seqs(
            self.query_seqs, self.bacterial_ref, method='blast')
        self.assertEqual(sorted(obs.index), sorted(self.bacterial_exp.index))
        self.assertEqual(sorted(missed.index), sorted(self.fungal_exp.index))

    def test_exclude_seqs_vsearch(self):
        obs, missed = exclude_seqs(
            self.query_seqs, self.fungal_ref, method='vsearch')
        self.assertEqual(sorted(obs.index), sorted(self.fungal_exp.index))
        self.assertEqual(
            sorted(missed.index), sorted(self.bacterial_exp.index))
