# ----------------------------------------------------------------------------
# Copyright (c) 2017, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import pandas as pd
import qiime2
from warnings import filterwarnings
from qiime2.plugin.testing import TestPluginBase
from q2_types.feature_data import DNAFASTAFormat

from q2_quality_control.quality_control import (exclude_seqs)


filterwarnings("ignore", category=UserWarning)


def _dnafastaformats_to_series(fasta):
    fasta = qiime2.Artifact.import_data("FeatureData[Sequence]", fasta)
    return fasta.view(pd.Series)


class QualityControlTestsBase(TestPluginBase):
    package = 'q2_quality_control.tests'


class SequenceQualityControlBase(QualityControlTestsBase):

    def setUp(self):
        super().setUp()

        def _load_DNAFASTAFormat(reads_fn):
            reads_fp = self.get_data_path(reads_fn)
            return DNAFASTAFormat(reads_fp, mode='r')

        self.query_seqs = _load_DNAFASTAFormat('query-sequences.fasta')
        self.bacterial_ref = _load_DNAFASTAFormat(
            'bacterial-ref-sequences.fasta')
        self.bacterial_exp = _dnafastaformats_to_series(
            _load_DNAFASTAFormat('bacterial-query-sequences.fasta'))
        self.fungal_ref = _load_DNAFASTAFormat('fungal-ref-sequences.fasta')
        self.fungal_exp = _dnafastaformats_to_series(
            _load_DNAFASTAFormat('fungal-query-sequences.fasta'))
        self.query_seqs_with_mismatch = _load_DNAFASTAFormat(
            'query-sequences-with-mismatch.fasta')
        self.query_seqs_short = _load_DNAFASTAFormat(
            'query-sequences-short.fasta')
        self.query_seqs_part_rand = _load_DNAFASTAFormat(
            'query-partially-random.fasta')


class ExcludeSeqsBase(object):

    method = None

    def test_exclude_seqs_bacterial_hit_fungal_miss(self):
        obs, missed = exclude_seqs(
            self.query_seqs, self.bacterial_ref, method=self.method)
        self.assertEqual(
            sorted(obs.index), sorted(self.bacterial_exp.index))
        self.assertEqual(
            sorted(missed.index), sorted(self.fungal_exp.index))

    def test_exclude_seqs_fungal_hit_bacterial_miss(self):
        obs, missed = exclude_seqs(
            self.query_seqs, self.fungal_ref, method=self.method)
        self.assertEqual(sorted(obs.index), sorted(self.fungal_exp.index))
        self.assertEqual(
            sorted(missed.index), sorted(self.bacterial_exp.index))

    def test_exclude_seqs_all_hit(self):
        obs, missed = exclude_seqs(
            self.query_seqs, self.query_seqs, method=self.method)
        self.assertEqual(sorted(obs.index), sorted(
            _dnafastaformats_to_series(self.query_seqs).index))
        self.assertEqual(sorted(missed.index), [])

    def test_exclude_seqs_all_miss(self):
        obs, missed = exclude_seqs(
            self.query_seqs_with_mismatch, self.fungal_ref, method=self.method)
        self.assertEqual(sorted(missed.index), sorted(
            _dnafastaformats_to_series(
                self.query_seqs_with_mismatch).index))
        self.assertEqual(sorted(obs.index), [])

    def test_exclude_seqs_97_perc_identity(self):
        obs, missed = exclude_seqs(
            self.query_seqs_with_mismatch, self.bacterial_ref,
            method=self.method)
        self.assertEqual(
            sorted(obs.index), ['2MISA', '2MISB'])
        self.assertEqual(
            sorted(missed.index), ['10MISA', '8MISA', '8MISB'])

    def test_exclude_seqs_96_perc_identity(self):
        obs, missed = exclude_seqs(
            self.query_seqs_with_mismatch, self.bacterial_ref,
            method=self.method, perc_identity=0.965)
        self.assertEqual(
            sorted(obs.index), ['2MISA', '2MISB', '8MISA', '8MISB'])
        self.assertEqual(
            sorted(missed.index), ['10MISA'])

    def test_exclude_seqs_99_perc_identity(self):
        obs, missed = exclude_seqs(
            self.query_seqs_with_mismatch, self.bacterial_ref,
            method=self.method, perc_identity=0.99)
        self.assertEqual(sorted(missed.index), sorted(
            _dnafastaformats_to_series(
                self.query_seqs_with_mismatch).index))
        self.assertEqual(sorted(obs.index), [])


class BlastTests(ExcludeSeqsBase, SequenceQualityControlBase):
    method = 'blast'


class VsearchTests(ExcludeSeqsBase, SequenceQualityControlBase):
    method = 'vsearch'


class SequenceQualityControlTests(SequenceQualityControlBase):

    def setUp(self):
        super().setUp()

    def test_exclude_seqs_high_evalue_low_perc_query_aligned_permissive(self):
        obs, missed = exclude_seqs(
            self.query_seqs_part_rand, self.bacterial_ref,
            method='blast', perc_identity=0.97, evalue=10000000000000000,
            perc_query_aligned=0.1)
        self.assertEqual(sorted(obs.index), sorted(
            _dnafastaformats_to_series(self.query_seqs_part_rand).index))
        self.assertEqual(sorted(missed.index), [])

    def test_exclude_seqs_blast_low_evalue_discards_weak_matches(self):
        obs, missed = exclude_seqs(
            self.query_seqs_part_rand, self.bacterial_ref,
            method='blast', perc_identity=0.97, evalue=10**-30,
            perc_query_aligned=0.1)
        self.assertEqual(
            sorted(obs.index), ['YAYIMATCH'])
        self.assertEqual(
            sorted(missed.index), ['RAND1', 'RAND2'])

    def test_exclude_seqs_short_seqs_miss_with_default_blast(self):
        obs, missed = exclude_seqs(
            self.query_seqs_short, self.bacterial_ref, method='blast')
        self.assertEqual(sorted(missed.index), sorted(
            _dnafastaformats_to_series(self.query_seqs_short).index))
        self.assertEqual(sorted(obs.index), [])

    def test_exclude_seqs_short_seqs_hit_with_default_vsearch(self):
        obs, missed = exclude_seqs(
            self.query_seqs_short, self.bacterial_ref, method='vsearch')
        self.assertEqual(sorted(obs.index), sorted(
            _dnafastaformats_to_series(self.query_seqs_short).index))
        self.assertEqual(sorted(missed.index), [])

    def test_exclude_seqs_short_seqs_hit_with_blastn_short(self):
        obs, missed = obs, missed = exclude_seqs(
            self.query_seqs_short, self.bacterial_ref,
            method='blastn-short', evalue=10000)
        self.assertEqual(sorted(obs.index), sorted(
            _dnafastaformats_to_series(self.query_seqs_short).index))
        self.assertEqual(sorted(missed.index), [])

    def test_exclude_seqs_short_seqs_miss_with_blastn_short_low_eval(self):
        obs, missed = obs, missed = exclude_seqs(
            self.query_seqs_short, self.bacterial_ref,
            method='blastn-short', perc_identity=0.01, evalue=10**-30)
        self.assertEqual(sorted(missed.index), sorted(
            _dnafastaformats_to_series(self.query_seqs_short).index))
        self.assertEqual(sorted(obs.index), [])
