# ----------------------------------------------------------------------------
# Copyright (c) 2017, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import numpy as np
import pandas as pd
import qiime2
from warnings import filterwarnings
from qiime2.plugin.testing import TestPluginBase
from q2_types.feature_data import DNAFASTAFormat

from q2_quality_control.quality_control import (
    exclude_seqs, evaluate_taxonomic_composition)
from q2_quality_control._utilities import (
    _evaluate_taxonomic_composition, _collapse_table,
    _drop_nan_zero_and_non_target_rows, _compute_per_level_accuracy,
    compute_taxon_accuracy, _tally_misclassifications,
    _identify_incorrect_classifications, _find_nearest_common_lineage)


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


class UtilitiesTests(QualityControlTestsBase):

    def test_drop_nan_zero_and_non_target_rows(self):
        test_df1 = pd.DataFrame({'a; b': [0., 0., 0.], 'b; c': [1., 0., 0.],
                                 'c; d': [1., np.nan, 1.]})
        test_df2 = pd.DataFrame(
            {'a; b': [1., 1.], 'b; c': [1., 1.], 'c; d': [1., 1.]})
        filtered_df = pd.DataFrame({'b;c': [1.], 'c;d': [1.]})
        new_df = _drop_nan_zero_and_non_target_rows(test_df1, test_df2)
        self.assertTrue(filtered_df.equals(new_df))

    def test_compute_taxon_accuracy(self):
        res = compute_taxon_accuracy(
            pd.Series({'a;b': 1, 'b;c': 1, 'c;d': 1}),
            pd.Series({'a;b': 1, 'b;c': 1, 'c;e': 1, 'd;e': 1}))
        self.assertEqual(res, (0.6666666666666666, 0.5))

    def test_compute_taxon_accuracy_no_matches(self):
        res = compute_taxon_accuracy(
            pd.Series({'a': 1, 'b': 1, 'c': 1}),
            pd.Series({'a;b': 1, 'b;c': 1, 'c;e': 1, 'd;e': 1}))
        self.assertEqual(res, (0., 0.))

    def test_compute_taxon_accuracy_all_match(self):
        res = compute_taxon_accuracy(
            pd.Series({'a;b': 1, 'b;c': 1, 'c;d': 1}),
            pd.Series({'a;b': 1, 'b;c': 1, 'c;d': 1}))
        self.assertEqual(res, (1., 1.))

    def test_collapse_table(self):
        old_table = pd.DataFrame(
            {'a;b;c;d;e': [1, 1, 1], 'a;b;f;g;h': [1, 1, 1],
             'a;b': [1, 1, 1]})
        new_table = _collapse_table(old_table, 2)
        self.assertEquals(set(new_table.columns), set(['a;b']))
        self.assertTrue(np.array_equal(
            new_table.values, np.array([[3], [3], [3]])))
        new_table = _collapse_table(old_table, 3)
        self.assertEquals(set(new_table.columns),
                          set(('a;b;__', 'a;b;c', 'a;b;f')))
        self.assertTrue(np.array_equal(
            new_table.values, np.array([[1, 1, 1], [1, 1, 1], [1, 1, 1]])))

    def test_identify_incorrect_classifications(self):
        res = _identify_incorrect_classifications(
            set(('a', 'c', 'd', 'e', 'f')), set(('a', 'b', 'c', 'd')))
        self.assertEquals(res, (set(('e', 'f')), set(('b'))))

    # matches should never reach this function, so expect underclassification
    def test_tally_misclassifications(self):
        res = _tally_misclassifications(
            set(('a;b;c', 'd;e;f;m', 'g;h', 'j', 'm;n;o')),
            set(('a;b;c', 'd;e;f', 'g;h;i', 'j;k;l')))
        self.assertEquals(
            res, (['d;e;f;m', 'm;n;o'], ['a;b;c', 'g;h', 'j'],
                  [-2, -1, 0, 1, 3]))

    def test_tally_misclassifications_no_match(self):
        res = _tally_misclassifications(
            set(('m;n;o', 'p;q;r', 'q;r;s', 'r;s;t;u')),
            set(('a;b;c', 'd;e;f', 'g;h;i', 'j;k;l')))
        self.assertEquals(res, (['m;n;o', 'p;q;r', 'q;r;s', 'r;s;t;u'], [],
                                [3, 3, 3, 4]))

    # matches should never reach this function, so expect underclassification
    def test_tally_misclassifications_all_match(self):
        res = _tally_misclassifications(
            set(('a;b;c', 'd;e;f', 'g;h;i', 'j;k;l')),
            set(('a;b;c', 'd;e;f', 'g;h;i', 'j;k;l')))
        self.assertEquals(
            res, ([], ['a;b;c', 'd;e;f', 'g;h;i', 'j;k;l'], [0, 0, 0, 0]))


class NCLUtilitiesTests(QualityControlTestsBase):

    def setUp(self):
        super().setUp()
        self.taxa = ['Aa;Bb;Cc;Dd;Ee;Ff;Gg', 'Aa;Bb;Cc;Dd;Ee;Ff;Hh',
                     'Aa;Bb;Cc;Dd;Ee;Ii;Jj', 'Aa;Bb;Cc;Dd;Kk;Ll;Mm']

    def test_find_nearest_common_lineage_match(self):
        res = _find_nearest_common_lineage('Aa;Bb;Cc;Dd;Ee;Ff;Gg', self.taxa)
        self.assertEqual(res, 0)

    def test_find_nearest_common_lineage_no_match(self):
        res = _find_nearest_common_lineage('Ab;Bb;Cb;Db;Eb;b;Gb', self.taxa)
        self.assertEqual(res, 7)

    def test_find_nearest_common_lineage_misclassification(self):
        res = _find_nearest_common_lineage('Aa;Bb;Cc;Dd;Kk;Ll;Mn', self.taxa)
        self.assertEqual(res, 1)

    def test_find_nearest_common_lineage_underclassification(self):
        res = _find_nearest_common_lineage('Aa;Bb;Cc;Dd;Ee;Ff', self.taxa)
        self.assertEqual(res, -1)

    def test_find_nearest_common_lineage_overclassification(self):
        res = _find_nearest_common_lineage('Aa;Bb;Cc;Dd;Ee;Ff;Gg;O', self.taxa)
        self.assertEqual(res, 1)


class EvaluateCompositionTests(QualityControlTestsBase):

    def setUp(self):
        super().setUp()
        self.exp = pd.DataFrame(
            {'k__Aa;p__Ba;c__Ca;o__Da;f__Ea;g__Fa;s__Ga': [0.15, 0.15, 0.15],
             'k__Ab;p__Bb;c__Cb;o__Db;f__Eb;g__Fb;s__Gb': [0.15, 0.15, 0.15],
             'k__Ac;p__Bc;c__Cc;o__Dc;f__Ec;g__Fc;s__Gc': [0.15, 0.15, 0.15],
             'k__Ad;p__Bd;c__Cd;o__Dd;f__Ed;g__Fd;s__Gd': [0.15, 0.15, 0.15],
             'k__Ae;p__Be;c__Ce;o__De;f__Ee;g__Fe;s__Ge': [0.15, 0.15, 0.15],
             'k__Af;p__Bf;c__Cf;o__Df;f__Ef;g__Ff': [0.25, 0.25, 0.25]},
            index=['s1', 's2', 's3'])
        self.obs = pd.DataFrame(
            {'k__Aa;p__Ba;c__Ca;o__Da;f__Ea;g__Fa;s__Ga': [0.10, 0.15, 0.15],
             'k__Ab;p__Bb;c__Cb;o__Db;f__Eb;g__Fb;s__Gb': [0.15, 0.10, 0.10],
             'k__Ac;p__Bc;c__Cc;o__Dc;f__Ec': [0.20, 0.17, 0.15],
             'k__Ad;p__Bd;c__Cd;o__Dd;f__Ed;g__Fd;s__Gd': [0.15, 0.18, 0.20],
             'k__Ae;p__Be;c__Ce;o__De;f__Ee;g__Fe': [0.12, 0.16, 0.15],
             'k__Af;p__Bf;c__Cf;o__Df;f__Ef;g__Ff;s__Gf': [0.20, 0.21, 0.25],
             'k__Ag;p__Bg;c__Cg;o__Dg;f__Eg;g__Fg;s__Gg': [0.08, 0.03, 0.0]},
            index=['s1', 's2', 's3'])

    # test that visualizer runs without fail; internal functions are all tested
    # with various utility tests, this just makes sure the plugin works.
    def test_plugin_evaluate_taxonomic_composition(self):
        evaluate_taxonomic_composition(
            self.temp_dir.name, self.exp, self.obs, 7)

    def test_compute_per_level_accuracy(self):
        res = _compute_per_level_accuracy(self.exp, self.obs, 7)
        self.assertTrue(np.array_equal(res[0].values, exp_res))
        self.assertEqual(res[1], exp_v)

    def test_evaluate_taxonomic_composition(self):
        res = _evaluate_taxonomic_composition(
            self.exp, self.obs, 7, palette='Set1',
            yvals='TAR,TDR,R,Observed / Expected Taxa')
        # results
        self.assertTrue(np.array_equal(res[0].values, exp_res))
        # false negative features
        fn = pd.DataFrame({'s1': [0.15, 0.15, 0.25], 's2': [0.15, 0.15, 0.25],
                           's3': [0.15, 0.15, 0.25]},
                          index=['k__Ac;p__Bc;c__Cc;o__Dc;f__Ec;g__Fc;s__Gc',
                                 'k__Ae;p__Be;c__Ce;o__De;f__Ee;g__Fe;s__Ge',
                                 'k__Af;p__Bf;c__Cf;o__Df;f__Ef;g__Ff'])
        self.assertTrue(res[1].equals(fn))
        # misclassified features
        mc = pd.DataFrame({'s1': [0.20, 0.08], 's2': [0.21, 0.03],
                           's3': [0.25, 0.0]},
                          index=['k__Af;p__Bf;c__Cf;o__Df;f__Ef;g__Ff;s__Gf',
                                 'k__Ag;p__Bg;c__Cg;o__Dg;f__Eg;g__Fg;s__Gg'])
        self.assertTrue(res[2].equals(mc))
        # underclassified features
        uc = pd.DataFrame({'s1': [0.20, 0.12], 's2': [0.17, 0.16],
                           's3': [0.15, 0.15]},
                          index=['k__Ac;p__Bc;c__Cc;o__Dc;f__Ec',
                                 'k__Ae;p__Be;c__Ce;o__De;f__Ee;g__Fe'])
        self.assertTrue(res[3].equals(uc))

    def test_evaluate_taxonomic_composition_invalid_yvals(self):
        with self.assertRaisesRegex(ValueError, "yvals must only"):
            _evaluate_taxonomic_composition(
                self.exp, self.obs, 7, yvals='something_nasty', palette='Set1')


exp_res = np.array(
    [['s1', 1, 7, 1.1666666666666667, 1.0, 0.8571428571428571,
      0.4711111111111111, 0.07555555555555556, 0.7424214437879759,
      0.05597710165934492, 0.19011627371391487],
     ['s2', 1, 7, 1.1666666666666667, 1.0, 0.8571428571428571,
      0.7355555555555554, 0.0377777777777778, 0.8984751466535625,
      0.005968233652573896, 0.16073596169717613],
     ['s3', 1, 7, 1.1666666666666667, 1.0, 1.0, 1.0, 0.0,
      0.9302605094190634, 0.0023753591719758017, 0.1763834207376394],
     ['s1', 2, 7, 1.1666666666666667, 1.0, 0.8571428571428571,
      0.4711111111111111, 0.07555555555555556, 0.7424214437879759,
      0.05597710165934492, 0.19011627371391487],
     ['s2', 2, 7, 1.1666666666666667, 1.0, 0.8571428571428571,
      0.7355555555555554, 0.0377777777777778, 0.8984751466535625,
      0.005968233652573896, 0.16073596169717613],
     ['s3', 2, 7, 1.1666666666666667, 1.0, 1.0, 1.0, 0.0,
      0.9302605094190634, 0.0023753591719758017, 0.1763834207376394],
     ['s1', 3, 7, 1.1666666666666667, 1.0, 0.8571428571428571,
      0.4711111111111111, 0.07555555555555556, 0.7424214437879759,
      0.05597710165934492, 0.19011627371391487],
     ['s2', 3, 7, 1.1666666666666667, 1.0, 0.8571428571428571,
      0.7355555555555554, 0.0377777777777778, 0.8984751466535625,
      0.005968233652573896, 0.16073596169717613],
     ['s3', 3, 7, 1.1666666666666667, 1.0, 1.0, 1.0, 0.0,
      0.9302605094190634, 0.0023753591719758017, 0.1763834207376394],
     ['s1', 4, 7, 1.1666666666666667, 1.0, 0.8571428571428571,
      0.4711111111111111, 0.07555555555555556, 0.7424214437879759,
      0.05597710165934492, 0.19011627371391487],
     ['s2', 4, 7, 1.1666666666666667, 1.0, 0.8571428571428571,
      0.7355555555555554, 0.0377777777777778, 0.8984751466535625,
      0.005968233652573896, 0.16073596169717613],
     ['s3', 4, 7, 1.1666666666666667, 1.0, 1.0, 1.0, 0.0,
      0.9302605094190634, 0.0023753591719758017, 0.1763834207376394],
     ['s1', 5, 7, 1.1666666666666667, 1.0, 0.8571428571428571,
      0.4711111111111111, 0.07555555555555556, 0.7424214437879759,
      0.05597710165934492, 0.19011627371391487],
     ['s2', 5, 7, 1.1666666666666667, 1.0, 0.8571428571428571,
      0.7355555555555554, 0.0377777777777778, 0.8984751466535625,
      0.005968233652573896, 0.16073596169717613],
     ['s3', 5, 7, 1.1666666666666667, 1.0, 1.0, 1.0, 0.0,
      0.9302605094190634, 0.0023753591719758017, 0.1763834207376394],
     ['s1', 6, 7, 1.1666666666666667, 0.8333333333333334,
      0.7142857142857143, 0.06000000000000002, 0.1175,
      0.07644707871564385, 0.8572192090908852, 0.3194787421201396],
     ['s2', 6, 7, 1.1666666666666667, 0.8333333333333334,
      0.7142857142857143, 0.3199999999999999, 0.08500000000000002,
      0.3604847272474663, 0.3803642430099776, 0.3380335289089925],
     ['s3', 6, 7, 1.1666666666666667, 0.8333333333333334,
      0.8333333333333334, 0.5499999999999999, 0.05625000000000001,
      0.5244044240850757, 0.18213394390544005, 0.3645773809037893],
     ['s1', 7, 7, 1.1666666666666667, 0.5, 0.42857142857142855,
      -0.5333333333333334, 0.15333333333333332, -0.6183185276802246,
      0.05671756838000683, 0.2396757068299673],
     ['s2', 7, 7, 1.1666666666666667, 0.5, 0.42857142857142855,
      -0.47333333333333333, 0.14733333333333334, -0.5108045859736967,
      0.13135717451117476, 0.2816518733787826],
     ['s3', 7, 7, 1.1666666666666667, 0.5, 0.5, -0.4333333333333333,
      0.14333333333333334, -0.41957319583913677, 0.2274066375244744,
      0.33145303002252235]], dtype=object)

exp_v = {
    1: {'exp': [
            0.14999999999999999, 0.14999999999999999, 0.14999999999999999,
            0.14999999999999999, 0.14999999999999999, 0.25, 0.0,
            0.14999999999999999, 0.14999999999999999, 0.14999999999999999,
            0.14999999999999999, 0.14999999999999999, 0.25, 0.0,
            0.14999999999999999, 0.14999999999999999, 0.14999999999999999,
            0.14999999999999999, 0.14999999999999999, 0.25, 0.0],
        'obs': [
            0.10000000000000001, 0.14999999999999999, 0.20000000000000001,
            0.14999999999999999, 0.12, 0.20000000000000001,
            0.080000000000000002, 0.14999999999999999, 0.10000000000000001,
            0.17000000000000001, 0.17999999999999999, 0.16,
            0.20999999999999999, 0.029999999999999999, 0.14999999999999999,
            0.10000000000000001, 0.14999999999999999, 0.20000000000000001,
            0.14999999999999999, 0.25, 0.0]},
    2: {'exp': [
            0.14999999999999999, 0.14999999999999999, 0.14999999999999999,
            0.14999999999999999, 0.14999999999999999, 0.25, 0.0,
            0.14999999999999999, 0.14999999999999999, 0.14999999999999999,
            0.14999999999999999, 0.14999999999999999, 0.25, 0.0,
            0.14999999999999999, 0.14999999999999999, 0.14999999999999999,
            0.14999999999999999, 0.14999999999999999, 0.25, 0.0],
        'obs': [
            0.10000000000000001, 0.14999999999999999, 0.20000000000000001,
            0.14999999999999999, 0.12, 0.20000000000000001,
            0.080000000000000002, 0.14999999999999999, 0.10000000000000001,
            0.17000000000000001, 0.17999999999999999, 0.16,
            0.20999999999999999, 0.029999999999999999,
            0.14999999999999999, 0.10000000000000001, 0.14999999999999999,
            0.20000000000000001, 0.14999999999999999, 0.25, 0.0]},
    3: {'exp': [
            0.14999999999999999, 0.14999999999999999, 0.14999999999999999,
            0.14999999999999999, 0.14999999999999999, 0.25, 0.0,
            0.14999999999999999, 0.14999999999999999, 0.14999999999999999,
            0.14999999999999999, 0.14999999999999999, 0.25, 0.0,
            0.14999999999999999, 0.14999999999999999, 0.14999999999999999,
            0.14999999999999999, 0.14999999999999999, 0.25, 0.0],
        'obs': [
            0.10000000000000001, 0.14999999999999999, 0.20000000000000001,
            0.14999999999999999, 0.12, 0.20000000000000001,
            0.080000000000000002, 0.14999999999999999, 0.10000000000000001,
            0.17000000000000001, 0.17999999999999999, 0.16,
            0.20999999999999999, 0.029999999999999999,
            0.14999999999999999, 0.10000000000000001, 0.14999999999999999,
            0.20000000000000001, 0.14999999999999999, 0.25, 0.0]},
    4: {'exp': [
            0.14999999999999999, 0.14999999999999999, 0.14999999999999999,
            0.14999999999999999, 0.14999999999999999, 0.25, 0.0,
            0.14999999999999999, 0.14999999999999999, 0.14999999999999999,
            0.14999999999999999, 0.14999999999999999, 0.25, 0.0,
            0.14999999999999999, 0.14999999999999999, 0.14999999999999999,
            0.14999999999999999, 0.14999999999999999, 0.25, 0.0],
        'obs': [
            0.10000000000000001, 0.14999999999999999, 0.20000000000000001,
            0.14999999999999999, 0.12, 0.20000000000000001,
            0.080000000000000002, 0.14999999999999999, 0.10000000000000001,
            0.17000000000000001, 0.17999999999999999, 0.16,
            0.20999999999999999, 0.029999999999999999,
            0.14999999999999999, 0.10000000000000001, 0.14999999999999999,
            0.20000000000000001, 0.14999999999999999, 0.25, 0.0]},
    5: {'exp': [
            0.14999999999999999, 0.14999999999999999, 0.14999999999999999,
            0.14999999999999999, 0.14999999999999999, 0.25, 0.0,
            0.14999999999999999, 0.14999999999999999, 0.14999999999999999,
            0.14999999999999999, 0.14999999999999999, 0.25, 0.0,
            0.14999999999999999, 0.14999999999999999, 0.14999999999999999,
            0.14999999999999999, 0.14999999999999999, 0.25, 0.0],
        'obs': [
            0.10000000000000001, 0.14999999999999999, 0.20000000000000001,
            0.14999999999999999, 0.12, 0.20000000000000001,
            0.080000000000000002, 0.14999999999999999, 0.10000000000000001,
            0.17000000000000001, 0.17999999999999999, 0.16,
            0.20999999999999999, 0.029999999999999999, 0.14999999999999999,
            0.10000000000000001, 0.14999999999999999,
            0.20000000000000001, 0.14999999999999999, 0.25, 0.0]},
    6: {'exp': [
            0.14999999999999999, 0.14999999999999999, 0.0, 0.14999999999999999,
            0.14999999999999999, 0.14999999999999999, 0.25, 0.0,
            0.14999999999999999, 0.14999999999999999, 0.0, 0.14999999999999999,
            0.14999999999999999, 0.14999999999999999, 0.25, 0.0,
            0.14999999999999999, 0.14999999999999999, 0.0, 0.14999999999999999,
            0.14999999999999999, 0.14999999999999999, 0.25, 0.0],
        'obs': [
            0.10000000000000001, 0.14999999999999999, 0.20000000000000001, 0.0,
            0.14999999999999999, 0.12, 0.20000000000000001,
            0.080000000000000002, 0.14999999999999999, 0.10000000000000001,
            0.17000000000000001, 0.0, 0.17999999999999999, 0.16,
            0.20999999999999999, 0.029999999999999999,
            0.14999999999999999, 0.10000000000000001, 0.14999999999999999, 0.0,
            0.20000000000000001, 0.14999999999999999, 0.25, 0.0]},
    7: {'exp': [
            0.14999999999999999, 0.14999999999999999, 0.0, 0.14999999999999999,
            0.14999999999999999, 0.0, 0.14999999999999999, 0.25, 0.0, 0.0,
            0.14999999999999999, 0.14999999999999999, 0.0, 0.14999999999999999,
            0.14999999999999999, 0.0, 0.14999999999999999, 0.25, 0.0, 0.0,
            0.14999999999999999, 0.14999999999999999, 0.0, 0.14999999999999999,
            0.14999999999999999, 0.0, 0.14999999999999999, 0.25, 0.0, 0.0],
        'obs': [
            0.10000000000000001, 0.14999999999999999, 0.20000000000000001, 0.0,
            0.14999999999999999, 0.12, 0.0, 0.0, 0.20000000000000001,
            0.080000000000000002, 0.14999999999999999, 0.10000000000000001,
            0.17000000000000001, 0.0, 0.17999999999999999, 0.16, 0.0, 0.0,
            0.20999999999999999, 0.029999999999999999, 0.14999999999999999,
            0.10000000000000001, 0.14999999999999999, 0.0, 0.20000000000000001,
            0.14999999999999999, 0.0, 0.0, 0.25, 0.0]}}
