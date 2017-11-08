# ----------------------------------------------------------------------------
# Copyright (c) 2017, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import numpy as np
import numpy.testing as npt
import pandas as pd
import qiime2
from warnings import filterwarnings
from qiime2.plugin.testing import TestPluginBase
from q2_types.feature_data import DNAFASTAFormat
import pandas.util.testing as pdt

from q2_quality_control.quality_control import (
    exclude_seqs, evaluate_composition, evaluate_seqs)
from q2_quality_control._utilities import (
    _evaluate_composition, _collapse_table, _drop_nans_zeros,
    _compute_per_level_accuracy, compute_taxon_accuracy,
    _tally_misclassifications, _identify_incorrect_classifications,
    _find_nearest_common_lineage, _interpret_metric_selection)
from q2_quality_control._evaluate_seqs import _evaluate_seqs


filterwarnings("ignore", category=UserWarning)


def _dnafastaformats_to_series(fasta):
    fasta = qiime2.Artifact.import_data("FeatureData[Sequence]", fasta)
    return fasta.view(pd.Series)


# test template for EvaluateSeqsTests
def load_evaluate_seqs(query_sequences, reference_sequences, exp_fp):
    results, alignments, g = _evaluate_seqs(
        query_sequences, reference_sequences, show_alignments=False)
    # need to cast to numeric to match dtypes that are interpreted in exp
    # as it is read in by read_csv
    results = results.apply(lambda x: pd.to_numeric(x, errors='ignore'))
    exp = pd.read_csv(exp_fp, sep='\t', index_col=0)
    pdt.assert_frame_equal(results, exp)


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

    def test_drop_nans_zeros(self):
        test_df1 = pd.DataFrame({'a; b': [0., 0., 0.], 'b; c': [1., 0., 0.],
                                 'c; d': [1., np.nan, 1.]})
        filtered_df = pd.DataFrame(
            {'b;c': [1., 0.], 'c;d': [1., 1.]}, index=[0, 2])
        new_df = _drop_nans_zeros(test_df1)
        pdt.assert_frame_equal(filtered_df, new_df)

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
        npt.assert_array_equal(new_table.values, np.array([[3], [3], [3]]))
        new_table = _collapse_table(old_table, 3)
        self.assertEquals(set(new_table.columns),
                          set(('a;b;__', 'a;b;c', 'a;b;f')))
        npt.assert_array_equal(
            new_table.values, np.array([[1, 1, 1], [1, 1, 1], [1, 1, 1]]))

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

    def test_interpret_metric_selection_valid(self):
        yvals = _interpret_metric_selection(
            True, True, False, True, False, False)
        self.assertEquals(yvals, ['TAR', 'TDR', 'r-squared'])

    def test_interpret_metric_selection_invalid(self):
        with self.assertRaisesRegex(ValueError, "At least one metric"):
            _interpret_metric_selection(
                False, False, False, False, False, False)


class EvaluateSeqsTests(SequenceQualityControlBase):

    # check that visualizer runs without fail. Following tests actually test
    # the nuts/bolts for expected values.
    def test_evaluate_seqs_plugin(self):
        evaluate_seqs(output_dir=self.temp_dir.name,
                      query_sequences=self.bacterial_ref,
                      reference_sequences=self.bacterial_ref,
                      show_alignments=True)

    def test_evaluate_seqs_identical(self):
        load_evaluate_seqs(
            self.bacterial_ref, self.bacterial_ref,
            self.get_data_path('test_evaluate_seqs_identical.tsv'))

    def test_evaluate_seqs_mismatches(self):
        load_evaluate_seqs(
            self.query_seqs_with_mismatch, self.bacterial_ref,
            self.get_data_path('test_evaluate_seqs_mismatches.tsv'))

    def test_evaluate_seqs_part_rand(self):
        load_evaluate_seqs(
            self.query_seqs_part_rand, self.bacterial_ref,
            self.get_data_path('test_evaluate_seqs_part_rand.tsv'))


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
        self.exp_results = pd.read_csv(
            self.get_data_path('exp-results.tsv'), sep='\t', index_col=0)
        self.exp = pd.DataFrame(
            {'k__Aa;p__Ba;c__Ca;o__Da;f__Ea;g__Fa;s__Ga': [0.15, 0.15, 0.15],
             'k__Ab;p__Bb;c__Cb;o__Db;f__Eb;g__Fb;s__Gb': [0.15, 0.15, 0.15],
             'k__Ac;p__Bc;c__Cc;o__Dc;f__Ec;g__Fc;s__Gc': [0.15, 0.15, 0.15],
             'k__Ad;p__Bd;c__Cd;o__Dd;f__Ed;g__Fd;s__Gd': [0.15, 0.15, 0.15],
             'k__Ae;p__Be;c__Ce;o__De;f__Ee;g__Fe;s__Ge': [0.15, 0.15, 0.15],
             'k__Af;p__Bf;c__Cf;o__Df;f__Ef;g__Ff': [0.25, 0.25, 0.25]},
            index=['s3', 's1', 's2'])
        self.obs = pd.DataFrame(
            {'k__Aa;p__Ba;c__Ca;o__Da;f__Ea;g__Fa;s__Ga': [0.10, 0.15, 0.15],
             'k__Ab;p__Bb;c__Cb;o__Db;f__Eb;g__Fb;s__Gb': [0.15, 0.10, 0.10],
             'k__Ac;p__Bc;c__Cc;o__Dc;f__Ec': [0.20, 0.17, 0.15],
             'k__Ad;p__Bd;c__Cd;o__Dd;f__Ed;g__Fd;s__Gd': [0.15, 0.18, 0.20],
             'k__Ae;p__Be;c__Ce;o__De;f__Ee;g__Fe': [0.12, 0.16, 0.15],
             'k__Af;p__Bf;c__Cf;o__Df;f__Ef;g__Ff;s__Gf': [0.20, 0.21, 0.25],
             'k__Ag;p__Bg;c__Cg;o__Dg;f__Eg;g__Fg;s__Gg': [0.08, 0.03, 0.0]},
            index=['s1', 's2', 's3'])
        self.false_neg = pd.DataFrame(
            {'s1': [0.15, 0.15, 0.25], 's2': [0.15, 0.15, 0.25],
             's3': [0.15, 0.15, 0.25]},
            index=['k__Ac;p__Bc;c__Cc;o__Dc;f__Ec;g__Fc;s__Gc',
                   'k__Ae;p__Be;c__Ce;o__De;f__Ee;g__Fe;s__Ge',
                   'k__Af;p__Bf;c__Cf;o__Df;f__Ef;g__Ff;__'])
        self.false_neg.index.name = 'Taxon'
        self.misclassified = pd.DataFrame(
            {'s1': [0.20, 0.08], 's2': [0.21, 0.03], 's3': [0.25, 0.0]},
            index=['k__Af;p__Bf;c__Cf;o__Df;f__Ef;g__Ff;s__Gf',
                   'k__Ag;p__Bg;c__Cg;o__Dg;f__Eg;g__Fg;s__Gg'])
        self.misclassified.index.name = 'Taxon'
        self.underclassified = pd.DataFrame(
            {'s1': [0.20, 0.12], 's2': [0.17, 0.16], 's3': [0.15, 0.15]},
            index=['k__Ac;p__Bc;c__Cc;o__Dc;f__Ec;__;__',
                   'k__Ae;p__Be;c__Ce;o__De;f__Ee;g__Fe;__'])
        self.underclassified.index.name = 'Taxon'
        self.exp_one_sample = pd.DataFrame(
            {'k__Aa;p__Ba;c__Ca;o__Da;f__Ea;g__Fa;s__Ga': [0.15],
             'k__Ab;p__Bb;c__Cb;o__Db;f__Eb;g__Fb;s__Gb': [0.15],
             'k__Ac;p__Bc;c__Cc;o__Dc;f__Ec;g__Fc;s__Gc': [0.15],
             'k__Ad;p__Bd;c__Cd;o__Dd;f__Ed;g__Fd;s__Gd': [0.15],
             'k__Ae;p__Be;c__Ce;o__De;f__Ee;g__Fe;s__Ge': [0.15],
             'k__Af;p__Bf;c__Cf;o__Df;f__Ef;g__Ff': [0.25]},
            index=['there_can_only_be_one'])
        self.metadata_one_sample = qiime2.MetadataCategory(pd.DataFrame(
            {'mock_id': ['there_can_only_be_one', 'there_can_only_be_one',
                         'there_can_only_be_one']},
            index=['s3', 's1', 's2'])['mock_id'])

    # test that visualizer runs without fail; internal functions are all tested
    # with various utility tests, this just makes sure the plugin works.
    def test_plugin_evaluate_composition(self):
        evaluate_composition(
            output_dir=self.temp_dir.name, expected_features=self.exp,
            observed_features=self.obs, depth=7)

    def test_compute_per_level_accuracy(self):
        metadata = {_id: _id for _id in self.obs.index}
        res = _compute_per_level_accuracy(self.exp, self.obs, metadata, 7)
        pdt.assert_frame_equal(res[0], self.exp_results)
        self.assertEqual(res[1], exp_v)

    def test_evaluate_composition(self):
        res = _evaluate_composition(
            self.exp, self.obs, depth=7, palette='Set1',
            plot_tar=True, plot_tdr=True, plot_r_value=True,
            plot_r_squared=True, plot_observed_features=True,
            plot_observed_features_ratio=True, metadata=None)
        pdt.assert_frame_equal(res[0], self.exp_results)
        pdt.assert_frame_equal(res[1], self.false_neg)
        pdt.assert_frame_equal(res[2], self.misclassified)
        pdt.assert_frame_equal(res[3], self.underclassified)

    def test_evaluate_composition_metadata_map_to_mock_sample(self):
        res = _evaluate_composition(
            self.exp_one_sample, self.obs, depth=7, palette='Set1',
            plot_tar=True, plot_tdr=True, plot_r_value=True,
            plot_r_squared=True, plot_observed_features=True,
            plot_observed_features_ratio=True,
            metadata=self.metadata_one_sample)
        pdt.assert_frame_equal(res[0], self.exp_results)
        # false_neg should contain only one column header since the map
        # contains one sample (rename to match column name in exp_one_sample)
        false_neg = self.false_neg[['s1']]
        false_neg.columns = ['there_can_only_be_one']
        false_neg.index.name = 'Taxon'
        pdt.assert_frame_equal(res[1], false_neg)
        pdt.assert_frame_equal(res[2], self.misclassified)
        pdt.assert_frame_equal(res[3], self.underclassified)

    def test_evaluate_composition_dont_test_all_levels(self):
        empty_expectations = pd.DataFrame(
            columns=['s1', 's2', 's3']).astype(float)
        empty_expectations.index.name = 'Taxon'
        mc = self.misclassified.loc[[
            'k__Ag;p__Bg;c__Cg;o__Dg;f__Eg;g__Fg;s__Gg']]
        mc.index = ['k__Ag;p__Bg;c__Cg;o__Dg;f__Eg']
        mc.index.name = 'Taxon'
        res = _evaluate_composition(
            self.exp, self.obs, depth=5, palette='Set1', plot_tar=True,
            plot_tdr=True, plot_r_value=True, plot_r_squared=True,
            plot_observed_features=True, plot_observed_features_ratio=True,
            metadata=None)
        pdt.assert_frame_equal(
            res[0], self.exp_results[self.exp_results['level'] < 6])
        pdt.assert_frame_equal(res[1], empty_expectations)
        pdt.assert_frame_equal(res[2], mc)
        pdt.assert_frame_equal(res[3], empty_expectations)

    def test_evaluate_composition_metadata_not_superset(self):
        incomplete_md = qiime2.MetadataCategory(pd.DataFrame(
            {'mock_id': ['there_can_only_be_one']}, index=['s3'])['mock_id'])
        with self.assertRaisesRegex(ValueError, "Missing samples in metadata"):
            _evaluate_composition(
                self.exp_one_sample, self.obs, depth=7, palette='Set1',
                plot_tar=True, plot_tdr=True, plot_r_value=True,
                plot_r_squared=True, plot_observed_features=True,
                plot_observed_features_ratio=True, metadata=incomplete_md)

    def test_evaluate_composition_metadata_values_not_subset(self):
        underrepresented_md = qiime2.MetadataCategory(pd.DataFrame(
            {'mock_id': ['there_can_only_be_one', 'what_is_this?',
                         'i_am_not_really_here']},
            index=['s3', 's1', 's2'])['mock_id'])
        with self.assertRaisesRegex(ValueError, "Missing samples in table"):
            _evaluate_composition(
                self.exp_one_sample, self.obs, depth=7, palette='Set1',
                plot_tar=True, plot_tdr=True, plot_r_value=True,
                plot_r_squared=True, plot_observed_features=True,
                plot_observed_features_ratio=True,
                metadata=underrepresented_md)

    def test_evaluate_composition_depth_too_high(self):
        with self.assertRaisesRegex(ValueError, "8 is larger than the max"):
            _evaluate_composition(
                self.exp, self.obs, depth=8, palette='Set1',
                plot_tar=True, plot_tdr=True, plot_r_value=True,
                plot_r_squared=True, plot_observed_features=True,
                plot_observed_features_ratio=True, metadata=None)


class EvaluateCompositionMockrobiotaDataTests(QualityControlTestsBase):

    def setUp(self):
        super().setUp()
        self.exp_results = pd.read_csv(
            self.get_data_path('mock-3-results.tsv'), sep='\t', index_col=0)
        self.exp = qiime2.Artifact.load(
            self.get_data_path('qc-mock-3-expected.qza')).view(pd.DataFrame)
        self.obs = qiime2.Artifact.load(
            self.get_data_path('qc-mock-3-observed.qza')).view(pd.DataFrame)

        self.false_neg = pd.DataFrame(
            {'HMPMockV1.1.Even1': [0.047619, 0.047619, 0.047619],
             'HMPMockV1.1.Even2': [0.047619, 0.047619, 0.047619],
             'HMPMockV1.2.Staggered1': [0.2143622714, 0.0214362274,
                                        0.0002143626],
             'HMPMockV1.2.Staggered2': [0.2143622714, 0.0214362274,
                                        0.0002143626]},
            index=['k__Bacteria;p__Firmicutes;c__Bacilli;o__Bacillales;'
                   'f__Staphylococcaceae;g__Staphylococcus;s__aureus',
                   'k__Bacteria;p__Firmicutes;c__Bacilli;o__Bacillales;'
                   'f__Staphylococcaceae;g__Staphylococcus;s__epidermidis',
                   'k__Bacteria;p__Thermi;c__Deinococci;o__Deinococcales;'
                   'f__Deinococcaceae;g__Deinococcus;s__'])
        self.false_neg.index.name = 'Taxon'
        self.misclassified = pd.DataFrame(
            {'HMPMockV1.1.Even1': [0.08634],
             'HMPMockV1.1.Even2': [0.0533176566813],
             'HMPMockV1.2.Staggered1': [0.],
             'HMPMockV1.2.Staggered2': [0.]},
            index=['k__Bacteria;p__[Thermi];c__Deinococci;o__Deinococcales;'
                   'f__Deinococcaceae;g__Deinococcus;s__'])
        self.misclassified.index.name = 'Taxon'
        self.underclassified = pd.DataFrame(
            {'HMPMockV1.1.Even1': [0.536876],
             'HMPMockV1.1.Even2': [0.577293],
             'HMPMockV1.2.Staggered1': [0.639295],
             'HMPMockV1.2.Staggered2': [0.666156]},
            index=['k__Bacteria;p__Firmicutes;c__Bacilli;o__Bacillales;'
                   'f__Staphylococcaceae;g__Staphylococcus;__'])
        self.underclassified.index.name = 'Taxon'
        self.metadata = qiime2.MetadataCategory(pd.DataFrame(
            {'mock_id': ['HMPMockV1.1.Even1', 'HMPMockV1.1.Even1',
                         'HMPMockV1.2.Staggered1', 'HMPMockV1.2.Staggered1']},
            index=['HMPMockV1.1.Even1', 'HMPMockV1.1.Even2',
                   'HMPMockV1.2.Staggered1', 'HMPMockV1.2.Staggered2']
            )['mock_id'])

    def test_evaluate_composition_mockrobiota(self):
        res = _evaluate_composition(
            self.exp, self.obs, depth=7, palette='Set1',
            plot_tar=True, plot_tdr=True, plot_r_value=True,
            plot_r_squared=True, plot_observed_features=True,
            plot_observed_features_ratio=True, metadata=None)
        pdt.assert_frame_equal(res[0], self.exp_results)
        pdt.assert_frame_equal(res[1], self.false_neg)
        pdt.assert_frame_equal(res[2], self.misclassified)
        pdt.assert_frame_equal(res[3], self.underclassified)

    def test_evaluate_composition_mockrobiota_metadata_map(self):
        res = _evaluate_composition(
            self.exp, self.obs, depth=7, palette='Set1',
            plot_tar=True, plot_tdr=True, plot_r_value=True,
            plot_r_squared=True, plot_observed_features=True,
            plot_observed_features_ratio=True, metadata=self.metadata)
        pdt.assert_frame_equal(res[0], self.exp_results)
        pdt.assert_frame_equal(res[1], self.false_neg)
        pdt.assert_frame_equal(res[2], self.misclassified)
        pdt.assert_frame_equal(res[3], self.underclassified)


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
