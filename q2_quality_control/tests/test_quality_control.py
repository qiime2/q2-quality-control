# ----------------------------------------------------------------------------
# Copyright (c) 2017-2021, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import numpy as np
import numpy.testing as npt
import pandas as pd
import qiime2
import biom
from warnings import filterwarnings
from qiime2.plugin.testing import TestPluginBase
from q2_types.feature_data import DNAFASTAFormat
import pandas.util.testing as pdt

from q2_quality_control.quality_control import (
    exclude_seqs, evaluate_composition, evaluate_seqs, evaluate_taxonomy)
from q2_quality_control._utilities import (
    _evaluate_composition, _collapse_table, _drop_nans_zeros,
    _compute_per_level_accuracy, compute_taxon_accuracy,
    _tally_misclassifications, _identify_incorrect_classifications,
    _find_nearest_common_lineage, _interpret_metric_selection,
    _match_samples_by_index, _validate_metadata_and_exp_table)
from q2_quality_control._evaluate_seqs import _evaluate_seqs
from q2_quality_control._evaluate_taxonomy import (
    _evaluate_taxonomy, _extract_taxa_names, _index_is_subset,
    _validate_indices_and_set_joining_mode)


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

    def setUp(self):
        super().setUp()


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
        self.query_seqs_left_justify = _load_DNAFASTAFormat(
            'query-sequences-left-justified.fasta')


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

    def test_exclude_seqs_left_justify_value_error(self):
        with self.assertRaisesRegex(ValueError, "left_justify is not "
                                                "compatible with "
                                                "method='blast'"):
            exclude_seqs(
                self.query_seqs_left_justify, self.bacterial_ref,
                method=self.method, left_justify=True)


class VsearchTests(ExcludeSeqsBase, SequenceQualityControlBase):
    method = 'vsearch'

    def test_exclude_seqs_left_justify_hits(self):
        obs, missed = exclude_seqs(
            self.query_seqs_left_justify, self.bacterial_ref,
            method=self.method,
            left_justify=True)
        self.assertCountEqual(
            sorted(obs.index),
            ['1111886-leftjustmatch', '1111882-leftjustwithindel',
             '1111879-leftjustmatch']
        )

    def test_exclude_seqs_left_justify_missed(self):
        obs, missed = exclude_seqs(
            self.query_seqs_left_justify, self.bacterial_ref,
            method=self.method, left_justify=True)
        self.assertEqual(sorted(missed.index), ['1111883-noleftjustmatch'])


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
        self.assertEqual(res, (0.5, 0.6666666666666666))

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
        self.assertEqual(set(new_table.columns), set(['a;b']))
        npt.assert_array_equal(new_table.values, np.array([[3], [3], [3]]))
        new_table = _collapse_table(old_table, 3)
        self.assertEqual(set(new_table.columns),
                         set(('a;b;__', 'a;b;c', 'a;b;f')))
        npt.assert_array_equal(
            new_table.values, np.array([[1, 1, 1], [1, 1, 1], [1, 1, 1]]))

    def test_identify_incorrect_classifications(self):
        res = _identify_incorrect_classifications(
            set(('a', 'c', 'd', 'e', 'f')), set(('a', 'b', 'c', 'd')))
        self.assertEqual(res, (set(('e', 'f')), set(('b'))))

    # matches should never reach this function, so expect underclassification
    def test_tally_misclassifications(self):
        res = _tally_misclassifications(
            set(('a;b;c', 'd;e;f;m', 'g;h', 'j', 'm;n;o')),
            set(('a;b;c', 'd;e;f', 'g;h;i', 'j;k;l')))
        self.assertEqual(
            res, (['d;e;f;m', 'm;n;o'], ['a;b;c', 'g;h', 'j'],
                  [-2, -1, 0, 1, 3]))

    def test_tally_misclassifications_no_match(self):
        res = _tally_misclassifications(
            set(('m;n;o', 'p;q;r', 'q;r;s', 'r;s;t;u')),
            set(('a;b;c', 'd;e;f', 'g;h;i', 'j;k;l')))
        self.assertEqual(res, (['m;n;o', 'p;q;r', 'q;r;s', 'r;s;t;u'], [],
                               [3, 3, 3, 4]))

    # matches should never reach this function, so expect underclassification
    def test_tally_misclassifications_all_match(self):
        res = _tally_misclassifications(
            set(('a;b;c', 'd;e;f', 'g;h;i', 'j;k;l')),
            set(('a;b;c', 'd;e;f', 'g;h;i', 'j;k;l')))
        self.assertEqual(
            res, ([], ['a;b;c', 'd;e;f', 'g;h;i', 'j;k;l'], [0, 0, 0, 0]))

    def test_interpret_metric_selection_valid(self):
        yvals = _interpret_metric_selection(
            True, True, False, True, False, False, False, False)
        self.assertEqual(yvals, ['TAR', 'TDR', 'r-squared'])

    def test_interpret_metric_selection_invalid(self):
        with self.assertRaisesRegex(ValueError, "At least one metric"):
            _interpret_metric_selection(
                False, False, False, False, False, False, False, False)

    def test_match_samples_by_index(self):
        df_a = pd.DataFrame({'a': (1., 2., 3.)}, index=['1', '2', '3'])
        df_b = pd.DataFrame({'a': (2., 3., 7.)}, index=['2', '3', '7'])
        df_c = pd.DataFrame({'a': (2., 3.)}, index=['2', '3'])
        df_a, df_b = _match_samples_by_index(df_a, df_b)
        pdt.assert_frame_equal(df_a.dropna(), df_c)
        pdt.assert_frame_equal(df_b.dropna(), df_c)

    def test_validate_metadata_and_exp_table_metadata_not_superset(self):
        incomplete_md = pd.Series(['there_can_only_be_one'], name='mock_id',
                                  index=pd.Index(['s3'], name='id'))
        with self.assertRaisesRegex(ValueError, "Missing samples in metadata"):
            _validate_metadata_and_exp_table(
                incomplete_md, exp_one_sample, obs)

    def test_validate_metadata_and_exp_table_metadata_values_not_subset(self):
        underrepresented_md = pd.Series(
            ['there_can_only_be_one', 'what_is_this?', 'i_am_not_really_here'],
            name='mock_id', index=pd.Index(['s3', 's1', 's2'], name='id'))
        with self.assertRaisesRegex(ValueError, "Missing samples in table"):
            _validate_metadata_and_exp_table(
                underrepresented_md, exp_one_sample, obs)

    # test for issue-25: metadata is superset, but values are subset of exp
    # should output new metadata that fits the table
    def test_validate_metadata_and_exp_table_md_superset_vals_subset(self):
        overabundant_md = pd.Series(
            ['s1', 's1', 's1', 's2', 's2'], name='mock_id',
            index=pd.Index(['s3', 's1', 's2', 'fake1', 'fake2'], name='id'))
        corrected_md = pd.Series(
            ['s1', 's1', 's1'], name='mock_id',
            index=pd.Index(['s1', 's2', 's3']))
        new_md, junk = _validate_metadata_and_exp_table(
            overabundant_md, exp, obs)
        pdt.assert_series_equal(new_md, corrected_md)


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
        self.exp = exp
        self.obs = obs
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
        self.exp_one_sample = exp_one_sample
        self.metadata_one_sample = qiime2.CategoricalMetadataColumn(
            pd.Series(['there_can_only_be_one',
                       'there_can_only_be_one',
                       'there_can_only_be_one'], name='mock_id',
                      index=pd.Index(['s3', 's1', 's2'], name='id')))

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

    # confirm accurate behavior when exp/obs vector lengths are equal to
    # one. See https://github.com/qiime2/q2-quality-control/issues/44
    def test_compute_per_level_accuracy_vector_length_one(self):
        res = _compute_per_level_accuracy(
            exp_one_kingdom, obs_one_kingdom, None, depth=1)
        # assert that all r2 values are nan (if r2 is nan, r will be also)
        self.assertTrue(np.isnan(res[0]['r-squared']).all())

    def test_evaluate_composition(self):
        res = _evaluate_composition(
            self.exp, self.obs, depth=7, palette='Set1',
            plot_tar=True, plot_tdr=True, plot_r_value=True,
            plot_r_squared=True, plot_bray_curtis=False,
            plot_jaccard=False, plot_observed_features=True,
            plot_observed_features_ratio=True, metadata=None)
        pdt.assert_frame_equal(res[0], self.exp_results)
        pdt.assert_frame_equal(res[1], self.false_neg)
        pdt.assert_frame_equal(res[2], self.misclassified)
        pdt.assert_frame_equal(res[3], self.underclassified)

    def test_evaluate_composition_metadata_map_to_mock_sample(self):
        res = _evaluate_composition(
            self.exp_one_sample, self.obs, depth=7, palette='Set1',
            plot_tar=True, plot_tdr=True, plot_r_value=True,
            plot_r_squared=True, plot_bray_curtis=False,
            plot_jaccard=False, plot_observed_features=True,
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
            plot_observed_features=True, plot_bray_curtis=False,
            plot_jaccard=False, plot_observed_features_ratio=True,
            metadata=None)
        pdt.assert_frame_equal(
            res[0], self.exp_results[self.exp_results['level'] < 6])
        pdt.assert_frame_equal(res[1], empty_expectations)
        pdt.assert_frame_equal(res[2], mc)
        pdt.assert_frame_equal(res[3], empty_expectations)

    def test_evaluate_composition_metadata_not_superset(self):
        incomplete_md = qiime2.CategoricalMetadataColumn(
                    pd.Series(['there_can_only_be_one'], name='mock_id',
                              index=pd.Index(['s3'], name='id')))
        with self.assertRaisesRegex(ValueError, "Missing samples in metadata"):
            _evaluate_composition(
                self.exp_one_sample, self.obs, depth=7, palette='Set1',
                plot_tar=True, plot_tdr=True, plot_r_value=True,
                plot_r_squared=True, plot_bray_curtis=False,
                plot_jaccard=False, plot_observed_features=True,
                plot_observed_features_ratio=True, metadata=incomplete_md)

    def test_evaluate_composition_metadata_values_not_subset(self):
        underrepresented_md = qiime2.CategoricalMetadataColumn(
            pd.Series(['there_can_only_be_one', 'what_is_this?',
                       'i_am_not_really_here'],
                      name='mock_id', index=pd.Index(['s3', 's1', 's2'],
                                                     name='id')))
        with self.assertRaisesRegex(ValueError, "Missing samples in table"):
            _evaluate_composition(
                self.exp_one_sample, self.obs, depth=7, palette='Set1',
                plot_tar=True, plot_tdr=True, plot_r_value=True,
                plot_r_squared=True, plot_bray_curtis=False,
                plot_jaccard=False, plot_observed_features=True,
                plot_observed_features_ratio=True,
                metadata=underrepresented_md)

    def test_evaluate_composition_depth_too_high(self):
        with self.assertRaisesRegex(ValueError, "8 is larger than the max"):
            _evaluate_composition(
                self.exp, self.obs, depth=8, palette='Set1',
                plot_tar=True, plot_tdr=True, plot_r_value=True,
                plot_r_squared=True, plot_bray_curtis=False,
                plot_jaccard=False, plot_observed_features=True,
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
        self.metadata = qiime2.CategoricalMetadataColumn(
            pd.Series(['HMPMockV1.1.Even1',
                       'HMPMockV1.1.Even1',
                       'HMPMockV1.2.Staggered1',
                       'HMPMockV1.2.Staggered1'], name='mock_id',
                      index=pd.Index(['HMPMockV1.1.Even1',
                                      'HMPMockV1.1.Even2',
                                      'HMPMockV1.2.Staggered1',
                                      'HMPMockV1.2.Staggered2'], name='id')))

    def test_evaluate_composition_mockrobiota(self):
        res = _evaluate_composition(
            self.exp, self.obs, depth=7, palette='Set1',
            plot_tar=True, plot_tdr=True, plot_r_value=True,
            plot_r_squared=True, plot_bray_curtis=True,
            plot_jaccard=True, plot_observed_features=True,
            plot_observed_features_ratio=True, metadata=None)
        pdt.assert_frame_equal(res[0], self.exp_results)
        pdt.assert_frame_equal(res[1], self.false_neg)
        pdt.assert_frame_equal(res[2], self.misclassified)
        pdt.assert_frame_equal(res[3], self.underclassified)

    def test_evaluate_composition_mockrobiota_metadata_map(self):
        res = _evaluate_composition(
            self.exp, self.obs, depth=7, palette='Set1',
            plot_tar=True, plot_tdr=True, plot_r_value=True,
            plot_r_squared=True, plot_bray_curtis=False,
            plot_jaccard=False, plot_observed_features=True,
            plot_observed_features_ratio=True, metadata=self.metadata)
        false_neg = self.false_neg[['HMPMockV1.1.Even1',
                                    'HMPMockV1.2.Staggered1']]
        pdt.assert_frame_equal(res[0], self.exp_results)
        pdt.assert_frame_equal(res[1], false_neg)
        pdt.assert_frame_equal(res[2], self.misclassified)
        pdt.assert_frame_equal(res[3], self.underclassified)


class EvaluateTaxonomyTests(QualityControlTestsBase):

    def setUp(self):
        super().setUp()
        self.exp_taxa = pd.read_csv(
            self.get_data_path('mock-3-expected-taxonomy.tsv'), sep='\t',
            index_col=0)
        self.obs_taxa = pd.read_csv(
            self.get_data_path('mock-3-observed-taxonomy.tsv'), sep='\t',
            index_col=0)
        self.obs_table = biom.load_table(
            self.get_data_path('mock-3-obs-table.biom'))
        self.prf_res_unw = pd.DataFrame.from_dict({
            'level': {0: 1, 1: 2, 2: 3, 3: 4, 4: 5, 5: 6, 6: 7},
            'Precision': {0: 0.96, 1: 0.96, 2: 0.96, 3: 0.96, 4: 0.96,
                          5: 0.9583333333333334, 6: 0.7368421052631579},
            'Recall': {0: 0.96, 1: 0.96, 2: 0.96, 3: 0.96, 4: 0.96,
                       5: 0.92, 6: 0.56},
            'F-measure': {0: 0.96, 1: 0.96, 2: 0.96, 3: 0.96, 4: 0.96,
                          5: 0.9387755102040817, 6: 0.6363636363636364}})[[
                             'level', 'Precision', 'Recall', 'F-measure']]

    # just check that visualizer works. Other tests check nuts/bolts
    def test_evaluate_taxonomy(self):
        evaluate_taxonomy(
            self.temp_dir.name, self.exp_taxa, self.obs_taxa, depth=7)

    def test_evaluate_taxonomy_basic(self):
        prf = _evaluate_taxonomy(
            self.exp_taxa, self.obs_taxa, require_exp_ids=False,
            require_obs_ids=False, feature_table=None, sample_id=None,
            level_range=range(0, 7))
        pdt.assert_frame_equal(prf, self.prf_res_unw)

    def test_evaluate_taxonomy_sample_id_no_table(self):
        prf = _evaluate_taxonomy(
            self.exp_taxa, self.obs_taxa, require_exp_ids=False,
            require_obs_ids=False, feature_table=None,
            sample_id='HMPMockV1.1.Even1', level_range=range(0, 7))
        prf_res_unw = self.prf_res_unw
        prf_res_unw['SampleID'] = 'HMPMockV1.1.Even1'
        pdt.assert_frame_equal(prf, prf_res_unw)

    def test_evaluate_taxonomy_sample_id_not_found_FAIL(self):
        with self.assertRaisesRegex(
                ValueError, 'Sample id not found in feature table: Peanut'):
            _evaluate_taxonomy(
                self.exp_taxa, self.obs_taxa, require_exp_ids=False,
                require_obs_ids=False, feature_table=self.obs_table,
                sample_id='Peanut')

    def test_evaluate_taxonomy_plus_table_no_sample_id(self):
        prf = _evaluate_taxonomy(
            self.exp_taxa, self.obs_taxa, require_exp_ids=False,
            require_obs_ids=False, feature_table=self.obs_table,
            sample_id=None, level_range=range(6, 7))
        pdt.assert_frame_equal(prf, pd.DataFrame.from_dict(
            {'level': {0: 7},
             'Precision': {0: 0.6377837950426994},
             'Recall': {0: 0.19487462344605203},
             'F-measure': {0: 0.2985326855267221}})[[
                'level', 'Precision', 'Recall', 'F-measure']])

    def test_evaluate_taxonomy_sample_id_plus_table(self):
        prf = _evaluate_taxonomy(
            self.exp_taxa, self.obs_taxa, require_exp_ids=False,
            require_obs_ids=False, feature_table=self.obs_table,
            sample_id='HMPMockV1.1.Even1', level_range=range(6, 7))
        pdt.assert_frame_equal(prf, pd.DataFrame.from_dict(
            {'level': {0: 7},
             'Precision': {0: 0.5523983315954119},
             'Recall': {0: 0.19757575757575757},
             'F-measure': {0: 0.2910514387748094},
             'SampleID': {0: 'HMPMockV1.1.Even1'}})[[
                'level', 'Precision', 'Recall', 'F-measure', 'SampleID']])

    def test_evaluate_taxonomy_plus_table_missing_features_in_table_FAIL(self):
        obs_table = self.obs_table.filter(
            ['47ad35356a9bfec68416d32e4f039021',
             'c4269a6e9bd66eca53e710c9f9d9ad4f',
             '7a8d29c59b803baaed9cc1f04ce0dc33',
             '86adb6193435090318cf24df07770d07'], axis='observation')

        with self.assertRaisesRegex(
                ValueError, "Feature ids not found in feature table"):
            _evaluate_taxonomy(
                self.exp_taxa, self.obs_taxa, require_exp_ids=False,
                require_obs_ids=False, feature_table=obs_table,
                sample_id=None, level_range=range(6, 7))

    def _test_extract_taxa_names_species(self):
        names = _extract_taxa_names(self.self.exp_taxa.iloc['Taxon'])
        self.assertEqual(names, ['s__', 's__', 's__sphaeroides', 's__acnes',
                                 's__aureus', '', 's__coli', 's__pylori',
                                 's__', 's__', 's__', 's__', 's__cereus',
                                 's__', 's__', 's__', 's__', 's__', 's__',
                                 's__', 's__', 's__monocytogenes',
                                 's__aeruginosa', 's__', 's__agalactiae'])

    def _test_extract_taxa_names_phylum(self):
        names = _extract_taxa_names(
            self.self.exp_taxa.iloc['Taxon'], level=slice(0, 1))
        self.assertEqual(names, ['k__Bacteria', 'k__Bacteria', 'k__Bacteria',
                                 'k__Bacteria', 'k__Bacteria', 'other',
                                 'k__Bacteria', 'k__Bacteria', 'k__Archaea',
                                 'k__Bacteria', 'k__Archaea', 'k__Bacteria',
                                 'k__Bacteria', 'k__Bacteria', 'k__Bacteria',
                                 'k__Bacteria', 'k__Bacteria', 'k__Bacteria',
                                 'k__Archaea', 'k__Bacteria', 'k__Bacteria',
                                 'k__Bacteria', 'k__Bacteria', 'k__Bacteria',
                                 'k__Bacteria'])

    def _test_extract_taxa_names_all_levels(self):
        names = _extract_taxa_names(
            self.self.exp_taxa['Taxon'], level=slice(0, 7))
        self.assertEqual(names, list(self.self.exp_taxa.iloc['Taxon']))

    def test_index_is_subset_equal_series(self):
        # if series1 is superset of or == series2, no problem
        _index_is_subset(self.exp_taxa, self.obs_taxa, 'observed')

    def test_index_is_subset_TRUE(self):
        _index_is_subset(self.exp_taxa, self.obs_taxa[:-5], 'observed')

    def test_index_is_subset_FALSE(self):
        # if series1 is subset of series2, raise error
        with self.assertRaisesRegex(
                ValueError, 'ids not found in observed ids'):
            _index_is_subset(self.exp_taxa[:-5], self.obs_taxa, 'observed')

    def test_validate_indices_and_set_joining_mode_equal(self):
        # series1 == series2 indices, no prob bob
        join_how = _validate_indices_and_set_joining_mode(
            self.exp_taxa, self.obs_taxa, require_exp_ids=True,
            require_obs_ids=True)
        self.assertEqual(join_how, 'inner')

    def test_validate_indices_and_set_joining_mode_inner(self):
        # any time require_exp_id == require_obs_ids, join_how = 'inner'
        join_how = _validate_indices_and_set_joining_mode(
            self.exp_taxa, self.obs_taxa, require_exp_ids=False,
            require_obs_ids=False)
        self.assertEqual(join_how, 'inner')

    def test_validate_indices_and_set_joining_mode_right(self):
        # exp < obs indices but require_obs_ids=False, pass
        join_how = _validate_indices_and_set_joining_mode(
            self.exp_taxa[-5:], self.obs_taxa, require_exp_ids=True,
            require_obs_ids=False)
        self.assertEqual(join_how, 'left')

    def test_validate_indices_and_set_joining_mode_left(self):
        # exp > obs indices but require_exp_ids=False, pass
        join_how = _validate_indices_and_set_joining_mode(
            self.exp_taxa, self.obs_taxa[-5:], require_exp_ids=False,
            require_obs_ids=True)
        self.assertEqual(join_how, 'right')

    def test_validate_indices_and_set_joining_mode_require_obs_id_FAIL(self):
        # series1 < series2 indices but require_obs_id=True, FAIL
        with self.assertRaisesRegex(
                ValueError, 'ids not found in expected ids'):
            _validate_indices_and_set_joining_mode(
                self.exp_taxa[-5:], self.obs_taxa, require_exp_ids=False,
                require_obs_ids=True)

    def test_validate_indices_and_set_joining_mode_require_exp_id_FAIL(self):
        # series1 > series2 indices but require_exp_ids=True, FAIL
        with self.assertRaisesRegex(
                ValueError, 'ids not found in observed ids'):
            _validate_indices_and_set_joining_mode(
                self.exp_taxa, self.obs_taxa[-5:], require_exp_ids=True,
                require_obs_ids=False)


exp = pd.DataFrame(
    {'k__Aa;p__Ba;c__Ca;o__Da;f__Ea;g__Fa;s__Ga': [0.15, 0.15, 0.15],
     'k__Ab;p__Bb;c__Cb;o__Db;f__Eb;g__Fb;s__Gb': [0.15, 0.15, 0.15],
     'k__Ac;p__Bc;c__Cc;o__Dc;f__Ec;g__Fc;s__Gc': [0.15, 0.15, 0.15],
     'k__Ad;p__Bd;c__Cd;o__Dd;f__Ed;g__Fd;s__Gd': [0.15, 0.15, 0.15],
     'k__Ae;p__Be;c__Ce;o__De;f__Ee;g__Fe;s__Ge': [0.15, 0.15, 0.15],
     'k__Af;p__Bf;c__Cf;o__Df;f__Ef;g__Ff': [0.25, 0.25, 0.25]},
    index=['s1', 's2', 's3'])

exp_one_sample = pd.DataFrame(
    {'k__Aa;p__Ba;c__Ca;o__Da;f__Ea;g__Fa;s__Ga': [0.15],
     'k__Ab;p__Bb;c__Cb;o__Db;f__Eb;g__Fb;s__Gb': [0.15],
     'k__Ac;p__Bc;c__Cc;o__Dc;f__Ec;g__Fc;s__Gc': [0.15],
     'k__Ad;p__Bd;c__Cd;o__Dd;f__Ed;g__Fd;s__Gd': [0.15],
     'k__Ae;p__Be;c__Ce;o__De;f__Ee;g__Fe;s__Ge': [0.15],
     'k__Af;p__Bf;c__Cf;o__Df;f__Ef;g__Ff': [0.25]},
    index=['there_can_only_be_one'])

obs = pd.DataFrame(
    {'k__Aa;p__Ba;c__Ca;o__Da;f__Ea;g__Fa;s__Ga': [0.10, 0.15, 0.15],
     'k__Ab;p__Bb;c__Cb;o__Db;f__Eb;g__Fb;s__Gb': [0.15, 0.10, 0.10],
     'k__Ac;p__Bc;c__Cc;o__Dc;f__Ec': [0.20, 0.17, 0.15],
     'k__Ad;p__Bd;c__Cd;o__Dd;f__Ed;g__Fd;s__Gd': [0.15, 0.18, 0.20],
     'k__Ae;p__Be;c__Ce;o__De;f__Ee;g__Fe': [0.12, 0.16, 0.15],
     'k__Af;p__Bf;c__Cf;o__Df;f__Ef;g__Ff;s__Gf': [0.20, 0.21, 0.25],
     'k__Ag;p__Bg;c__Cg;o__Dg;f__Eg;g__Fg;s__Gg': [0.08, 0.03, 0.0]},
    index=['s1', 's2', 's3'])

exp_one_kingdom = pd.DataFrame(
    {'k__A;p__B;c__C;o__Da;f__Ea;g__Fa;s__Ga': [0.15, 0.15, 0.15],
     'k__A;p__B;c__C;o__Db;f__Eb;g__Fb;s__Gb': [0.15, 0.15, 0.15],
     'k__A;p__B;c__Cc;o__Dc;f__Ec;g__Fc;s__Gc': [0.15, 0.15, 0.15],
     'k__A;p__Bd;c__Cd;o__Dd;f__Ed;g__Fd;s__Gd': [0.15, 0.15, 0.15],
     'k__A;p__Be;c__Ce;o__De;f__Ee;g__Fe;s__Ge': [0.15, 0.15, 0.15],
     'k__A;p__Bf;c__Cf;o__Df;f__Ef;g__Ff': [0.25, 0.25, 0.25]},
    index=['s1', 's2', 's3'])

obs_one_kingdom = pd.DataFrame(
    {'k__A;p__B;c__C;o__Da;f__Ea;g__Fa;s__Ga': [0.10, 0.15, 0.15],
     'k__A;p__B;c__C;o__Db;f__Eb;g__Fb;s__Gb': [0.15, 0.10, 0.10],
     'k__A;p__B;c__Cc;o__Dc;f__Ec': [0.20, 0.17, 0.15],
     'k__A;p__Bd;c__Cd;o__Dd;f__Ed;g__Fd;s__Gd': [0.15, 0.18, 0.20],
     'k__A;p__Be;c__Ce;o__De;f__Ee;g__Fe': [0.12, 0.16, 0.15],
     'k__A;p__Bf;c__Cf;o__Df;f__Ef;g__Ff;s__Gf': [0.20, 0.21, 0.25],
     'k__A;p__Bg;c__Cg;o__Dg;f__Eg;g__Fg;s__Gg': [0.08, 0.03, 0.0]},
    index=['s1', 's2', 's3'])

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
