# ----------------------------------------------------------------------------
# Copyright (c) 2017, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from os.path import join
from scipy.stats import linregress
import seaborn as sns
from itertools import cycle
import matplotlib.patches as mpatches
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import q2_taxa
import q2templates
import pkg_resources


TEMPLATES = pkg_resources.resource_filename('q2_quality_control', 'assets')


def _load_metadata(metadata):
    metadata = metadata.to_series()
    metadata = pd.to_numeric(metadata, errors='ignore')
    return metadata


def _validate_metadata_is_superset(metadata, table):
    metadata_ids = set(metadata.index)
    table_ids = set(table.index.tolist())
    if not table_ids.issubset(metadata_ids):
        raise ValueError('Missing samples in metadata: %r' %
                         table_ids.difference(metadata_ids))


def _validate_metadata_values_are_subset(metadata, table):
    # pull unique exp IDs (metadata vals) from metadata for comparison
    metadata_ids = set(metadata.unique())
    table_ids = set(table.index.tolist())
    if not metadata_ids.issubset(table_ids):
        raise ValueError('Missing samples in table: %r' %
                         table_ids.difference(metadata_ids))


def _interpret_metric_selection(plot_tar, plot_tdr, plot_r_value,
                                plot_r_squared, plot_observed_features,
                                plot_observed_features_ratio):
    # convert metric boolean choices to list of column names
    yvals = []
    for var, val in zip(
            (plot_tar, plot_tdr, plot_r_value, plot_r_squared,
             plot_observed_features, plot_observed_features_ratio),
            ('TAR', 'TDR', 'r-value', 'r-squared', 'Observed Taxa',
             'Observed / Expected Taxa')):
        if var:
            yvals.append(val)

    # at least one metric must be plotted
    if len(yvals) < 1:
        raise ValueError("At least one metric must be plotted.")

    return yvals


def _evaluate_composition(exp, obs, depth, palette, metadata, plot_tar,
                          plot_tdr, plot_r_value, plot_r_squared,
                          plot_observed_features,
                          plot_observed_features_ratio):

    yvals = _interpret_metric_selection(
        plot_tar, plot_tdr, plot_r_value, plot_r_squared,
        plot_observed_features, plot_observed_features_ratio)

    # If metadata are passed, validate and convert to series
    if metadata is not None:
        metadata = _load_metadata(metadata)
        # validate that metadata ids are superset of obs and values (exp ids)
        # are subset of exp
        _validate_metadata_is_superset(metadata, obs)
        _validate_metadata_values_are_subset(metadata, exp)
        # DROP NANS/ZERO ABUNDANCE FEATURES
        obs = _drop_nans_zeros(obs)
        exp = _drop_nans_zeros(exp)

    # if no metadata are passed, we assume that sample IDs correspond directly
    # between obs and exp tables
    else:
        # DROP MISSING SAMPLES/NANS/ZERO ABUNDANCE FEATURES
        exp, obs = _match_samples_by_index(exp, obs)
        obs = _drop_nans_zeros(obs)
        exp = _drop_nans_zeros(exp)

    # TAR/TDR for obs vs. exp at each level
    results, vectors = _compute_per_level_accuracy(exp, obs, metadata, depth)

    # regplot of taxa at top level
    composition_regression = _regplot_from_dict(vectors, palette=palette)

    # pointplot on obs vs. exp at maximum level (coming off of the prior loop)
    score_plot = _pointplot_multiple_y(
        results, xval='level', yvals=yvals, palette=palette)

    # list taxa that are obs but not exp
    obs_collapsed = _collapse_table(obs, depth)
    exp_collapsed = _collapse_table(exp, depth)
    obs_features = set(obs_collapsed.columns)
    exp_features = set(exp_collapsed.columns)
    fp_features, fn_features = _identify_incorrect_classifications(
        obs_features, exp_features)

    misclassifications, underclassifications, mismatches = \
        _tally_misclassifications(fp_features, exp_features)

    # PLOT depth of mismatch at maximum level as histogram
    mismatch_histogram = _plot_histogram(mismatches)

    # generate tables of false postive/negative taxa. Sort alphabetically.
    fn_features = exp_collapsed[list(fn_features)].T.sort_index()
    misclassifications = obs_collapsed[misclassifications].T.sort_index()
    underclassifications = obs_collapsed[underclassifications].T.sort_index()

    return (results, fn_features, misclassifications, underclassifications,
            composition_regression, score_plot, mismatch_histogram)


def _match_samples_by_index(df, df2):
    # drop all rows (samples) in df that do not match df2, and vice versa
    df = df.loc[df2.index]
    df2 = df2.loc[df.index]
    return df, df2


def _drop_nans_zeros(df):
    # replace nan with zero
    df = df.fillna(0)
    # drop rows / cols with all zero values
    df = df.loc[(df.sum(axis=1) != 0), (df.sum(axis=0) != 0)]
    # remove spaces from rownames (taxonomy classifiers are
    # inconsistent about adding/removing spaces so let's just be safe)
    replacements = {t: t.replace(" ", "") for t in df.columns}
    return df.rename(columns=replacements, inplace=False)


def _compute_per_level_accuracy(exp, obs, metadata, depth):
    results = []
    vectors = {}
    for level in range(1, depth + 1):
        vectors[level] = {'exp': [], 'obs': []}
        # collapse taxonomy strings to level
        exp_collapsed = _collapse_table(exp, level)
        obs_collapsed = _collapse_table(obs, level)
        # compute stats for each sample individually
        for sample in obs_collapsed.index:
            result = [sample, level]
            # if metadata are passed, map exp sample ID to value in metadata
            if metadata is not None:
                exp_id = metadata[sample]
            else:
                exp_id = sample
            # concatenate obs/exp observations to align features
            joined_table = pd.concat(
                [exp_collapsed.loc[exp_id],
                 obs_collapsed.loc[sample]], axis=1).fillna(0)
            # split joined table apart again for computing stats
            exp_vector = joined_table.iloc[:, 0]
            obs_vector = joined_table.iloc[:, 1]
            exp_features = exp_vector[exp_vector != 0]
            obs_features = obs_vector[obs_vector != 0]
            # Count observed taxa
            observed_feature_count = len(obs_features)
            observed_feature_ratio = (
                observed_feature_count / len(exp_features))
            result.extend([observed_feature_count, observed_feature_ratio])
            # compute TAR/TDR
            result.extend(compute_taxon_accuracy(exp_features, obs_features))
            # compute linear least-squares regression results
            result.extend(linregress(exp_vector, obs_vector))
            results.append(result)
            # store vectors for constructing regplots
            vectors[level]['exp'].extend(exp_vector)
            vectors[level]['obs'].extend(obs_vector)
    results = pd.DataFrame(
        results, columns=['sample', 'level', 'Observed Taxa',
                          'Observed / Expected Taxa', 'TAR', 'TDR', 'Slope',
                          'Intercept', 'r-value', 'P value', 'Std Err'])
    results['r-squared'] = results['r-value']**2
    return results, vectors


# ported and modified from tax-credit with permission of nbokulich
def compute_taxon_accuracy(exp, obs):
    actual_obs_ids = set(exp.index)
    expected_obs_ids = set(obs.index)

    tp = len(actual_obs_ids & expected_obs_ids)
    fp = len(actual_obs_ids - expected_obs_ids)
    fn = len(expected_obs_ids - actual_obs_ids)

    if tp > 0:
        p = tp / (tp + fp)
        r = tp / (tp + fn)
    else:
        p, r = 0, 0

    return p, r


def _find_nearest_common_lineage(feature, exp_features):
    feature_depth = len(feature.split(';'))
    # slice off empty labels
    feature = feature.rstrip(';_').split(';')
    for i in range(0, len(feature), 1):
        # start with None (i.e., don't remove any elements from list)
        # this will check for underclassifications first
        if i == 0:
            i = None
        else:
            i = -i
        for f in exp_features:
            if ';'.join(feature[:i]) in f:
                # underclassified features will be substring of exp feature(s)
                # after stripping
                if i is None:
                    # return difference in length between obs feature and exp
                    # match (negative for underclassification)
                    return len(feature) - len(f.split(';'))
                # misclassified features only match after scraping tip lineages
                else:
                    # return positive number: degree of misclassification
                    # (includes overclassification)
                    return -i
    # if no matches are found anywhere, return full length (no common lineages)
    return feature_depth


def _pointplot_multiple_y(results, xval, yvals, palette):
    colors = cycle(sns.color_palette(palette, n_colors=len(yvals)))
    fig, axes = plt.subplots(1)
    handles = []
    for score in yvals:
        color = next(colors)
        sns.pointplot(data=results, x=xval, y=score, ax=axes, color=color)
        handles.append(mpatches.Patch(color=color, label=score))
    axes.set_ylabel('Score')
    axes.set_xlabel('Taxonomic level')
    axes.legend(
        handles=handles, bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    return fig


def _regplot_from_dict(vectors, palette):
    yvals = vectors.keys()
    colors = cycle(sns.color_palette(palette, n_colors=len(yvals)))
    fig, axes = plt.subplots(1)
    handles = []
    for level in yvals:
        color = next(colors)
        sns.regplot(np.array(vectors[level]['exp']),
                    np.array(vectors[level]['obs']), ax=axes, color=color)
        handles.append(mpatches.Patch(color=color, label=level))
    axes.legend(
        handles=handles, bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)

    # Plot arbitrary line with slope of 1 (true ratio)
    x0, x1 = axes.axes.get_xlim()
    y0, y1 = axes.axes.get_ylim()
    lims = [min(x0, y0), max(x1, y1)]
    axes.plot(lims, lims, ':k')

    plt.xlabel('Expected abundance')
    plt.ylabel('Observed abundance')

    return fig


def _plot_histogram(mismatches):
    fig, axes = plt.subplots(1)
    n, bins, patches = plt.hist(mismatches, align='left')
    plt.ylabel('Count')
    plt.xlabel('Distance to nearest expected feature')
    return fig


def _collapse_table(table, level):
    return q2_taxa._method.collapse(table, pd.Series(
        table.columns, index=table.columns, name='Taxon'), level)


def _identify_incorrect_classifications(obs_features, exp_features):
    fp_features = obs_features - exp_features
    fn_features = exp_features - obs_features
    return fp_features, fn_features


def _tally_misclassifications(fp_features, exp_features):
    misclassifications = []
    underclassifications = []
    mismatches = []
    if len(fp_features) > 0:
        # find and report underclassifications and over/misclassification
        for feature in fp_features:
            mismatch = _find_nearest_common_lineage(feature, exp_features)
            mismatches.append(mismatch)
            if mismatch > 0:
                misclassifications.append(feature)
            else:
                underclassifications.append(feature)
    return (sorted(misclassifications), sorted(underclassifications),
            sorted(mismatches))


def _visualize(output_dir, results, false_negative_features=None,
               misclassifications=None, underclassifications=None,
               composition_regression=None, score_plot=None,
               mismatch_histogram=None):

    pd.set_option('display.max_colwidth', -1)

    # save results
    results.to_csv(join(output_dir, 'results.tsv'), sep='\t')
    results = q2templates.df_to_html(results, index=False)

    if false_negative_features is not None:
        false_negative_features.to_csv(join(
            output_dir, 'false_negative_features.tsv'), sep='\t')
        false_negative_features = q2templates.df_to_html(
            false_negative_features, index=True)

    if misclassifications is not None:
        misclassifications.to_csv(join(
            output_dir, 'misclassifications.tsv'), sep='\t')
        misclassifications = q2templates.df_to_html(
            misclassifications, index=True)

    if underclassifications is not None:
        underclassifications.to_csv(join(
            output_dir, 'underclassifications.tsv'), sep='\t')
        underclassifications = q2templates.df_to_html(
            underclassifications, index=True)

    if composition_regression is not None:
        composition_regression.savefig(join(
            output_dir, 'composition_regression.png'), bbox_inches='tight')
        composition_regression.savefig(join(
            output_dir, 'composition_regression.pdf'), bbox_inches='tight')

    if score_plot is not None:
        score_plot.savefig(join(
            output_dir, 'score_plot.png'), bbox_inches='tight')
        score_plot.savefig(join(
            output_dir, 'score_plot.pdf'), bbox_inches='tight')

    if mismatch_histogram is not None:
        mismatch_histogram.savefig(join(
            output_dir, 'mismatch_histogram.png'), bbox_inches='tight')
        mismatch_histogram.savefig(join(
            output_dir, 'mismatch_histogram.pdf'), bbox_inches='tight')

    index = join(TEMPLATES, 'index.html')
    q2templates.render(index, output_dir, context={
        'title': 'evaluate_composition',
        'results': results,
        'false_negative_features': false_negative_features,
        'misclassifications': misclassifications,
        'underclassifications': underclassifications,
        'composition_regression': composition_regression,
        'score_plot': score_plot,
        'mismatch_histogram': mismatch_histogram,
    })
