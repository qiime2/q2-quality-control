# ----------------------------------------------------------------------------
# Copyright (c) 2017-2021, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import pandas as pd


def _evaluate_taxonomy(exp_taxa, obs_taxa, require_exp_ids=True,
                       require_obs_ids=True, feature_table=None,
                       sample_id=None, level_range=range(0, 7)):
    # validate inputs
    join_how = _validate_indices_and_set_joining_mode(
        exp_taxa, obs_taxa, require_exp_ids, require_obs_ids)

    # merge tables
    taxa = exp_taxa.join(
        obs_taxa, how=join_how, lsuffix='_exp', rsuffix='_obs')

    # extract class weights from feature table (if one is supplied)
    weights = _extract_frequencies_from_feature_table(
        taxa, feature_table, sample_id)

    # calculate precision/recall
    prf = _prf_to_dataframe(taxa['Taxon_exp'], taxa['Taxon_obs'],
                            sample_weight=weights, level_range=level_range)

    if sample_id is not None:
        prf['SampleID'] = sample_id

    return prf


# modified from tax-credit with permission of nbokulich
def _compute_prf(exp, obs, l_range=range(0, 7), sample_weight=None):
    p, r, f = {}, {}, {}
    # iterate over multiple taxonomic levels
    for lvl in l_range:
        lvl = lvl + 1
        _obs = _extract_taxa_names(obs, level=slice(0, lvl))
        _exp = _extract_taxa_names(exp, level=slice(0, lvl))
        p[lvl], r[lvl], f[lvl] = _precision_recall_fscore(
            _exp, _obs, sample_weight=sample_weight)

    return p, r, f


# modified from tax-credit with permission of nbokulich
def _extract_taxa_names(inlist, level=slice(6, 7), delim=';', stripchars=None):
    '''Extract taxon names at a given level from list of taxon names

    stripchars: str
        Default, None, will strip leading and trailing whitespace. Set to ""
        to turn off character stripping.
    '''
    # Truncate taxonomies and pass to set
    name_list = [delim.join(line.split(delim)[level]).strip(stripchars)
                 for line in inlist]
    return name_list


# ported from tax-credit with permission of nbokulich
def _precision_recall_fscore(exp, obs, sample_weight=None):
    # precision, recall, fscore, calculated using microaveraging
    if sample_weight is None:
        sample_weight = [1]*len(exp)
    tp, fp, fn = 0, 0, 0
    for e, o, w in zip(exp, obs, sample_weight):
        if o == e:
            # tp for the true class, the rest are tn
            tp += w
        elif e.startswith(o) or o in (
                'Unclassified', 'Unassigned', 'No blast hit', 'other'):
            # fn for the true class
            fn += w
            # no fp for the predicted class, because it was right to some level
            # the rest are tn
        else:
            # fp the the predicted class
            fp += w
            # fn for the true class, the rest are tn
            fn += w

    # avoid divide by zero error. If no true positives, all scores = 0
    if tp == 0:
        return 0, 0, 0
    else:
        p = tp / (tp + fp)
        r = tp / (tp + fn)
        f = 2.*p*r / (p + r)

    return p, r, f


def _index_is_subset(series1, series2, name):
    ix1 = series1.index
    ix2 = series2.index
    if set(ix1) < set(ix2):
        raise ValueError(
            'Observed and expected ids do not match. Missing '
            'ids not found in {0} ids: {1}'.format(
                name, set(ix1) - (set(ix2))))


def _prf_to_dataframe(exp, obs, sample_weight=None, level_range=range(0, 7)):
    p, r, f = _compute_prf(
        exp, obs, l_range=level_range, sample_weight=sample_weight)
    prf = pd.DataFrame([p, r, f]).T
    prf = prf.reset_index()
    prf.columns = ['level', 'Precision', 'Recall', 'F-measure']
    return prf


def _validate_indices_and_set_joining_mode(exp_taxa, obs_taxa, require_exp_ids,
                                           require_obs_ids):
    if require_exp_ids:
        join_how = 'left'
        _index_is_subset(obs_taxa, exp_taxa, 'observed')

    if require_obs_ids:
        join_how = 'right'
        _index_is_subset(exp_taxa, obs_taxa, 'expected')

    if require_exp_ids == require_obs_ids:
        join_how = 'inner'

    return join_how


def _extract_frequencies_from_feature_table(taxa, feature_table, sample_id):
    if feature_table is not None:
        taxa_ix = taxa.index
        diff = set(taxa_ix) - set(feature_table.ids(axis='observation'))
        if len(diff) != 0:
            raise ValueError('Feature ids not found in feature table: '
                             '{0}'.format(diff))
        table = feature_table.filter(taxa_ix, axis='observation')
        table = table.sort_order(taxa_ix, axis='observation')

        if sample_id is not None:
            # confirm that sample_id is in table
            if sample_id not in table.ids(axis='sample'):
                raise ValueError('Sample id not found in feature table: '
                                 '{0}'.format(sample_id))
            weights = table.data(sample_id, axis='sample')

        else:
            weights = table.sum(axis='observation')

    else:
        weights = None

    return weights
