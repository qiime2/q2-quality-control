# ----------------------------------------------------------------------------
# Copyright (c) 2017-2022, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
from itertools import combinations
import pandas as pd
from qiime2.plugin.util import transform
import qiime2
import biom
import numpy as np
import re
from q2_quality_control._stats import DecontamScoreFormat
from q2_types.feature_table import FeatureTable, Frequency
from qiime2.plugin import (Str, Plugin, Choices, Range, Float, Int, Bool,
                           MetadataColumn, Visualization,Categorical, Citations, TypeMap,
                           Visualization, TypeMatch, Metadata, Collection, List)


def decontam_identify_batches(ctx, table, metadata,
                              split_column,
                              method,
                              filter_empty_features=None,
                              freq_concentration_column=None,
                              prev_control_column=None,
                              prev_control_indicator=None,
                              threshold=0.1,
                              weighted=True,
                              bin_size=0.02):

    #Gets Actions for pipeline and sets up results arrays
    decon = ctx.get_action('quality_control', 'decontam_identify')
    decon_score_viz = ctx.get_action('quality_control', 'decontam_score_viz')
    spliter = ctx.get_action('feature_table', 'split')
    results = []

    #New Work Flow
    split_tables, = spliter(table=table, metadata=metadata.get_column(split_column),
                            filter_empty_features=filter_empty_features)
    table_col = split_tables.collection
    table_dic = dict(table_col)

    # deltes NA sample subset
    if ('NA' in table_dic.keys()):
        del table_dic["NA"]

    split_tables_dict = {}
    decon_results = {}
    decon_viz_table_df = pd.DataFrame()
    decon_viz_decon_df = pd.DataFrame()
    for key in table_dic:
        #Adds subtable to outpur object
        keyer = 'ASV-table_' + str(split_column) + '-' + key
        split_tables_dict[keyer] = table_dic[key]

        temp_table = split_tables_dict[keyer]
        table_df = temp_table.view(pd.DataFrame)
        blank_df = pd.DataFrame(np.nan, index=[0], columns=table_df.columns)
        new_index_value = 'this_is_new_table_' + split_column + '_' + key
        blank_df.index = [new_index_value] + blank_df.index[1:].tolist()
        decon_viz_table_df = pd.concat([decon_viz_table_df, blank_df, temp_table.view(pd.DataFrame)], axis=0)

        #adds decon_scores from sub table to output object
        temp_results, = decon(table=split_tables_dict[keyer], metadata=metadata, method=method,
                    freq_concentration_column=freq_concentration_column,
                    prev_control_column=prev_control_column,
                    prev_control_indicator=prev_control_indicator)
        decon_results[('decon-scores_' + str(split_column) + '-' + key)] = temp_results

        table_df = temp_results.view(pd.DataFrame)
        prefix = split_column + '_' + key + '_JORDENRABASCO_'
        # Add the prefix to all index values
        table_df.index = [prefix + str(index) for index in table_df.index]

        blank_df = pd.DataFrame(np.nan, index=[0], columns=table_df.columns)
        new_index_value = 'this_is_new_table_' + split_column + '_' + key
        blank_df.index = [new_index_value] + blank_df.index[1:].tolist()
        decon_viz_decon_df = pd.concat([decon_viz_decon_df, blank_df, table_df], axis=0)

    #adds score viz to output object
    decon_viz_table_df.index.name = 'SampleID'  # Set SampleID as the index
    q2_feature_table = qiime2.Artifact.import_data('FeatureTable[Frequency]', decon_viz_table_df)
    decon_viz_decon_df.index.name = '#OTUID'
    decontam_table = qiime2.Artifact.import_data('FeatureData[DecontamScore]', decon_viz_decon_df)

    temp_viz_results, = decon_score_viz(decontam_scores=decontam_table,
                    table=q2_feature_table, threshold=threshold,
                    weighted=weighted, bin_size=bin_size)
        #decon_viz_results[('decon-scores-hist_' + str(split_column) + '-' + key)] = temp_viz_results

    #split_tables.collection = new_table_dic
    results.append(split_tables_dict)
    results.append(decon_results)
    results.append(temp_viz_results)
    #return (split_tables_dict, decon_results, *decon_viz_results)
    return tuple(results)

