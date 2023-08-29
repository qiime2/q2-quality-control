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
from q2_types.feature_table import FeatureTable, Frequency

def decontam_identify_batches(ctx, table, metadata,
                              split_columns,
                              filter_empty_features,
                              method,
                              freq_concentration_column,
                              prev_control_column,
                              prev_control_indicator):

    #Gets Actions for pipeline and sets up results arrays
    decon = ctx.get_action('quality_control', 'decontam_identify')
    spliter = ctx.get_action('feature_table', 'split')
    results = []

    #Generates list of all possilbe subset combinations and sorts them from longest to shortest
    split_column_arr = split_columns.split(' ')
    list_combinations = list()
    for n in range(len(split_column_arr) + 1):
        list_combinations += list(combinations(split_column_arr, n))
    list_combinations.pop(0)
    list_combinations = sorted(list_combinations, key=lambda l: (len(l), l), reverse=True)
    print(list_combinations)

    #Intitates action algorithm
    split_tables_dict = {'original': table}
    already_done_arr = []
    for split_col in list_combinations:
        subject_table_dict = {'original': table}
        print("---------------------------------------------------------")
        if split_col not in already_done_arr:
            print(split_col)
            for inter_col in split_col:
                temp_dict = {}
                for subject_table_key in subject_table_dict.keys():
                    subject_metadata = metadata

                    # Gets apporpriate table and metadata for splitter method
                    subject_table = subject_table_dict[subject_table_key]
                    df = subject_table.view(pd.DataFrame)
                    metadata_df = subject_metadata.to_dataframe()
                    metadata_df = metadata_df[metadata_df.index.isin(df.index)]

                    #checks if the subset created a table with no entries
                    subject_metadata = qiime2.Metadata(metadata_df)
                    split_tables, = spliter(table=subject_table, metadata=subject_metadata.get_column(inter_col),
                                       filter_empty_features=filter_empty_features)
                    table_col = split_tables.collection
                    table_dic = dict(table_col)

                    #deltes NA sample subset
                    if('NA' in table_dic.keys()):
                        del table_dic["NA"]

                    #checks for base case vs exponential case
                    if subject_table_key == 'original':
                        subject_table_key = ''
                    else:
                        subject_table_key = subject_table_key + '_'

                    #updates nomenclature
                    for key in table_dic:
                        keyer = subject_table_key + inter_col + '-' + key
                        temp_dict[keyer] = table_dic[key]
                subject_table_dict = temp_dict
                split_tables_dict.update(subject_table_dict)
                print(split_tables_dict.keys())
            index = 0
            #print(split_col)
            #print(already_done_arr)
            while(index <= (len(split_col)-1)):
                temp_arr = split_col[0:((len(split_col)-1-index))]
                already_done_arr.append(temp_arr)
                index = index + 1
            #print(already_done_arr)
            print("---------------------------------------------------------")

    new_table_dic = {}
    decon_results = {}
    for keyer in split_tables_dict.keys():
        new_key_feature = ('ASV-table_' + str(keyer))
        new_table_dic[new_key_feature] = split_tables_dict[keyer]
        sub_table = new_table_dic[new_key_feature]
        if method == 'combined':
            temp_results, = decon(table=sub_table, metadata=metadata, method=method,
                    freq_concentration_column=freq_concentration_column,
                    prev_control_column=prev_control_column,
                    prev_control_indicator=prev_control_indicator)
            decon_results[('decon-scores_' + str(keyer))] = temp_results
        elif method == 'prevalence':
            temp_results, = decon(table=sub_table, metadata=metadata, method=method,
                    prev_control_column=prev_control_column,
                    prev_control_indicator=prev_control_indicator)
            decon_results[('decon-scores_' + str(keyer))] = temp_results
        else:
            temp_results, = decon(table=sub_table, metadata=metadata, method=method,
                    freq_concentration_column=freq_concentration_column)
            decon_results[('decon-scores_' + str(keyer))] = temp_results

    #split_tables.collection = new_table_dic
    results.append(new_table_dic)
    results.append(decon_results)
    return tuple(results)


def decontam_fast(ctx, table,
                      metadata,
                      method,
                      freq_concentration_column,
                      prev_control_column,
                      prev_control_indicator, threshold):
    decon = ctx.get_action('quality_control', 'decontam_identify')
    decon_rem = ctx.get_action('quality_control', 'decontam_remove')
    results = []
    if method == 'combinded':
        identify_temp_output, = decon(table=table, metadata=metadata, method=method,
                            freq_concentration_column=freq_concentration_column,
                            prev_control_column=prev_control_column,
                            prev_control_indicator=prev_control_indicator)
        results.append(identify_temp_output)
        results += decon_rem(decontam_scores=identify_temp_output, table=table, threshold=threshold)
    elif method == 'prevalence':
        identify_temp_output, = decon(table=table, metadata=metadata, method=method,
                            prev_control_column=prev_control_column,
                            prev_control_indicator=prev_control_indicator)
        results.append(identify_temp_output)
        results += decon_rem(decontam_scores=identify_temp_output, table=table, threshold=threshold)
    else:
        identify_temp_output, = decon(table=table, metadata=metadata, method=method,
                            freq_concentration_column=freq_concentration_column)
        results.append(identify_temp_output)
        results += decon_rem(decontam_scores=identify_temp_output, table=table, threshold=threshold)
    return tuple(results)