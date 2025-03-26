# ----------------------------------------------------------------------------
# Copyright (c) 2017-2025, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
def decontam_identify_batches(ctx, table, metadata,
                              split_column,
                              method,
                              rep_seqs=None,
                              filter_empty_features=None,
                              freq_concentration_column=None,
                              prev_control_column=None,
                              prev_control_indicator=None,
                              threshold=0.1,
                              weighted=True,
                              bin_size=0.02):

    # Gets Actions for pipeline and sets up results arrays
    decon = ctx.get_action('quality_control', 'decontam_identify')
    decon_score_viz = ctx.get_action('quality_control', 'decontam_score_viz')
    spliter = ctx.get_action('feature_table', 'split')

    # New Work Flow
    split_tables, = spliter(table=table,
                            metadata=metadata.get_column(split_column),
                            filter_empty_features=filter_empty_features)
    table_col = split_tables.collection
    table_dic = dict(table_col)

    # deletes NA sample subset
    if ('NA' in table_dic.keys()):
        del table_dic["NA"]

    split_tables_dict = {}
    decon_results = {}
    for key in table_dic:
        # Adds subtable to outpur object
        keyer = str(split_column) + '-' + key
        split_tables_dict[keyer] = table_dic[key]
        # adds decon_scores from sub table to output object
        temp_results, = decon(
            table=split_tables_dict[keyer],
            metadata=metadata, method=method,
            freq_concentration_column=freq_concentration_column,
            prev_control_column=prev_control_column,
            prev_control_indicator=prev_control_indicator)
        decon_results[keyer] = temp_results
    temp_viz_results, = decon_score_viz(
        decontam_scores=decon_results,
        rep_seqs=rep_seqs,
        table=split_tables_dict, threshold=threshold,
        weighted=weighted, bin_size=bin_size)

    return split_tables_dict, decon_results, temp_viz_results
