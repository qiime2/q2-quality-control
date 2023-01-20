# ----------------------------------------------------------------------------
# Copyright (c) 2016-2022, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import os.path
import pkg_resources

import decimal
import q2templates
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import qiime2


_BOOLEAN = (lambda x: type(x) is bool, 'True or False')

TEMPLATES = pkg_resources.resource_filename('q2_decontam._threshold_graph',
                                            'assets')
def decontam_score_viz(output_dir, decon_identify_table: qiime2.Metadata, asv_or_otu_table: pd.DataFrame, threshold: float=0.1, weighted: bool=True, bin_size: float=0.02):

    df = decon_identify_table.to_dataframe()
    values = df['p'].tolist()
    values = np.array(values)
    temp = asv_or_otu_table.sum(axis='columns')
    read_nums = np.array(temp.tolist())

    contam_asvs = 0
    true_asvs = 0
    contam_reads = 0
    true_reads = 0
    index = 0
    for val in values:
        if val < threshold:
            contam_asvs = contam_asvs + 1
            contam_reads = contam_reads + read_nums[index]
        else:
            true_asvs = true_asvs + 1
            true_reads = true_reads + read_nums[index]
        index = index + 1
    binwidth = bin_size
    bin_diff = threshold % binwidth
    temp_dec = decimal.Decimal(str(binwidth))
    num_dec = abs(temp_dec.as_tuple().exponent)
    bin_diff = round(bin_diff,num_dec)
    bins = np.concatenate([
        np.arange((0.0-(binwidth*2)), (1.0+(binwidth*2)), binwidth)
    ])
    if(weighted == True):
        y_lab = 'Number of Reads'
        blue_lab = "True Reads"
        red_lab = "Contaminant Reads"
        h, bins, patches = plt.hist(values, bins, weights=np.array(temp.tolist()))
        plt.yscale('log')
    else:
        y_lab = 'number of ASVs'
        blue_lab = "True ASVs"
        red_lab = "Contaminant ASVs"
        h, bins, patches = plt.hist(values, bins)

    plt.xlim(0.0, 1.0)
    plt.xlabel('score value')
    plt.ylabel(y_lab)

    if bin_diff == 0:
        plt.setp([p for p, b in zip(patches, bins) if b < (threshold)], color='r', edgecolor="white",
                 label=red_lab)
        plt.setp([p for p, b in zip(patches, bins) if b >= (threshold)], color='b', edgecolor="white",
                 label=blue_lab)
    else:
        plt.setp([p for p, b in zip(patches, bins) if b == (threshold - bin_diff)],
                 color='m', edgecolor="white")
        plt.setp([p for p, b in zip(patches, bins) if b < (threshold-bin_diff)], color='r', edgecolor="white",
                 label=red_lab)
        plt.setp([p for p, b in zip(patches, bins) if b > (threshold)], color='b', edgecolor="white",
                 label=blue_lab)

    plt.axvline(threshold, ymin=-.1,ymax=1.1 ,color='k', linestyle='dashed', linewidth=1, label="Threshold")
    handles, labels = plt.gca().get_legend_handles_labels()
    by_label = dict(zip(labels, handles))
    plt.legend(by_label.values(), by_label.keys(), loc="upper left", framealpha=1)

    percent_reads = (100*float(contam_reads)/float((contam_reads+true_reads)))
    percent_asvs = (100*float(contam_asvs)/float((contam_asvs+true_asvs)))

    for ext in ['png', 'svg']:
        img_fp = os.path.join(output_dir, 'identify-table-histogram.%s' % ext)
        plt.savefig(img_fp)
    index_fp = os.path.join(TEMPLATES, 'index.html')

    if (weighted == True):
        q2templates.render(index_fp, output_dir, context={'contamer': str("{:,}".format(int(contam_reads))), 'truer': str("{:,}".format(int(true_reads))), 'percenter': str("%.2f" % percent_reads),
                                                          'contam_label': str(red_lab), 'true_label': str(blue_lab)})
    else:
        q2templates.render(index_fp, output_dir, context={'contamer': str("{:,}".format(int(contam_asvs))), 'truer': str("{:,}".format(int(true_asvs))), 'percenter': str("%.2f" % percent_asvs),
                                                          'contam_label': str(red_lab), 'true_label': str(blue_lab)})

