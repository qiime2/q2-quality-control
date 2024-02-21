# ----------------------------------------------------------------------------
# Copyright (c) 2016-2022, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import os.path
import pkg_resources
import shutil

import decimal
import q2templates
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import qiime2
from q2_types.feature_data import DNAIterator
import skbio

_PER_NUM = (lambda x: 1 >= x >= 0, 'between 0 and 1')
_BOOLEAN = (lambda x: type(x) is bool, 'True or False')
_PRESENT = (lambda x: x is not None, '')
# Better to choose to skip, than to implicitly ignore things that KeyError
_SKIP = (lambda x: True, '')
_valid_inputs = {
    'output_dir': _SKIP,
    'table': _SKIP,
    'decontam_scores': _SKIP,
    'threshold': _PER_NUM,
    'bin_size': _PER_NUM,
    'weighted': _BOOLEAN,
    'rep_seqs': _SKIP
}




def _check_inputs(**kwargs):
    for param, arg in kwargs.items():
        check_is_valid, explanation = _valid_inputs[param]
        if not check_is_valid(arg):
            raise ValueError('Argument to %r was %r, should be %s.'
                             % (param, arg, explanation))


TEMPLATES = pkg_resources.resource_filename(
    'q2_quality_control._threshold_graph', 'assets')
_blast_url_template = ("http://www.ncbi.nlm.nih.gov/BLAST/Blast.cgi?"
                       "ALIGNMENT_VIEW=Pairwise&PROGRAM=blastn&DATABASE"
                       "=nt&CMD=Put&QUERY=%s")

def decontam_score_viz(output_dir, decontam_scores: qiime2.Metadata,
                       rep_seqs: qiime2.Metadata,
                       table: pd.DataFrame, threshold: float = 0.1,
                       weighted: bool = True, bin_size: float = 0.02):
    _check_inputs(**locals())
    # initalizes dictionaries for iteration
    table_dict = dict(table)
    decontam_scores_dict = dict(decontam_scores)
    rep_seqs_df = rep_seqs.to_dataframe()

    # intializes arrays to pass data to the html
    image_paths_arr = []
    subset_key_arr = []
    contam_val_arr = []
    true_val_arr = []
    unknown_val_arr = []
    percent_val_arr = []
    gray_lab_arr = []
    red_lab_arr = []
    blue_lab_arr = []
    data_arr = []
    display_seqs_arr = []
    true_fasta_dest = []
    contam_fasta_dest = []
    nan_fasta_dest = []

    for key in table_dict.keys():
        table = table_dict[key]
        decontam_scores = decontam_scores_dict[key]
        df = decontam_scores.to_dataframe()
        if df['p'].isna().all():
            raise ValueError("No p-values exist for the data provided.")

        read_nums = table.sum(axis='rows')

        #indicies of nan p values
        nan_indices = df[df['p'].isna()].index.tolist()
        na_rep_seqs = rep_seqs_df[rep_seqs_df.index.isin(nan_indices)]['Sequence'].tolist()

        p_vals = df['p'].dropna()
        filt_read_nums = read_nums[p_vals.index]

        contams = (p_vals < threshold)


        contam_indices = contams.index[contams == True].tolist()
        contam_rep_seqs = rep_seqs_df[rep_seqs_df.index.isin(contam_indices)]['Sequence'].tolist()

        true_indices = contams.index[contams == False].tolist()
        true_rep_seqs = rep_seqs_df[rep_seqs_df.index.isin(true_indices)]['Sequence'].tolist()

        true_asvs_iter = DNAIterator((skbio.DNA(true_rep_seqs[i],
                              metadata={'id': true_indices[i]})
                    for i in range(0, len(true_indices))))

        na_asvs_iter = DNAIterator((skbio.DNA(na_rep_seqs[i],
                              metadata={'id': nan_indices[i]})
                    for i in range(0, len(nan_indices))))

        contam_asvs_iter = DNAIterator((skbio.DNA(contam_rep_seqs[i],
                                                metadata={'id': contam_indices[i]})
                                      for i in range(0, len(contam_indices))))
        

        #temp area
        display_sequences = set()
        sequences = {}
        if(len(table_dict.keys()) > 1):
            true_dest = str(key) + '_true_ASV_seqs.fasta'
            contam_dest = str(key) + '_contam_ASV_seqs.fasta'
            nan_dest = str(key) + '_na_ASV_seqs.fasta'
        else:
            true_dest = 'true_ASV_seqs.fasta'
            contam_dest = 'contam_ASV_seqs.fasta'
            nan_dest = 'na_ASV_seqs.fasta'
        with open(os.path.join(output_dir, true_dest), 'w') as fh:
            for sequence in true_asvs_iter:
                skbio.io.write(sequence, format='fasta', into=fh)
                str_seq = str(sequence)
                display_sequences.add(sequence.metadata['id'])
                sequences[sequence.metadata['id']] \
                    = {'url': _blast_url_template % str_seq,
                       'seq': str_seq,
                       'contam_or_naw': 'True_ASV',
                       'p_val': df.loc[sequence.metadata['id'], 'p'],
                       'read_nums': read_nums.loc[sequence.metadata['id']]}
        with open(os.path.join(output_dir, contam_dest), 'w') as fh:
            for sequence in contam_asvs_iter:
                skbio.io.write(sequence, format='fasta', into=fh)
                str_seq = str(sequence)
                display_sequences.add(sequence.metadata['id'])
                sequences[sequence.metadata['id']] \
                    = {'url': _blast_url_template % str_seq,
                       'seq': str_seq,
                       'contam_or_naw': 'Contaminant_ASV',
                       'p_val': df.loc[sequence.metadata['id'], 'p'],
                       'read_nums': read_nums.loc[sequence.metadata['id']]}
        with open(os.path.join(output_dir, nan_dest), 'w') as fh:
            for sequence in na_asvs_iter:
                skbio.io.write(sequence, format='fasta', into=fh)
                str_seq = str(sequence)
                display_sequences.add(sequence.metadata['id'])
                sequences[sequence.metadata['id']] \
                    = {'url': _blast_url_template % str_seq,
                       'seq': str_seq,
                       'contam_or_naw': 'Unknown',
                       'p_val': df.loc[sequence.metadata['id'], 'p'],
                       'read_nums': read_nums.loc[sequence.metadata['id']]}
        #temp area end

        contam_asvs = contams.sum()
        true_asvs = len(contams) - contam_asvs
        unknown_asvs = len(df['p']) - true_asvs - contam_asvs
        percent_asvs = contam_asvs / (contam_asvs + true_asvs) * 100

        contam_reads = filt_read_nums[contams[contams].index].sum()
        true_reads = filt_read_nums.sum() - contam_reads
        unknown_reads = read_nums.sum() - true_reads - contam_reads
        percent_reads = contam_reads / (contam_reads + true_reads) * 100

        binwidth = bin_size
        bin_diff = threshold - (binwidth * int(threshold/binwidth))
        temp_dec = decimal.Decimal(str(binwidth))
        num_dec = abs(temp_dec.as_tuple().exponent)
        temper_dec = decimal.Decimal(str(threshold))
        numer_dec = abs(temper_dec.as_tuple().exponent)
        if numer_dec > num_dec:
            num_dec = numer_dec
        bin_diff = round(bin_diff, num_dec)
        threshold = round(threshold, num_dec)
        lower_bound = round(((threshold - bin_diff)), num_dec)
        bins = np.concatenate([
            np.arange((0.0-(binwidth*2)), (1.0+(binwidth*2)), binwidth)
        ])

        if weighted is True:
            y_lab = 'Number of Reads'
            blue_lab = "True Reads"
            red_lab = "Contaminant Reads"
            gray_lab = "Unknown Reads"
            contam_val = contam_reads
            true_val = true_reads
            unknown_val = unknown_reads
            percent_val = percent_reads
            h, bins, patches = plt.hist(p_vals, bins, weights=filt_read_nums)
            plt.yscale('log')
        else:
            y_lab = 'number of ASVs'
            blue_lab = "True ASVs"
            red_lab = "Contaminant ASVs"
            gray_lab = "Unknown ASVs"
            contam_val = contam_asvs
            true_val = true_asvs
            unknown_val = unknown_asvs
            percent_val = percent_asvs
            h, bins, patches = plt.hist(p_vals, bins)

        plt.xlim(0.0, 1.0)
        plt.xlabel('Score Value')
        plt.ylabel(y_lab)
        arr_bins = list(bins)
        rounded_bins = [round(number, num_dec) for number in arr_bins]
        if threshold in rounded_bins:
            plt.setp([p for p, b in zip(patches, bins) if b < threshold],
                     color='r', edgecolor="white", label=red_lab)
            plt.setp([p for p, b in zip(patches, bins) if b >= threshold],
                     color='b', edgecolor="white", label=blue_lab)
        else:
            plt.setp([p for p, b in zip(patches, bins)
                      if b < threshold and b > (threshold - bin_size)],
                     color='m',
                     edgecolor="white")
            plt.setp([p for p, b in zip(patches, bins)
                      if b < lower_bound], color='r', edgecolor="white",
                     label=red_lab)
            plt.setp([p for p, b in zip(patches, bins)
                      if b >= threshold], color='b', edgecolor="white",
                     label=blue_lab)

        plt.axvline(threshold, ymin=-.1, ymax=1.1, color='k',
                    linestyle='dashed', linewidth=1, label="Threshold")
        handles, labels = plt.gca().get_legend_handles_labels()
        by_label = dict(zip(labels, handles))
        plt.legend(by_label.values(), by_label.keys(),
                   loc="upper left", framealpha=1)

        subset_key_arr.append(key)
        image_prefix = key + '-'
        for ext in ['png', 'svg']:
            img_fp = os.path.join(output_dir,
                                  image_prefix +
                                  'identify-table-histogram.%s' % ext)
            if ext == 'png':
                image_paths_arr.append(
                    './' + image_prefix + 'identify-table-histogram.png'
                )
            plt.savefig(img_fp)
        plt.clf()

        # increments arrays for passing to html
        contam_val_arr.append("{:.0f}".format(contam_val))
        true_val_arr.append("{:.0f}".format(true_val))
        unknown_val_arr.append("{:.0f}".format(unknown_val))
        percent_val_arr.append("%.2f" % percent_val)
        gray_lab_arr.append(gray_lab)
        red_lab_arr.append(red_lab)
        blue_lab_arr.append(blue_lab)
        data_arr.append(sequences)
        display_seqs_arr.append(display_sequences)
        true_fasta_dest.append(true_dest)
        contam_fasta_dest.append(contam_dest)
        nan_fasta_dest.append(nan_dest)

        


    index_fp = os.path.join(TEMPLATES, 'index.html')
    q2templates.render(index_fp, output_dir, context={
            'contamer': contam_val_arr,
            'truer': true_val_arr,
            'unknownr': unknown_val_arr,
            'percenter': percent_val_arr,
            'unknown_label': gray_lab_arr,
            'contam_label': red_lab_arr,
            'true_label': blue_lab_arr,
            'image_paths': image_paths_arr,
            'subset_id': subset_key_arr,
            'data_arr': data_arr,
            'display_sequences_arr': display_seqs_arr,
            'true_fastas': true_fasta_dest,
            'contam_fastas': contam_fasta_dest,
            'na_fastas': nan_fasta_dest,
    })
    js = os.path.join(
        TEMPLATES, 'js', 'tsorter.min.js')
    os.mkdir(os.path.join(output_dir, 'js'))
    shutil.copy(js, os.path.join(output_dir, 'js', 'tsorter.min.js'))
