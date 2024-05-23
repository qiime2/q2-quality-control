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
from q2_types.feature_data import DNAIterator

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


#checks inputs
def _check_inputs(**kwargs):
    for param, arg in kwargs.items():
        check_is_valid, explanation = _valid_inputs[param]
        if not check_is_valid(arg):
            raise ValueError('Argument to %r was %r, should be %s.'
                             % (param, arg, explanation))


TEMPLATES = pkg_resources.resource_filename(
    'q2_quality_control._threshold_graph', 'assets')

#generates seqeunce table and fasta files for each assignment
def write_table_fastas(output_dir, dest, sequences, seq_list, desig, decontam_scores, read_nums, table):
    _blast_url_template = ("http://www.ncbi.nlm.nih.gov/BLAST/Blast.cgi?"
                           "ALIGNMENT_VIEW=Pairwise&PROGRAM=blastn&DATABASE"
                           "=nt&CMD=Put&QUERY=%s")
    with open(os.path.join(output_dir, dest), 'w') as fh:
        for seq in seq_list:
            seq.write(fh)
            str_seq = str(seq)
            sequences[seq.metadata['id']] \
                = {'url': _blast_url_template % str_seq,
                   'seq': str_seq,
                   'contam_or_naw': desig,
                   'p_val': decontam_scores.loc[seq.metadata['id'], 'p'],
                   'read_nums': read_nums.loc[seq.metadata['id']],
                   'prevalence': (
                           table[seq.metadata['id']] != 0).sum()}
    return sequences

#main algorithm
def decontam_score_viz(output_dir, decontam_scores: pd.DataFrame,
                       table: pd.DataFrame,
                       rep_seqs: DNAIterator = None,
                       threshold: float = 0.1,
                       weighted: bool = True, bin_size: float = 0.02):

    #checks inputs
    _check_inputs(**locals())

    # initalizes dictionaries for iteration
    table_dict = dict(table)
    decontam_scores_dict = dict(decontam_scores)

    # Sets rep seq flags if rep_seq indciator is >1 then the seqeunces are printed in table otherwise no sequences are printed and temp lists
    rep_seq_indicator = ["Are there rep seqs?"]

    # intializes arrays to pass data to the html
    image_paths_arr = [] #array for image paths for render on template (length 1 when running base decontam-score-viz)
    subset_key_arr = [] #array for table subset id when runnning batches (length 1 when running base decontam-score-viz)
    contam_val_arr = [] #contaminant features count (length 1 when running base decontam-score-viz)
    true_val_arr = [] #non-contaminant features count (length 1 when running base decontam-score-viz)
    unknown_val_arr = [] #NA features count (length 1 when running base decontam-score-viz)
    percent_val_arr = [] #% contaminant features (length 1 when running base decontam-score-viz)
    red_lab_arr = [] # contaminant feature/read label  (length 1 when running base decontam-score-viz)
    blue_lab_arr = [] #non-contaminant feature/read label (length 1 when running base decontam-score-viz)
    data_arr = [] # information for sequencing table (only render in base decontam-score-viz)
    true_fasta_dest = [] #destination for true/nan seq fastas (only render if rep-seqs not None)
    contam_fasta_dest = [] #destination for contaminant feature fastas (only render if rep-seqs not None)
    sorted_key_arr = [] # sorted feature ids (only render in base decontam-score-viz)
    feature_or_read_arr = [] # feature or read label for graph rendering (length 1 when running base decontam-score-viz)

    temp_list =[]
    if rep_seqs is not None:
        for seq in rep_seqs:
            temp_list.append(seq)

    #iterates through tables and keys of ASV tables and decontam score tables
    for key in table_dict.keys():
        table = table_dict[key]
        decontam_scores = decontam_scores_dict[key]
        if decontam_scores['p'].isna().all():
            raise ValueError("No p-values exist for the data provided.")

        #pull out relevant data from objs
        read_nums = table.sum(axis='rows')
        p_vals = decontam_scores['p'].dropna()
        filt_read_nums = read_nums[p_vals.index]

        #start ASV contaminant differentiation
        contams = (p_vals < threshold)

        #parses index for contaminant identification
        true_indices = contams[contams == False].index
        contam_indices = contams[contams == True].index
        nan_indices = decontam_scores[decontam_scores['p'].isna()].index.tolist()


        #if rep reqs are not found then dummy variables are initalized otherwise individual seqeunce objects are inalized for true seqe, NA seqs, and contaminant seqs
        contam_rep_seqs=[]
        true_rep_seqs=[]
        if rep_seqs is not None:
            for seq in temp_list:
                if seq.metadata['id'] in contam_indices:
                    contam_rep_seqs.append(seq)
                    continue
                if (seq.metadata['id'] in true_indices) or (seq.metadata['id'] in nan_indices):
                    true_rep_seqs.append(seq)
        else:
            rep_seq_indicator.append("Nope there are not")
            true_rep_seqs = []
            contam_rep_seqs = []

        #initzalized seqeunces for display in table and fasta downloads
        sequences = {}
        if len(table_dict.keys()) > 1:
            true_dest = str(key) + '_non_contam.fasta'
            contam_dest = str(key) + '_contam.fasta'
        else:
            true_dest = 'non_contam.fasta'
            contam_dest = 'contam.fasta'


        #generate repseq table and fasta for non contaminants
        sequences = write_table_fastas(output_dir, true_dest, sequences, true_rep_seqs, "Non-Contaminant", decontam_scores, read_nums, table)

        #generate repseq table and fasta for contaminants
        sequences = write_table_fastas(output_dir, contam_dest, sequences, contam_rep_seqs, "Contaminant", decontam_scores, read_nums, table)

        #sorts sequences to be highest read nums first
        sorted_keys = sorted(
            sequences, key=lambda x: sequences[x]['read_nums'], reverse=True)

        #calculates percentage of contaminant asvs, true and contaminant asv feature nums
        contam_asvs = contams.sum()
        true_asvs = len(contams) - contam_asvs
        unknown_asvs = len(decontam_scores['p']) - true_asvs - contam_asvs
        percent_asvs = contam_asvs / (
                contam_asvs + true_asvs + unknown_asvs) * 100
        true_asvs = unknown_asvs + true_asvs
        contam_reads = filt_read_nums[contams[contams].index].sum()
        true_reads = filt_read_nums.sum() - contam_reads
        unknown_reads = read_nums.sum() - true_reads - contam_reads
        percent_reads = contam_reads / (
                contam_reads + true_reads + unknown_reads) * 100
        true_reads = unknown_reads + true_reads

        #bin width and different color calculations for histogram
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

        #weighted histogram determination for graph
        if weighted is True:
            y_lab = 'Number of Reads'
            blue_lab = "Non-Contaminant Reads"
            red_lab = "Contaminant Reads"
            feature_or_read = "Reads"
            contam_val = contam_reads
            true_val = true_reads
            unknown_val = unknown_reads
            percent_val = percent_reads
            h, bins, patches = plt.hist(p_vals, bins, weights=filt_read_nums)
            plt.yscale('log')
        else:
            y_lab = 'Number of Features'
            blue_lab = "Non-Contaminant Features"
            red_lab = "Contaminant Features"
            feature_or_read = "Features"
            contam_val = contam_asvs
            true_val = true_asvs
            unknown_val = unknown_asvs
            percent_val = percent_asvs
            h, bins, patches = plt.hist(p_vals, bins)

        plt.xlim(0.0, 1.0)
        plt.xlabel('Score')
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

        #threshold line plotting
        plt.axvline(threshold, ymin=-.1, ymax=1.1, color='k',
                    linestyle='dashed', linewidth=1, label="Threshold")
        handles, labels = plt.gca().get_legend_handles_labels()
        by_label = dict(zip(labels, handles))
        plt.legend(by_label.values(), by_label.keys(),
                   loc="upper left", framealpha=1)

        #code for saving plot as PNG for render, svg code kept for future potential additon of .svg file
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
        red_lab_arr.append(red_lab)
        blue_lab_arr.append(blue_lab)
        data_arr.append(sequences)
        true_fasta_dest.append(true_dest)
        contam_fasta_dest.append(contam_dest)
        sorted_key_arr.append(sorted_keys)
        feature_or_read_arr.append(feature_or_read)

    #intializes template and passes data arrays to template for render
    index_fp = os.path.join(TEMPLATES, 'index.html')
    q2templates.render(index_fp, output_dir, context={
            'contamer': contam_val_arr,
            'truer': true_val_arr,
            'unknownr': unknown_val_arr,
            'percenter': percent_val_arr,
            'contam_label': red_lab_arr,
            'true_label': blue_lab_arr,
            'image_paths': image_paths_arr,
            'subset_id': subset_key_arr,
            'data_arr': data_arr,
            'true_fastas': true_fasta_dest,
            'contam_fastas': contam_fasta_dest,
            'rep_seq_indicator': rep_seq_indicator,
            'table_keys_arr': sorted_key_arr,
            'feat_or_read': feature_or_read_arr,
    })
    #sortable JS code 
    js = os.path.join(
        TEMPLATES, 'js', 'tsorter.min.js')
    os.mkdir(os.path.join(output_dir, 'js'))
    shutil.copy(js, os.path.join(output_dir, 'js', 'tsorter.min.js'))
