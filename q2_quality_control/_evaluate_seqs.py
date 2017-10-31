# ----------------------------------------------------------------------------
# Copyright (c) 2017, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import tempfile
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import pkg_resources
import itertools

from ._utilities import _plot_histogram
from ._blast import _blast, _perc_coverage, _run_command


def _evaluate_seqs(query_sequences, reference_sequences):
    cmd = _blast(query_sequences, reference_sequences, None, 0)
    # modify command to output seq alignment information
    cmd[6] = ('6 qseqid sseqid pident length mismatch gapopen qstart qend '
              'sstart send evalue bitscore qlen qseq sseq btop')
    # generate report of blast alignment results
    alignments = _generate_alignments(cmd)

    # plot histogram of mismatches between obs seqs and top hit in exp seqs
    g = _plot_histogram(pd.to_numeric(alignments['Mismatches']))


def _generate_alignments(cmd):
    '''Run command line subprocess and extract hits.'''
    with tempfile.NamedTemporaryFile() as output:
        cmd = cmd + [output.name]
        _run_command(cmd)
        return _generate_alignment_results(output.name)


def _generate_alignment_results(blast_results):
    '''blast_results: output of blastn in outfmt == 6'''
    alignments = []
    for line in blast_results:
        q = line.split('\t')
        # parse blast traceback (btop) alignment results
        btop = _parse_blast_traceback(_split_numeric(q[-1]))
        # caluclate percent coverage based on qend, qstart, and qlen
        perc_coverage = _perc_coverage(q[7], q[6], q[-4])
        # join alignment information into a single linebreak-delimited string
        aln = "\n".join([q[-3], btop, q[-2]])
        # append blast results to report: remove btop, add perc_coverage
        alignments.append(tuple(q[:-3] + perc_coverage + aln))
    alignments = pd.DataFrame(
        alignments, columns=[
            'Query id', 'Subject id', 'Percent Identity', 'Alignment Length',
            'Mismatches', 'Gaps', 'Alignment Start (query)',
            'Alignment end (query)', 'Alignment Start (subject)',
            'Alignment end (subject)', 'E Value', 'Bit Score', 'Query length',
            'Percent Coverage', 'Pairwise Alignment'])
    alignments.set_index(['Query id'])
    return alignments


def _split_numeric(s):
    '''Split string into list of substrings grouped by numeric/non-numeric
    values
    '''
    return ["".join(x) for _, x in itertools.groupby(s, key=str.isdigit)]


def _parse_blast_traceback(btop):
    results = []
    for s in btop:
        # numbers indicate length of match
        if s.isdigit():
            results.append('|' * int(s))
        # characters indicate base mismatches. E.g., AG indicates A on query,
        # G on subject at next position; A- indicates A on query, gap in
        # subject.
        else:
            results.append(' ' * int(len(s) / 2))
    return "".join(results)
