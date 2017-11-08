# ----------------------------------------------------------------------------
# Copyright (c) 2017, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import tempfile
import pandas as pd

from ._utilities import _plot_histogram
from ._blast import _blast, _perc_coverage, _run_command


def _evaluate_seqs(query_sequences, reference_sequences, show_alignments):
    cmd = _blast(query_sequences, reference_sequences, None, 0)
    # modify command to output seq alignment information
    cmd[6] = ('6 qseqid sseqid pident length mismatch gapopen qstart qend '
              'sstart send evalue bitscore qlen qseq sseq')
    # generate report of blast alignment results
    results, alignments = _generate_alignments(cmd, show_alignments)

    # plot histogram of mismatches between obs seqs and top hit in exp seqs
    g = _plot_histogram(pd.to_numeric(results['Mismatches']))

    return results, alignments, g


def _generate_alignments(cmd, show_alignments):
    '''Run command line subprocess and extract hits.'''
    with tempfile.NamedTemporaryFile() as output:
        cmd = cmd + [output.name]
        _run_command(cmd)
        return _generate_alignment_results(output.name, show_alignments)


def _generate_alignment_results(blast_results, show_alignments):
    '''blast_results: output of blastn in outfmt == 6'''
    alignments = []
    results = []
    with open(blast_results, "r") as inputfile:
        for line in inputfile:
            line = line.rstrip()
            # qseqid sseqid pident length mismatch gapopen qstart qend
            # sstart send evalue bitscore qlen qseq sseq
            (qseqid, sseqid, pident, length, mismatch, gapopen, qstart, qend,
             sstart, send, evalue, bitscore, qlen, qseq, sseq) = \
                line.split('\t')
            # caluclate percent coverage based on qend, qstart, and qlen
            perc_coverage = _perc_coverage(qend, qstart, qlen)
            # add query and hit seqs to alignments
            if show_alignments:
                alignments.append((qseqid, 'query', qseq))
                alignments.append((qseqid, sseqid, sseq))
            # append blast results to report: remove seqs, add perc_coverage
            results.append((
                qseqid, sseqid, pident, length, mismatch, gapopen, qstart,
                qend, sstart, send, evalue, bitscore, qlen, perc_coverage))

    results = pd.DataFrame(results, columns=[
            'Query id', 'Subject id', 'Percent Identity', 'Alignment Length',
            'Mismatches', 'Gaps', 'Alignment Start (query)',
            'Alignment end (query)', 'Alignment Start (subject)',
            'Alignment end (subject)', 'E Value', 'Bit Score', 'Query length',
            'Percent Coverage'])

    if show_alignments:
        alignments = pd.DataFrame(
            alignments, columns=['Query id', 'Subject id', 'seq'])
        # convert to multiindex
        alignments = alignments.set_index(['Query id', 'Subject id'])
    else:
        alignments = None

    return results, alignments
