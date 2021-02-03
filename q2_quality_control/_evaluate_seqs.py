# ----------------------------------------------------------------------------
# Copyright (c) 2017-2021, QIIME 2 development team.
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
    cmd[6] = '6 qseqid sseqid pident length mismatch gapopen qstart qend ' \
             'sstart send evalue bitscore qlen qseq sseq'
    # generate report of blast alignment results
    results, alignments = _generate_alignments(cmd, show_alignments)

    # plot histogram of mismatches between obs seqs and top hit in exp seqs
    g = _plot_histogram(pd.to_numeric(results['Mismatches']))

    return results, alignments, g


def _generate_alignments(cmd, show_alignments):
    '''Run command line subprocess and extract hits.'''
    with tempfile.NamedTemporaryFile() as output:
        cmd.append(output.name)
        _run_command(cmd)
        return _generate_alignment_results(output.name, show_alignments)


def _generate_alignment_results(blast_results, show_alignments):
    '''blast_results: output of blastn in outfmt == 6'''
    results = pd.read_csv(blast_results, sep='\t', names=[
        'Query id', 'Subject id', 'Percent Identity', 'Alignment Length',
        'Mismatches', 'Gaps', 'Alignment Start (query)',
        'Alignment end (query)', 'Alignment Start (subject)',
        'Alignment end (subject)', 'E Value', 'Bit Score', 'Query length',
        'qseq', 'sseq'])
    # caluclate percent coverage based on qend, qstart, and qlen
    results['Percent Coverage'] = results.apply(lambda x: _perc_coverage(
        x['Alignment end (query)'], x['Alignment Start (query)'],
        x['Query length']), axis=1)

    if show_alignments:
        alignments = []
        for i, d in results.iterrows():
            alignments.append((d['Query id'], 'query', d['qseq']))
            alignments.append((d['Query id'], d['Subject id'], d['sseq']))
        alignments = pd.DataFrame(
            alignments, columns=['Query id', 'Subject id', 'seq'])
        # convert to multiindex
        alignments = alignments.set_index(['Query id', 'Subject id'])
    else:
        alignments = None

    results = results.drop(['qseq', 'sseq'], axis=1)

    return results, alignments
