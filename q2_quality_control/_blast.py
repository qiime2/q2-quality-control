# ----------------------------------------------------------------------------
# Copyright (c) 2017, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import tempfile
import subprocess


def _search_seqs(query_sequences, reference_sequences, evalue,
                 perc_identity, threads, perc_query_aligned, method):
    if method == 'blast':
        # blast uses float format but vsearch uses int for perc_identity
        perc_identity = perc_identity * 100
        cmd = _blast(
            query_sequences, reference_sequences, evalue, perc_identity)
    elif method == 'blastn-short':
        # blast uses float format but vsearch uses int for perc_identity
        perc_identity = perc_identity * 100
        cmd = _blastn_short(
            query_sequences, reference_sequences, evalue, perc_identity)
    elif method == 'vsearch':
        cmd = _vsearch(
            query_sequences, reference_sequences, perc_identity, threads)
    return _generate_assignments(cmd, perc_query_aligned)


def _blast(query_sequences, reference_sequences, evalue, perc_identity):
    seqs_fp = str(query_sequences)
    ref_fp = str(reference_sequences)
    cmd = ['blastn', '-query', seqs_fp, '-strand',
           'both', '-outfmt', '6 qseqid sseqid qlen qstart qend', '-subject',
           ref_fp, '-perc_identity', str(perc_identity),
           '-max_target_seqs', '1']
    if evalue is not None:
        cmd.extend(['-evalue', str(evalue)])
    cmd.append('-out')
    return cmd


def _blastn_short(query_sequences, reference_sequences, evalue,
                  perc_identity):
    # Should have identical settings to blast, but adjust word size and filter
    cmd = _blast(query_sequences, reference_sequences, evalue, perc_identity)
    cmd = cmd[:-1] + ['-word_size', '7', '-dust', 'no', '-out']
    return cmd


def _vsearch(query_sequences, reference_sequences, perc_identity, threads):
    seqs_fp = str(query_sequences)
    ref_fp = str(reference_sequences)
    cmd = ['vsearch', '--usearch_global', seqs_fp, '--id', str(perc_identity),
           '--strand', 'both', '--maxaccepts', '1', '--maxrejects', '0',
           '--db', ref_fp, '--threads', str(threads),
           '--userfields', 'query+target+ql+qlo+qhi', '--userout']
    return cmd


def _generate_assignments(cmd, perc_query_aligned):
    '''Run command line subprocess and extract hits.'''
    with tempfile.NamedTemporaryFile() as output:
        cmd = cmd + [output.name]
        _run_command(cmd)
        hits = _extract_hits(output.name, perc_query_aligned)
        return hits


def _extract_hits(blast_output, perc_query_aligned):
    '''import observed assignments in blast6 or blast7 format, return list of
    query IDs receiving hits.
    blast_output: path or list
        Taxonomy observation map in blast format 6 or 7. Each line consists of
        taxonomy assignments of a query sequence in tab-delimited format:
            <query_id>    <assignment_id>   <...other columns are ignored>
    '''
    hits = set()
    with open(blast_output, "r") as inputfile:
        # grab query IDs from each line (only queries with hits are listed)
        for line in inputfile:
            # ignore comment lines and blank lines
            if not line.startswith('#') and line != "":
                query_id, subject_id, query_len, start, end = line.split('\t')
                # check how much of alignment covers query
                perc_coverage = _perc_coverage(end, start, query_len)
                # check for minimum perc_query_aligned
                # if vsearch fails to find assignment, it reports '*' as the
                # accession ID, so we will not count those IDs as hits.
                if perc_coverage >= perc_query_aligned and subject_id != '*':
                    hits.add(query_id)
    return hits


def _perc_coverage(end, start, query_len):
    # query start is one-based relative to start of sequence
    # and hence we add 1 to adjust for comparison vs. length.
    # E.g., alignment of two identical 10-nt seqs will yield:
    # start = 1, end = 15. Hence 15 - 1 + 1 = 15 nt full
    # length of alignment.
    return (float(end) - float(start) + 1) / float(query_len)


# Replace this function with QIIME2 API for wrapping commands/binaries,
# pending https://github.com/qiime2/qiime2/issues/224
def _run_command(cmd, verbose=True):
    if verbose:
        print("Running external command line application. This may print "
              "messages to stdout and/or stderr.")
        print("The command being run is below. This command cannot "
              "be manually re-run as it will depend on temporary files that "
              "no longer exist.")
        print("\nCommand:", end=' ')
        print(" ".join(cmd), end='\n\n')
    subprocess.run(cmd, check=True)
