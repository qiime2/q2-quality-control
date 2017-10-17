# ----------------------------------------------------------------------------
# Copyright (c) 2017, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import tempfile
import subprocess
import pandas as pd


def _search_seqs(feature_sequences, reference_sequences, evalue,
                 perc_identity, threads, method):
    if method == 'blast':
        # blast uses float format but vsearch uses int for perc_identity
        perc_identity = perc_identity * 100
        cmd = _blast(
            feature_sequences, reference_sequences, evalue, perc_identity)
    elif method == 'vsearch':
        cmd = _vsearch(
            feature_sequences, reference_sequences, perc_identity, threads)
    return _generate_assignments(cmd)


def _blast(feature_sequences, reference_sequences, evalue, perc_identity):
    seqs_fp = str(feature_sequences)
    ref_fp = str(reference_sequences)
    cmd = ['blastn', '-query', seqs_fp, '-evalue', str(evalue), '-strand',
           'both', '-outfmt', '7', '-subject', ref_fp, '-perc_identity',
           str(perc_identity), '-max_target_seqs', '1', '-out']
    return cmd


def _vsearch(feature_sequences, reference_sequences, perc_identity, threads):
    seqs_fp = str(feature_sequences)
    ref_fp = str(reference_sequences)
    cmd = ['vsearch', '--usearch_global', seqs_fp, '--id', str(perc_identity),
           '--strand', 'both', '--maxaccepts', '1', '--maxrejects', '0',
           '--output_no_hits', '--db', ref_fp, '--threads', str(threads),
           '--blast6out']
    return cmd


def _generate_assignments(cmd):
    '''Run command line subprocess and extract hits.'''
    with tempfile.NamedTemporaryFile() as output:
        cmd = cmd + [output.name]
        _run_command(cmd)
        hits = _extract_hits(output.name)
        result = pd.DataFrame(hits, index=hits, columns=['Feature ID'])
        result.index.name = 'Feature ID'
        return result


def _extract_hits(blast_output):
    '''import observed assignments in blast6 or blast7 format, return list of
    query IDs receiving hits.
    blast_output: path or list
        Taxonomy observation map in blast format 6 or 7. Each line consists of
        taxonomy assignments of a query sequence in tab-delimited format:
            <query_id>    <assignment_id>   <...other columns are ignored>
    '''
    with open(blast_output, "r") as inputfile:
        # grab query IDs from each line (only queries with hits are listed)
        hits = {line.split('\t')[0] for line in inputfile
                # ignore comment lines and blank lines
                if not line.startswith('#')
                and line != ""
                # if vsearch fails to find assignment, it reports '*' as the
                # accession ID, so we will not count those IDs as hits.
                and line.split('\t')[1] != '*'}
    return list(hits)


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
