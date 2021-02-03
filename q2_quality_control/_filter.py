# ----------------------------------------------------------------------------
# Copyright (c) 2017-2021, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import os
import shutil
import tempfile

from q2_types.feature_data import DNAFASTAFormat
from q2_types.bowtie2 import Bowtie2IndexDirFmt
from q2_types.per_sample_sequences import (
    CasavaOneEightSingleLanePerSampleDirFmt,
)

from ._utilities import _run_command


# samtools flags
# -f 4 keeps only single alignments that are unmapped
# -f 12 keeps only paired alignments with both reads unmapped
# -F 256 removes reads that are not primary alignment
# -F 260 removes reads that are not primary alignment or unmapped
# -F 268 removes reads that are not primary alignment or unmapped
# or pair is unmapped.
KEEP_UNMAPPED_SINGLE = '4'
KEEP_UNMAPPED_PAIRED = '12'
REMOVE_SECONDARY_ALIGNMENTS = '256'
REMOVE_SECONDARY_OR_UNMAPPED_SINGLE = '260'
REMOVE_SECONDARY_OR_UNMAPPED_PAIRED = '268'

_filter_defaults = {
    'n_threads': 1,
    'mode': 'local',
    'sensitivity': 'sensitive',
    'exclude_seqs': True,
    'ref_gap_open_penalty': 5,
    'ref_gap_ext_penalty': 3,
}


def bowtie2_build(sequences: DNAFASTAFormat,
                  n_threads: str = 1) -> Bowtie2IndexDirFmt:
    database = Bowtie2IndexDirFmt()
    build_cmd = ['bowtie2-build', '--threads', str(n_threads),
                 str(sequences), str(database.path / 'db')]
    _run_command(build_cmd)
    return database


def filter_reads(
        demultiplexed_sequences: CasavaOneEightSingleLanePerSampleDirFmt,
        database: Bowtie2IndexDirFmt,
        n_threads: int = _filter_defaults['n_threads'],
        mode: str = _filter_defaults['mode'],
        sensitivity: str = _filter_defaults['sensitivity'],
        ref_gap_open_penalty: str = _filter_defaults['ref_gap_open_penalty'],
        ref_gap_ext_penalty: str = _filter_defaults['ref_gap_ext_penalty'],
        exclude_seqs: str = _filter_defaults['exclude_seqs']) \
            -> CasavaOneEightSingleLanePerSampleDirFmt:
    filtered_seqs = CasavaOneEightSingleLanePerSampleDirFmt()
    df = demultiplexed_sequences.manifest
    fastq_paths = [record[1:] for record in df.itertuples()]

    for fwd, rev in fastq_paths:
        _bowtie2_filter(fwd, rev, filtered_seqs, database, n_threads, mode,
                        sensitivity, ref_gap_open_penalty, ref_gap_ext_penalty,
                        exclude_seqs)
    return filtered_seqs


def _bowtie2_filter(f_read, r_read, outdir, database, n_threads, mode,
                    sensitivity, ref_gap_open_penalty, ref_gap_ext_penalty,
                    exclude_seqs):
    if mode == 'local':
        mode = '--{0}-{1}'.format(sensitivity, mode)
    else:
        mode = '--' + sensitivity
    rfg_setting = '{0},{1}'.format(ref_gap_open_penalty, ref_gap_ext_penalty)

    with tempfile.NamedTemporaryFile() as sam_f:
        samfile_output_path = sam_f.name
        with tempfile.NamedTemporaryFile() as bam_f:
            bamfile_output_path = bam_f.name

            # align to reference with bowtie
            bowtie_cmd = ['bowtie2', '-p', str(n_threads), mode,
                          '--rfg', rfg_setting,
                          '-x', str(database.path / database.get_basename())]
            if r_read is not None:
                bowtie_cmd += ['-1', f_read, '-2', r_read]
            else:
                bowtie_cmd += ['-U', f_read]
            bowtie_cmd += ['-S', samfile_output_path]
            _run_command(bowtie_cmd)

            # Filter alignment and convert to BAM with samtools
            if exclude_seqs:
                sam_flags = ['-F', REMOVE_SECONDARY_ALIGNMENTS,
                             '-f', KEEP_UNMAPPED_SINGLE]
                if r_read is not None:
                    sam_flags[-1] = KEEP_UNMAPPED_PAIRED
            else:
                sam_flags = ['-F', REMOVE_SECONDARY_OR_UNMAPPED_SINGLE]
                if r_read is not None:
                    sam_flags[-1] = REMOVE_SECONDARY_OR_UNMAPPED_PAIRED
            samtools_command = ['samtools', 'view', '-b', samfile_output_path,
                                '-o', bamfile_output_path, *sam_flags,
                                '-@', str(n_threads - 1)]
            _run_command(samtools_command)
            # sort BAM file by read name so pairs are ordered
            if r_read is not None:
                with tempfile.NamedTemporaryFile() as sort_f:
                    bamfile_sorted_output_path = sort_f.name
                    sort_command = [
                        'samtools', 'sort', '-n', '-@', str(n_threads - 1),
                        '-o', bamfile_sorted_output_path, bamfile_output_path]
                    _run_command(sort_command)
                    shutil.copyfile(
                        bamfile_sorted_output_path, bamfile_output_path)

            # Convert to FASTQ with samtools
            fwd = str(outdir.path / os.path.basename(f_read))
            _reads = ['-1', fwd]
            if r_read is not None:
                rev = str(outdir.path / os.path.basename(r_read))
                _reads += ['-2', rev]
            # -s /dev/null excludes singletons
            # -0 /dev/null excludes supplementary and secondary reads
            # -n keeps samtools from altering header IDs!
            convert_command = [
                'samtools', 'fastq', *_reads, '-0', '/dev/null',
                '-s', '/dev/null', '-n', bamfile_output_path]
            _run_command(convert_command)
