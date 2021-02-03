# ----------------------------------------------------------------------------
# Copyright (c) 2017-2021, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import gzip
import itertools
import unittest

from q2_types.per_sample_sequences import (
    SingleLanePerSampleSingleEndFastqDirFmt,
    SingleLanePerSamplePairedEndFastqDirFmt,
    FastqGzFormat,
)
from qiime2 import Artifact
from .test_quality_control import QualityControlTestsBase


class TestBowtie2Build(QualityControlTestsBase):

    # This test just makes sure this runs without error, which will include
    # validating the contents.
    def test_build(self):
        genome = Artifact.load(self.get_data_path('sars2-genome.qza'))
        self.plugin.methods['bowtie2_build'](genome)


seq_ids_that_map = ['SARS2:6:73:941:1973#', 'SARS2:6:73:231:3321#',
                    'SARS2:6:73:233:3421#', 'SARS2:6:73:552:2457#',
                    'SARS2:6:73:567:7631#']
seq_id_that_does_not_map = 'SARS2:6:73:356:9806#'


class TestFilterSingle(QualityControlTestsBase):

    def setUp(self):
        super().setUp()

        self.demuxed_art = Artifact.load(self.get_data_path('single-end.qza'))
        self.indexed_genome = Artifact.load(
            self.get_data_path('sars2-indexed.qza'))

    def test_filter_single_exclude_seqs(self):
        obs_art, = self.plugin.methods['filter_reads'](
            self.demuxed_art, self.indexed_genome, exclude_seqs=True)
        obs = obs_art.view(SingleLanePerSampleSingleEndFastqDirFmt)
        obs_seqs = obs.sequences.iter_views(FastqGzFormat)
        for _, obs_fp in obs_seqs:
            with gzip.open(str(obs_fp), 'rt') as obs_fh:
                # Iterate over expected and observed reads, side-by-side
                for records in itertools.zip_longest(*[obs_fh] * 4):
                    (obs_seq_h, obs_seq, _, obs_qual) = records
                    # Make sure seqs that map to genome were removed
                    obs_id = obs_seq_h.strip('@/012\n')
                    self.assertTrue(obs_id not in seq_ids_that_map)
                    self.assertTrue(obs_id in seq_id_that_does_not_map)

    def test_filter_single_keep_seqs(self):
        obs_art, = self.plugin.methods['filter_reads'](
            self.demuxed_art, self.indexed_genome, exclude_seqs=False)
        obs = obs_art.view(SingleLanePerSampleSingleEndFastqDirFmt)
        obs_seqs = obs.sequences.iter_views(FastqGzFormat)
        for _, obs_fp in obs_seqs:
            with gzip.open(str(obs_fp), 'rt') as obs_fh:
                # Iterate over expected and observed reads, side-by-side
                for records in itertools.zip_longest(*[obs_fh] * 4):
                    (obs_seq_h, obs_seq, _, obs_qual) = records
                    # Make sure seqs that do not map to genome were removed
                    obs_id = obs_seq_h.strip('@/012\n')
                    self.assertTrue(obs_id in seq_ids_that_map)
                    self.assertTrue(obs_id not in seq_id_that_does_not_map)


class TestFilterPaired(QualityControlTestsBase):

    def setUp(self):
        super().setUp()

        self.demuxed_art = Artifact.load(self.get_data_path('paired-end.qza'))
        self.indexed_genome = Artifact.load(
            self.get_data_path('sars2-indexed.qza'))

    def test_filter_paired_exclude_seqs(self):
        obs_art, = self.plugin.methods['filter_reads'](
            self.demuxed_art, self.indexed_genome, exclude_seqs=True)
        obs = obs_art.view(SingleLanePerSamplePairedEndFastqDirFmt)
        obs_seqs = obs.sequences.iter_views(FastqGzFormat)
        for _, obs_fp in obs_seqs:
            with gzip.open(str(obs_fp), 'rt') as obs_fh:
                # Iterate over expected and observed reads, side-by-side
                for records in itertools.zip_longest(*[obs_fh] * 4):
                    (obs_seq_h, obs_seq, _, obs_qual) = records
                    # Make sure seqs that map to genome were removed
                    obs_id = obs_seq_h.strip('@/012\n')
                    self.assertTrue(obs_id not in seq_ids_that_map)
                    self.assertTrue(obs_id in seq_id_that_does_not_map)

    def test_filter_paired_keep_seqs(self):
        obs_art, = self.plugin.methods['filter_reads'](
            self.demuxed_art, self.indexed_genome, exclude_seqs=False)
        obs = obs_art.view(SingleLanePerSamplePairedEndFastqDirFmt)
        obs_seqs = obs.sequences.iter_views(FastqGzFormat)
        for _, obs_fp in obs_seqs:
            with gzip.open(str(obs_fp), 'rt') as obs_fh:
                # Iterate over expected and observed reads, side-by-side
                for records in itertools.zip_longest(*[obs_fh] * 4):
                    (obs_seq_h, obs_seq, _, obs_qual) = records
                    # Make sure seqs that do not map to genome were removed
                    obs_id = obs_seq_h.strip('@/012\n')
                    self.assertTrue(obs_id in seq_ids_that_map)
                    self.assertTrue(obs_id not in seq_id_that_does_not_map)


if __name__ == '__main__':
    unittest.main()
