# ----------------------------------------------------------------------------
# Copyright (c) 2017-2025, QIIME 2 development team.
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


def _iter_ids(dirfmt):
    for _, obs_fp in dirfmt.sequences.iter_views(FastqGzFormat):
        with gzip.open(str(obs_fp), 'rt') as obs_fh:
            for records in itertools.zip_longest(*[obs_fh] * 4):
                if records[0] is None:
                    break
                obs_seq_h = records[0]
                obs_id = obs_seq_h.strip('@/012\n')
                yield obs_id


class TestFilterSingle(QualityControlTestsBase):

    def setUp(self):
        super().setUp()

        self.demuxed_art = Artifact.load(self.get_data_path('single-end.qza'))
        self.indexed_genome = Artifact.load(
            self.get_data_path('sars2-indexed.qza'))

    def test_filter_single_exclude_seqs(self):
        obs_art, comp_art, singleton_art = self.plugin.methods['filter_reads'](
            self.demuxed_art, self.indexed_genome, exclude_seqs=True)
        obs = obs_art.view(SingleLanePerSampleSingleEndFastqDirFmt)
        comp = comp_art.view(SingleLanePerSampleSingleEndFastqDirFmt)
        singleton = singleton_art.view(SingleLanePerSampleSingleEndFastqDirFmt)

        obs_ids = list(_iter_ids(obs))
        self.assertNotEqual(len(obs_ids), 0)
        for obs_id in obs_ids:
            # Make sure seqs that map to genome were removed
            self.assertTrue(obs_id not in seq_ids_that_map)
            self.assertTrue(obs_id in seq_id_that_does_not_map)

        comp_ids = list(_iter_ids(comp))
        self.assertNotEqual(len(comp_ids), 0)
        for comp_id in comp_ids:
            self.assertTrue(comp_id in seq_ids_that_map)
            self.assertTrue(comp_id not in seq_id_that_does_not_map)

        singleton_ids = list(_iter_ids(singleton))
        for singleton_id in singleton_ids:
            self.assertTrue(singleton_id not in obs_ids)

    def test_filter_single_keep_seqs(self):
        obs_art, comp_art, singleton_art = self.plugin.methods['filter_reads'](
            self.demuxed_art, self.indexed_genome, exclude_seqs=False)
        obs = obs_art.view(SingleLanePerSampleSingleEndFastqDirFmt)
        comp = comp_art.view(SingleLanePerSampleSingleEndFastqDirFmt)
        singleton = singleton_art.view(SingleLanePerSampleSingleEndFastqDirFmt)

        obs_ids = list(_iter_ids(obs))
        self.assertNotEqual(len(obs_ids), 0)
        for obs_id in obs_ids:
            # Make sure seqs that do not map to genome were removed
            self.assertTrue(obs_id in seq_ids_that_map)
            self.assertTrue(obs_id not in seq_id_that_does_not_map)

        comp_ids = list(_iter_ids(comp))
        self.assertNotEqual(len(comp_ids), 0)
        for comp_id in comp_ids:
            self.assertTrue(comp_id not in seq_ids_that_map)
            self.assertTrue(comp_id in seq_id_that_does_not_map)

        singleton_ids = list(_iter_ids(singleton))
        for singleton_id in singleton_ids:
            self.assertTrue(singleton_id not in obs_ids)


class TestFilterPaired(QualityControlTestsBase):

    def setUp(self):
        super().setUp()

        self.demuxed_art = Artifact.load(self.get_data_path('paired-end.qza'))
        self.indexed_genome = Artifact.load(
            self.get_data_path('sars2-indexed.qza'))

    def test_filter_paired_exclude_seqs(self):
        obs_art, comp_art, singleton_art = self.plugin.methods['filter_reads'](
            self.demuxed_art, self.indexed_genome, exclude_seqs=True)
        obs = obs_art.view(SingleLanePerSamplePairedEndFastqDirFmt)
        comp = comp_art.view(SingleLanePerSamplePairedEndFastqDirFmt)
        singleton = singleton_art.view(SingleLanePerSampleSingleEndFastqDirFmt)

        obs_ids = list(_iter_ids(obs))
        self.assertNotEqual(len(obs_ids), 0)
        for obs_id in obs_ids:
            # Make sure seqs that map to genome were removed
            self.assertTrue(obs_id not in seq_ids_that_map)
            self.assertTrue(obs_id in seq_id_that_does_not_map)

        comp_ids = list(_iter_ids(comp))
        self.assertNotEqual(len(comp_ids), 0)
        for comp_id in comp_ids:
            self.assertTrue(comp_id in seq_ids_that_map)
            self.assertTrue(comp_id not in seq_id_that_does_not_map)

        singleton_ids = list(_iter_ids(singleton))
        for singleton_id in singleton_ids:
            self.assertTrue(singleton_id not in obs_ids)

    def test_filter_paired_keep_seqs(self):
        obs_art, comp_art, singleton_art = self.plugin.methods['filter_reads'](
            self.demuxed_art, self.indexed_genome, exclude_seqs=False)
        obs = obs_art.view(SingleLanePerSamplePairedEndFastqDirFmt)
        comp = comp_art.view(SingleLanePerSamplePairedEndFastqDirFmt)
        singleton = singleton_art.view(SingleLanePerSampleSingleEndFastqDirFmt)

        obs_ids = list(_iter_ids(obs))
        self.assertNotEqual(len(obs_ids), 0)
        for obs_id in obs_ids:
            # Make sure seqs that do not map to genome were removed
            self.assertTrue(obs_id in seq_ids_that_map)
            self.assertTrue(obs_id not in seq_id_that_does_not_map)

        comp_ids = list(_iter_ids(comp))
        self.assertNotEqual(len(comp_ids), 0)
        for comp_id in comp_ids:
            self.assertTrue(comp_id not in seq_ids_that_map)
            self.assertTrue(comp_id in seq_id_that_does_not_map)

        singleton_ids = list(_iter_ids(singleton))
        for singleton_id in singleton_ids:
            self.assertTrue(singleton_id not in obs_ids)


if __name__ == '__main__':
    unittest.main()
