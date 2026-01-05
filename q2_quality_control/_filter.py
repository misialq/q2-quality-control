# ----------------------------------------------------------------------------
# Copyright (c) 2017-2026, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import os
import shutil
import tempfile
import gzip

from qiime2.plugin import get_available_cores
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


def bowtie2_build(
    sequences: DNAFASTAFormat, n_threads: int = 1
) -> Bowtie2IndexDirFmt:
    if n_threads == 0:
        n_threads = get_available_cores()

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
    exclude_seqs: str = _filter_defaults['exclude_seqs']
) -> (
        CasavaOneEightSingleLanePerSampleDirFmt,
        CasavaOneEightSingleLanePerSampleDirFmt,
        CasavaOneEightSingleLanePerSampleDirFmt):
    if n_threads == 0:
        n_threads = get_available_cores()

    filtered_seqs = CasavaOneEightSingleLanePerSampleDirFmt()
    other_seqs = CasavaOneEightSingleLanePerSampleDirFmt()
    flag0_seqs = CasavaOneEightSingleLanePerSampleDirFmt()
    df = demultiplexed_sequences.manifest
    fastq_paths = [record[1:] for record in df.itertuples()]

    for fwd, rev in fastq_paths:
        _bowtie2_filter(fwd, rev, filtered_seqs, other_seqs, flag0_seqs,
                        database, n_threads, mode, sensitivity,
                        ref_gap_open_penalty, ref_gap_ext_penalty,
                        exclude_seqs)
    return filtered_seqs, other_seqs, flag0_seqs


def _bowtie2_filter(f_read, r_read, outdir_keep, outdir_other, outdir_flag0,
                    database, n_threads, mode, sensitivity,
                    ref_gap_open_penalty, ref_gap_ext_penalty, exclude_seqs):
    if mode == 'local':
        mode = '--{0}-{1}'.format(sensitivity, mode)
    else:
        mode = '--' + sensitivity
    rfg_setting = '{0},{1}'.format(ref_gap_open_penalty, ref_gap_ext_penalty)

    with tempfile.NamedTemporaryFile() as sam_f:
        samfile_output_path = sam_f.name
        with tempfile.NamedTemporaryFile() as bam_keep, \
                tempfile.NamedTemporaryFile() as bam_other, \
                tempfile.NamedTemporaryFile() as flag0_keep_tmp, \
                tempfile.NamedTemporaryFile() as flag0_other_tmp:
            bam_keep_path = bam_keep.name
            bam_other_path = bam_other.name
            flag0_keep_path = flag0_keep_tmp.name
            flag0_other_path = flag0_other_tmp.name

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

            # Filter alignment with samtools to create separate keep and
            # complement BAM files. A second call to ``samtools view`` with the
            # reverse flag selection ensures the two outputs form a proper
            # complement of the original input reads.
            if exclude_seqs:
                keep_flags = ['-F', REMOVE_SECONDARY_ALIGNMENTS,
                              '-f', KEEP_UNMAPPED_SINGLE]
                if r_read is not None:
                    keep_flags[-1] = KEEP_UNMAPPED_PAIRED
                other_flags = ['-F', REMOVE_SECONDARY_OR_UNMAPPED_SINGLE]
                if r_read is not None:
                    other_flags[-1] = REMOVE_SECONDARY_OR_UNMAPPED_PAIRED
            else:
                keep_flags = ['-F', REMOVE_SECONDARY_OR_UNMAPPED_SINGLE]
                if r_read is not None:
                    keep_flags[-1] = REMOVE_SECONDARY_OR_UNMAPPED_PAIRED
                other_flags = ['-F', REMOVE_SECONDARY_ALIGNMENTS,
                               '-f', KEEP_UNMAPPED_SINGLE]
                if r_read is not None:
                    other_flags[-1] = KEEP_UNMAPPED_PAIRED

            samtools_keep = [
                'samtools', 'view', '-b', '-o', bam_keep_path,
                *keep_flags, '-@', str(n_threads - 1), samfile_output_path
            ]
            _run_command(samtools_keep)
            samtools_other = [
                'samtools', 'view', '-b', '-o', bam_other_path,
                *other_flags, '-@', str(n_threads - 1), samfile_output_path
            ]
            _run_command(samtools_other)
            # sort BAM files by read name so pairs are ordered
            if r_read is not None:
                with tempfile.NamedTemporaryFile() as sort_f1, \
                        tempfile.NamedTemporaryFile() as sort_f2:
                    sort_keep = [
                        'samtools', 'sort', '-n', '-@', str(n_threads - 1),
                        '-o', sort_f1.name, bam_keep_path]
                    _run_command(sort_keep)
                    shutil.copyfile(sort_f1.name, bam_keep_path)
                    sort_other = [
                        'samtools', 'sort', '-n', '-@', str(n_threads - 1),
                        '-o', sort_f2.name, bam_other_path]
                    _run_command(sort_other)
                    shutil.copyfile(sort_f2.name, bam_other_path)

            # Convert BAMs to FASTQ with samtools
            fwd_keep = str(outdir_keep.path / os.path.basename(f_read))
            fwd_other = str(outdir_other.path / os.path.basename(f_read))
            flag0_fp = str(outdir_flag0.path / os.path.basename(f_read))
            if r_read is None:
                convert_keep = [
                    'samtools', 'fastq', '-0', flag0_keep_path,
                    '-s', '/dev/null', '-@', str(n_threads - 1), '-n',
                    bam_keep_path]
                _run_command(convert_keep)
                convert_other = [
                    'samtools', 'fastq', '-0', flag0_other_path,
                    '-s', '/dev/null', '-@', str(n_threads - 1), '-n',
                    bam_other_path]
                _run_command(convert_other)
            else:
                rev_keep = str(outdir_keep.path / os.path.basename(r_read))
                rev_other = str(outdir_other.path / os.path.basename(r_read))
                convert_keep = [
                    'samtools', 'fastq', '-0', flag0_keep_path,
                    '-1', fwd_keep, '-2', rev_keep, '-s', '/dev/null',
                    '-@', str(n_threads - 1), '-n', bam_keep_path]
                _run_command(convert_keep)
                convert_other = [
                    'samtools', 'fastq', '-0', flag0_other_path,
                    '-1', fwd_other, '-2', rev_other, '-s', '/dev/null',
                    '-@', str(n_threads - 1), '-n', bam_other_path]
                _run_command(convert_other)

            # combine READ_OTHER outputs from both conversions
            with open(flag0_keep_path, 'rb') as fk, \
                    open(flag0_other_path, 'rb') as fo, \
                    gzip.open(flag0_fp, 'wb') as out_fh:
                shutil.copyfileobj(fk, out_fh)
                shutil.copyfileobj(fo, out_fh)
            # -s /dev/null excludes singletons
            # -n keeps samtools from altering header IDs!
