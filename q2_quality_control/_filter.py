# ----------------------------------------------------------------------------
# Copyright (c) 2017-2025, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import gzip
import os
import shutil
import tempfile

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
) -> (CasavaOneEightSingleLanePerSampleDirFmt,
      CasavaOneEightSingleLanePerSampleDirFmt,
      CasavaOneEightSingleLanePerSampleDirFmt):
    if n_threads == 0:
        n_threads = get_available_cores()

    filtered_seqs = CasavaOneEightSingleLanePerSampleDirFmt()
    complement_seqs = CasavaOneEightSingleLanePerSampleDirFmt()
    singleton_seqs = CasavaOneEightSingleLanePerSampleDirFmt()
    df = demultiplexed_sequences.manifest
    fastq_paths = [record[1:] for record in df.itertuples()]

    for fwd, rev in fastq_paths:
        _bowtie2_filter(
            fwd, rev, filtered_seqs, complement_seqs, singleton_seqs,
            database, n_threads, mode, sensitivity, ref_gap_open_penalty,
            ref_gap_ext_penalty, exclude_seqs)
    return filtered_seqs, complement_seqs, singleton_seqs


def _bowtie2_filter(f_read, r_read, filtered_outdir, complement_outdir,
                    singleton_outdir, database, n_threads, mode, sensitivity,
                    ref_gap_open_penalty, ref_gap_ext_penalty, exclude_seqs):
    if mode == 'local':
        mode = '--{0}-{1}'.format(sensitivity, mode)
    else:
        mode = '--' + sensitivity
    rfg_setting = '{0},{1}'.format(ref_gap_open_penalty, ref_gap_ext_penalty)

    with tempfile.NamedTemporaryFile() as sam_f:
        samfile_output_path = sam_f.name
        with tempfile.NamedTemporaryFile() as primary_bam_f:
            primary_bam_path = primary_bam_f.name

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

            # Remove secondary alignments then split into filtered/complement
            # so outputs are strict complements of each other.
            samtools_primary_cmd = [
                'samtools', 'view', '-b', '-F', REMOVE_SECONDARY_ALIGNMENTS,
                samfile_output_path, '-o', primary_bam_path,
                '-@', str(max(n_threads - 1, 0))]
            _run_command(samtools_primary_cmd)

            with tempfile.NamedTemporaryFile() as filtered_bam_f, \
                    tempfile.NamedTemporaryFile() as complement_bam_f:
                filtered_bam_path = filtered_bam_f.name
                complement_bam_path = complement_bam_f.name

                filter_flags = _get_filter_flags(exclude_seqs,
                                                 r_read is not None)
                samtools_filter_cmd = [
                    'samtools', 'view', '-b', primary_bam_path,
                    '-o', filtered_bam_path, '-U', complement_bam_path,
                    *filter_flags, '-@', str(max(n_threads - 1, 0))]
                _run_command(samtools_filter_cmd)

                # sort BAM files by read name so pairs are ordered
                if r_read is not None:
                    _sort_bam_by_name(filtered_bam_path, n_threads)
                    _sort_bam_by_name(complement_bam_path, n_threads)

                # Convert filtered/complement to FASTQ with samtools
                filtered_fwd = str(filtered_outdir.path /
                                   os.path.basename(f_read))
                complement_fwd = str(complement_outdir.path /
                                     os.path.basename(f_read))
                singleton_path = str(singleton_outdir.path /
                                     os.path.basename(f_read))

                if r_read is None:
                    _write_single_end_fastq(
                        filtered_bam_path, filtered_fwd, n_threads)
                    _write_single_end_fastq(
                        complement_bam_path, complement_fwd, n_threads)
                    _ensure_empty_fastq(singleton_path)
                    _log_counts_single_end(
                        f_read, filtered_fwd, complement_fwd, singleton_path)
                    return

                filtered_rev = str(filtered_outdir.path /
                                   os.path.basename(r_read))
                complement_rev = str(complement_outdir.path /
                                     os.path.basename(r_read))

                filtered_singleton_paths = _singleton_temp_paths(
                    singleton_path)
                complement_singleton_paths = _singleton_temp_paths(
                    singleton_path)

                _write_paired_end_fastq(
                    filtered_bam_path, filtered_fwd, filtered_rev,
                    filtered_singleton_paths, n_threads)
                _write_paired_end_fastq(
                    complement_bam_path, complement_fwd, complement_rev,
                    complement_singleton_paths, n_threads)

                _concat_files(
                    filtered_singleton_paths + complement_singleton_paths,
                    singleton_path)
                _cleanup_files(
                    filtered_singleton_paths + complement_singleton_paths)

                _log_counts_paired(
                    f_read, filtered_fwd, filtered_rev, complement_fwd,
                    complement_rev, singleton_path)


def _get_filter_flags(exclude_seqs, is_paired):
    if exclude_seqs:
        if is_paired:
            return ['-f', KEEP_UNMAPPED_PAIRED]
        return ['-f', KEEP_UNMAPPED_SINGLE]
    if is_paired:
        return ['-F', KEEP_UNMAPPED_PAIRED]
    return ['-F', KEEP_UNMAPPED_SINGLE]


def _sort_bam_by_name(bam_path, n_threads):
    with tempfile.NamedTemporaryFile() as sort_f:
        bamfile_sorted_output_path = sort_f.name
        sort_command = [
            'samtools', 'sort', '-n', '-@', str(max(n_threads - 1, 0)),
            '-o', bamfile_sorted_output_path, bam_path]
        _run_command(sort_command)
        shutil.copyfile(bamfile_sorted_output_path, bam_path)


def _write_single_end_fastq(bam_path, output_path, n_threads):
    convert_command = [
        'samtools', 'fastq', '-0', output_path, '-@',
        str(max(n_threads - 1, 0)), '-n', bam_path]
    _run_command(convert_command)


def _write_paired_end_fastq(bam_path, fwd_path, rev_path, singleton_paths,
                            n_threads):
    # -n keeps samtools from altering header IDs!
    # -0 captures reads not part of a pair; -s captures orphaned mates.
    convert_command = [
        'samtools', 'fastq', '-0', singleton_paths[0],
        '-1', fwd_path, '-2', rev_path, '-s', singleton_paths[1],
        '-@', str(max(n_threads - 1, 0)), '-n', bam_path]
    _run_command(convert_command)


def _singleton_temp_paths(output_path):
    suffix = ''.join(os.path.basename(output_path).split('.', 1)[1:])
    if suffix:
        suffix = f".{suffix}"
    else:
        suffix = '.fastq'
    temp_paths = []
    for _ in range(2):
        temp_fh = tempfile.NamedTemporaryFile(suffix=suffix, delete=False)
        temp_paths.append(temp_fh.name)
        temp_fh.close()
    return temp_paths


def _concat_files(paths, dest_path):
    wrote_data = False
    with open(dest_path, 'wb') as out_fh:
        for path in paths:
            if not os.path.exists(path):
                continue
            wrote_data = True
            with open(path, 'rb') as in_fh:
                shutil.copyfileobj(in_fh, out_fh)
    if not wrote_data and str(dest_path).endswith('.gz'):
        with gzip.open(dest_path, 'wb'):
            pass


def _cleanup_files(paths):
    for path in paths:
        try:
            os.remove(path)
        except FileNotFoundError:
            pass


def _ensure_empty_fastq(path):
    if str(path).endswith('.gz'):
        with gzip.open(path, 'wb'):
            pass
    else:
        open(path, 'wb').close()


def _count_fastq_reads(path):
    if not os.path.exists(path):
        return 0
    opener = gzip.open if str(path).endswith('.gz') else open
    with opener(path, 'rt') as fh:
        return sum(1 for _ in fh) // 4


def _log_counts_single_end(f_read, filtered_path, complement_path,
                           singleton_path):
    print('filter_reads output counts for %s' % os.path.basename(f_read))
    print('filtered_sequences: %d reads -> %s' % (
        _count_fastq_reads(filtered_path), filtered_path))
    print('complement_sequences: %d reads -> %s' % (
        _count_fastq_reads(complement_path), complement_path))
    print('singleton_sequences: %d reads -> %s' % (
        _count_fastq_reads(singleton_path), singleton_path))


def _log_counts_paired(f_read, filtered_fwd, filtered_rev, complement_fwd,
                       complement_rev, singleton_path):
    print('filter_reads output counts for %s' % os.path.basename(f_read))
    print('filtered_sequences forward: %d reads -> %s' % (
        _count_fastq_reads(filtered_fwd), filtered_fwd))
    print('filtered_sequences reverse: %d reads -> %s' % (
        _count_fastq_reads(filtered_rev), filtered_rev))
    print('complement_sequences forward: %d reads -> %s' % (
        _count_fastq_reads(complement_fwd), complement_fwd))
    print('complement_sequences reverse: %d reads -> %s' % (
        _count_fastq_reads(complement_rev), complement_rev))
    print('singleton_sequences: %d reads -> %s' % (
        _count_fastq_reads(singleton_path), singleton_path))
