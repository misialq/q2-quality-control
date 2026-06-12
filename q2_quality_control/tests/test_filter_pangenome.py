# ----------------------------------------------------------------------------
# Copyright (c) 2017-2026, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import filecmp
import gzip
import io
import itertools
import os
import shutil
import tempfile
import unittest
import zipfile
from unittest.mock import patch, ANY, call, MagicMock

import requests
from qiime2 import Artifact
from q2_types.per_sample_sequences import (
    SingleLanePerSampleSingleEndFastqDirFmt,
    FastqGzFormat,
)

from q2_quality_control._filter_pangenome import (
    _fetch_and_extract_grch38,
    _verify_md5,
    _extract_fasta_from_gfa,
    _fetch_and_extract_pangenome,
    filter_reads_human_pangenome,
    _combine_fasta_files,
    EBI_SERVER_URL,
    NCBI_DATASETS_URL,
    GRCH38_GENOME_REL_FP,
    construct_human_pangenome_index,
)
from q2_quality_control._utilities import _run_command as real_run
from qiime2.plugin.testing import TestPluginBase

from qiime2.plugins.quality_control.pipelines import (
    construct_human_pangenome_index as construct_action,
    filter_reads_human_pangenome as filter_action
)


class TestPangenomeFiltering(TestPluginBase):
    package = "q2_quality_control.tests"

    @patch("q2_quality_control._filter_pangenome.requests.get")
    def test_fetch_and_extract_grch38(self, mock_get):
        dest_dir = self.temp_dir.name

        # zip the on-disk fixture that mirrors the NCBI Datasets download
        dataset_dir = self.get_data_path("pangenome/grch38_dataset")
        buffer = io.BytesIO()
        with zipfile.ZipFile(buffer, "w") as zf:
            for root, _, files in os.walk(dataset_dir):
                for name in files:
                    abs_fp = os.path.join(root, name)
                    zf.write(abs_fp, os.path.relpath(abs_fp, dataset_dir))
        zip_bytes = buffer.getvalue()

        # the mocked streaming HTTP response (used as a context manager)
        mock_response = MagicMock()
        mock_response.url = "https://api.ncbi.nlm.nih.gov/datasets/download"
        mock_response.iter_content.return_value = [zip_bytes]
        mock_get.return_value.__enter__.return_value = mock_response

        _fetch_and_extract_grch38(dest_dir)

        # the NCBI Datasets API is queried with the expected URL and params
        mock_get.assert_called_once_with(
            NCBI_DATASETS_URL,
            params={
                "include_annotation_type": ["GENOME_FASTA"],
                "hydrated": "FULLY_HYDRATED",
            },
            stream=True,
        )
        mock_response.raise_for_status.assert_called_once_with()

        # the verified genome is renamed to grch38.fasta with its contents
        # intact, and the original extracted file is gone (it was moved)
        grch38_fp = os.path.join(dest_dir, "grch38.fasta")
        self.assertTrue(os.path.exists(grch38_fp))
        self.assertTrue(
            filecmp.cmp(
                os.path.join(dataset_dir, GRCH38_GENOME_REL_FP), grch38_fp,
                shallow=False,
            ),
            "Extracted genome does not match the expected fixture",
        )
        self.assertFalse(
            os.path.exists(os.path.join(dest_dir, GRCH38_GENOME_REL_FP))
        )

    @patch("q2_quality_control._filter_pangenome._run_command")
    @patch("q2_quality_control._filter_pangenome.requests.get")
    def test_fetch_and_extract_grch38_download_error(self, mock_get, mock_run):
        mock_response = MagicMock()
        mock_response.iter_content.side_effect = (
            requests.exceptions.ChunkedEncodingError("connection dropped")
        )
        mock_get.return_value.__enter__.return_value = mock_response

        with self.assertRaisesRegex(Exception, "The download failed"):
            _fetch_and_extract_grch38(self.temp_dir.name)

        # extraction should never be reached if the download fails
        mock_run.assert_not_called()

    @patch("q2_quality_control._filter_pangenome._run_command")
    @patch("q2_quality_control._filter_pangenome.requests.get")
    def test_fetch_and_extract_grch38_http_error(self, mock_get, mock_run):
        mock_response = MagicMock()
        mock_response.raise_for_status.side_effect = (
            requests.exceptions.HTTPError("500 Server Error")
        )
        mock_get.return_value.__enter__.return_value = mock_response

        with self.assertRaisesRegex(Exception, "The download failed"):
            _fetch_and_extract_grch38(self.temp_dir.name)

        mock_run.assert_not_called()

    def test_verify_md5_valid(self):
        dataset_dir = self.get_data_path("pangenome/grch38_dataset")
        genome_fp = os.path.join(dataset_dir, GRCH38_GENOME_REL_FP)
        checksum_fp = os.path.join(dataset_dir, "md5sum.txt")

        _verify_md5(genome_fp, checksum_fp, GRCH38_GENOME_REL_FP)

    def test_verify_md5_mismatch(self):
        dataset_dir = self.get_data_path("pangenome/grch38_dataset")
        checksum_fp = os.path.join(dataset_dir, "md5sum.txt")
        # verify a different file against the genome's expected checksum
        other_fp = self.get_data_path("pangenome/pangenome.fasta")

        with self.assertRaisesRegex(
            ValueError, "does not match the expected MD5 hash"
        ):
            _verify_md5(other_fp, checksum_fp, GRCH38_GENOME_REL_FP)

    def test_verify_md5_missing_key(self):
        dataset_dir = self.get_data_path("pangenome/grch38_dataset")
        genome_fp = os.path.join(dataset_dir, GRCH38_GENOME_REL_FP)
        checksum_fp = os.path.join(dataset_dir, "md5sum.txt")

        with self.assertRaisesRegex(
            ValueError, "No checksum found for missing.fna"
        ):
            _verify_md5(genome_fp, checksum_fp, "missing.fna")

    @patch("subprocess.run")
    @patch("os.remove")
    def test_extract_from_gfa(self, p1, p2):
        fasta_fp = os.path.join(self.temp_dir.name, "some_fasta.fa")
        _extract_fasta_from_gfa("/some/gfa", fasta_fp)

        p2.assert_called_once_with(
            ["gfatools", "gfa2fa", "/some/gfa"], stdout=ANY
        )
        p1.assert_called_once_with("/some/gfa")

    @patch("subprocess.run", side_effect=OSError)
    def test_extract_from_gfa_error(self, p1):
        fasta_fp = os.path.join(self.temp_dir.name, "some_fasta.fa")
        with self.assertRaisesRegex(Exception, "Failed to extract"):
            _extract_fasta_from_gfa("/some/gfa", fasta_fp)

    def test_extract_from_gfa_end_to_end(self):
        # _extract_fasta_from_gfa deletes the input, so
        # work on a copy in the temp dir.
        gfa_fp = os.path.join(self.temp_dir.name, "pangenome.gfa")
        shutil.copy(self.get_data_path("pangenome/pangenome.gfa"), gfa_fp)
        obs = os.path.join(self.temp_dir.name, "from_gfa.fasta")

        _extract_fasta_from_gfa(gfa_fp, obs)

        self.assertTrue(
            filecmp.cmp(
                self.get_data_path("pangenome/pangenome.fasta"), obs,
                shallow=False,
            ),
            "Extracted FASTA does not match the expected output",
        )
        # the source GFA should be removed after a successful conversion
        self.assertFalse(os.path.exists(gfa_fp))

    @patch("q2_quality_control._filter_pangenome._run_command")
    def test_fetch_and_extract_pangenome(self, p1):
        _fetch_and_extract_pangenome("/some/where")

        dest_fp = os.path.join("/some/where", os.path.basename(EBI_SERVER_URL))
        p1.assert_has_calls(
            [
                call(["wget", EBI_SERVER_URL, "-q", "-O", dest_fp]),
                call(["gunzip", dest_fp]),
            ]
        )

    @patch(
        "q2_quality_control._filter_pangenome._run_command",
        side_effect=OSError
    )
    def test_fetch_and_extract_pangenome_error(self, p1):
        with self.assertRaisesRegex(Exception, "Unable to connect"):
            _fetch_and_extract_pangenome("/some/where")

    def test_combine_fasta_files_single(self):
        fname1 = "grch38"
        file1 = os.path.join(self.temp_dir.name, f"{fname1}.fasta")
        shutil.copy(self.get_data_path(f"pangenome/{fname1}.fasta"), file1)
        obs = os.path.join(self.temp_dir.name, "out.fasta")

        _combine_fasta_files(file1, fasta_out_fp=obs)

        self.assertTrue(
            filecmp.cmp(self.get_data_path(f"pangenome/{fname1}.fasta"), obs),
            "Files are not identical",
        )

    def test_combine_fasta_files_multi(self):
        fname1, fname2 = "pangenome", "grch38"
        file1 = os.path.join(self.temp_dir.name, f"{fname1}.fasta")
        file2 = os.path.join(self.temp_dir.name, f"{fname2}.fasta")
        shutil.copy(self.get_data_path(f"pangenome/{fname1}.fasta"), file1)
        shutil.copy(self.get_data_path(f"pangenome/{fname2}.fasta"), file2)
        obs = os.path.join(self.temp_dir.name, "out.fasta")

        _combine_fasta_files(file1, file2, fasta_out_fp=obs)

        self.assertTrue(
            filecmp.cmp(self.get_data_path("pangenome/combined.fasta"), obs),
            "Files are not identical",
        )

    @patch("subprocess.run", side_effect=OSError)
    def test_combine_fasta_files_error(self, p1):
        obs = os.path.join(self.temp_dir.name, "out.fasta")

        with self.assertRaisesRegex(Exception, "Failed to add the /fake/file"):
            _combine_fasta_files("/fake/file", fasta_out_fp=obs)

    @patch("q2_quality_control._filter_pangenome._fetch_and_extract_pangenome")
    @patch("q2_quality_control._filter_pangenome._fetch_and_extract_grch38")
    @patch("q2_quality_control._filter_pangenome._extract_fasta_from_gfa")
    def test_construct_pangenome_index(
        self, mock_extract_fasta, mock_fetch_grch38, mock_fetch_pangenome
    ):
        # we don't use the temp_dir from the test class as its content
        # would get deleted within the context that is being tested -
        # we need to be able to read files from that directory after
        # the function being tested exits
        temp_dir = tempfile.mkdtemp()

        # we construct our own context so that we can control each "action"
        # being retrieved from it
        ctx = MagicMock()
        ctx.get_action.return_value = MagicMock()
        mock_build_index_result = MagicMock()
        ctx.get_action("quality_control", "bowtie2_build").return_value = (
            mock_build_index_result,
        )
        ctx.make_artifact.return_value = MagicMock()

        # prepare some files that will be used by _combined_fasta_files
        open(os.path.join(temp_dir, "pangenome.gfa"), "w").close()
        shutil.copy(
            self.get_data_path("pangenome/grch38.fasta"),
            os.path.join(temp_dir, "grch38.fasta"),
        )
        shutil.copy(
            self.get_data_path("pangenome/pangenome.fasta"),
            os.path.join(temp_dir, "pangenome.fasta"),
        )

        with patch(
            "tempfile.TemporaryDirectory",
            return_value=MagicMock(
                name="TemporaryDirectory",
                __enter__=lambda x: temp_dir,
                __exit__=lambda x, y, z, w: None,
            ),
        ):
            generated_index = construct_human_pangenome_index(
                ctx=ctx, threads=1
            )

            # Assertions
            ctx.get_action.assert_any_call("quality_control", "bowtie2_build")

            mock_fetch_pangenome.assert_called_once_with(temp_dir)
            mock_fetch_grch38.assert_called_once_with(temp_dir)
            mock_extract_fasta.assert_called_once_with(
                os.path.join(temp_dir, "pangenome.gfa"),
                os.path.join(temp_dir, "pangenome.fasta"),
            )

            self.assertTrue(
                filecmp.cmp(
                    self.get_data_path("pangenome/combined.fasta"),
                    os.path.join(temp_dir, "combined.fasta"),
                ),
                "Files are not identical",
            )

            ctx.make_artifact.assert_called_once_with(
                "FeatureData[Sequence]",
                os.path.join(temp_dir, "combined.fasta")
            )
            ctx.get_action(
                "quality_control", "bowtie2_build"
            ).assert_has_calls(
                [call(sequences=ANY, n_threads=1)], any_order=True
            )

            self.assertIsNotNone(generated_index)

        # clean up
        shutil.rmtree(temp_dir)

    def test_construct_pangenome_index_end_to_end(self):
        # Runs the whole pipeline through the QIIME 2 framework with a real
        # ctx: gunzip, unzip, md5 verification, gfatools, seqtk, make_artifact
        # and a real bowtie2_build all execute. Only the two network downloads
        # (the pangenome `wget` and the NCBI `requests.get`) are stubbed with
        # local fixtures.

        # a real NCBI-style zip (genome + md5sum.txt) for the GRCh38 download
        dataset_dir = self.get_data_path("pangenome/grch38_dataset")
        buffer = io.BytesIO()
        with zipfile.ZipFile(buffer, "w") as zf:
            for root, _, files in os.walk(dataset_dir):
                for name in files:
                    abs_fp = os.path.join(root, name)
                    zf.write(abs_fp, os.path.relpath(abs_fp, dataset_dir))
        grch38_zip = buffer.getvalue()

        gfa_fixture = self.get_data_path("pangenome/pangenome.gfa")

        def fake_run(cmd, **kwargs):
            # the only network call routed through _run_command is the
            # pangenome `wget`; stage a local gzipped GFA in its place and
            # let everything else (gunzip, unzip) run for real
            if cmd[0] == "wget":
                dest_fp = cmd[cmd.index("-O") + 1]
                with open(gfa_fixture, "rb") as src, \
                        gzip.open(dest_fp, "wb") as dst:
                    shutil.copyfileobj(src, dst)
            else:
                real_run(cmd, **kwargs)

        mock_response = MagicMock()
        mock_response.url = "https://api.ncbi.nlm.nih.gov/datasets/download"
        mock_response.iter_content.return_value = [grch38_zip]

        with patch(
            "q2_quality_control._filter_pangenome.requests.get"
        ) as mock_get, patch(
            "q2_quality_control._filter_pangenome._run_command",
            side_effect=fake_run,
        ):
            mock_get.return_value.__enter__.return_value = mock_response

            (index,) = construct_action(threads=1)

        # the pipeline completed and produced a real Bowtie2 index artifact;
        # reaching here means every (non-network) step ran successfully
        self.assertEqual(str(index.type), "Bowtie2Index")
        index_files = os.listdir(str(index._archiver.data_dir))
        self.assertTrue(
            any(f.endswith((".bt2", ".bt2l")) for f in index_files),
            f"No bowtie2 index files found in {index_files}",
        )

    def test_filter_reads_pangenome(self):
        # we construct our own context so that we can control each "action"
        # being retrieved from it
        ctx = MagicMock()
        ctx.get_action.return_value = MagicMock()
        mock_filtered_reads_result = MagicMock()
        ctx.get_action("quality_control", "filter_reads").return_value = (
            mock_filtered_reads_result,
        )
        mock_index = MagicMock()
        ctx.get_action(
            "quality_control", "construct_human_pangenome_index"
        ).return_value = (mock_index,)

        reads = MagicMock()

        filtered_reads, generated_index = filter_reads_human_pangenome(
            ctx=ctx, reads=reads, index=None, threads=4
        )

        # Assertions
        ctx.get_action.assert_any_call(
            "quality_control", "construct_human_pangenome_index"
        )
        ctx.get_action.assert_any_call("quality_control", "filter_reads")

        ctx.get_action("quality_control", "filter_reads").assert_has_calls(
            [
                call(
                    demultiplexed_sequences=reads,
                    database=generated_index,
                    exclude_seqs=True,
                    n_threads=4,
                    mode="local",
                    ref_gap_open_penalty=5,
                    ref_gap_ext_penalty=3,
                )
            ],
            any_order=True,
        )
        self.assertIsNotNone(generated_index)

    def test_filter_reads_pangenome_with_index(self):
        # distinct mock per action so we can verify that the index
        # construction action is *not* invoked when an index is provided
        mock_construct = MagicMock(name="construct_index")
        mock_filter = MagicMock(name="filter_reads")
        mock_filtered_reads_result = MagicMock()
        mock_filter.return_value = (mock_filtered_reads_result,)

        actions = {
            "construct_human_pangenome_index": mock_construct,
            "filter_reads": mock_filter,
        }
        ctx = MagicMock()
        ctx.get_action.side_effect = lambda plugin, action: actions[action]

        reads = MagicMock()
        provided_index = MagicMock()

        filtered_reads, returned_index = filter_reads_human_pangenome(
            ctx=ctx, reads=reads, index=provided_index, threads=2
        )

        # the index is supplied, so no new index is constructed
        mock_construct.assert_not_called()

        # the provided index is used directly and returned unchanged
        mock_filter.assert_called_once_with(
            demultiplexed_sequences=reads,
            database=provided_index,
            exclude_seqs=True,
            n_threads=2,
            mode="local",
            ref_gap_open_penalty=5,
            ref_gap_ext_penalty=3,
        )
        self.assertIs(returned_index, provided_index)
        self.assertIs(filtered_reads, mock_filtered_reads_result)

    def test_filter_reads_human_pangenome_end_to_end(self):
        # real e2e (no network): filter real reads against a real bowtie2
        # index through the framework. The action hardcodes exclude_seqs=True,
        # so reads mapping to the reference are removed and the unmapped read
        # is kept.

        seq_ids_that_map = [
            "SARS2:6:73:941:1973#", "SARS2:6:73:231:3321#",
            "SARS2:6:73:233:3421#", "SARS2:6:73:552:2457#",
            "SARS2:6:73:567:7631#",
        ]
        seq_id_that_does_not_map = "SARS2:6:73:356:9806#"

        reads = Artifact.load(self.get_data_path("single-end.qza"))
        index = Artifact.load(self.get_data_path("sars2-indexed.qza"))

        filtered_reads, reference_index = filter_action(
            reads=reads, index=index, threads=1
        )

        # the index is returned as the reference_index output
        self.assertEqual(str(reference_index.type), "Bowtie2Index")

        obs = filtered_reads.view(SingleLanePerSampleSingleEndFastqDirFmt)
        observed_ids = []
        for _, obs_fp in obs.sequences.iter_views(FastqGzFormat):
            with gzip.open(str(obs_fp), "rt") as obs_fh:
                lines = obs_fh.readlines()
                self.assertNotEqual(len(lines), 0)
                for records in itertools.zip_longest(*[iter(lines)] * 4):
                    obs_id = records[0].strip("@/012\n")
                    observed_ids.append(obs_id)
                    # mapping reads must have been excluded
                    self.assertNotIn(obs_id, seq_ids_that_map)

        # the one non-mapping read survived the filter
        self.assertIn(seq_id_that_does_not_map, observed_ids)


if __name__ == "__main__":
    unittest.main()
