# ----------------------------------------------------------------------------
# Copyright (c) 2017-2025, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------
import unittest
import pandas as pd
import qiime2
import qiime2.plugin.util
import biom
from qiime2.plugin.testing import TestPluginBase
from q2_quality_control.decontam import (decontam_identify)
from skbio.sequence import DNA

from q2_quality_control._threshold_graph._visualizer import (
    decontam_score_viz)

import os
import tempfile


class TestIdentify(TestPluginBase):
    package = 'q2_quality_control.tests'

    def setUp(self):
        super().setUp()
        table = qiime2.Artifact.load(
            self.get_data_path('expected/decon_default_ASV_table.qza'))
        self.asv_table = table.view(qiime2.Metadata).to_dataframe()
        self.metadata_input = qiime2.Metadata.load(
            self.get_data_path('expected/test_metadata.tsv'))

    def test_prevalence(self):
        exp_table = pd.read_csv(
            self.get_data_path('expected/prevalence-score-table.tsv'),
            sep='\t', index_col=0)
        temp_transposed_table = exp_table.transpose()
        temp_transposed_table = temp_transposed_table.dropna()
        exp_table = temp_transposed_table.transpose()
        df_output_feature_table = decontam_identify(
            table=self.asv_table,
            metadata=self.metadata_input,
            method='prevalence',
            prev_control_column='Sample_or_Control',
            prev_control_indicator='Control Sample')
        df_output_feature_table = df_output_feature_table.round(decimals=6)
        exp_table = exp_table.round(decimals=6)
        with tempfile.TemporaryDirectory() as temp_dir_name:
            test_biom_fp = os.path.join(temp_dir_name, 'test_output.tsv')
            expected_biom_fp = os.path.join(temp_dir_name,
                                            'expected_output.tsv')
            df_output_feature_table.to_csv(test_biom_fp, sep="\t")
            exp_table.to_csv(expected_biom_fp, sep="\t")
            with open(test_biom_fp) as fh:
                test_table = biom.Table.from_tsv(fh, None, None, None)
            with open(expected_biom_fp) as th:
                expecter_table = biom.Table.from_tsv(th, None, None, None)

            self.assertEqual(test_table, expecter_table)

    def test_frequency(self):
        exp_table = pd.read_csv(
            self.get_data_path('expected/frequency-score-table.tsv'),
            sep='\t', index_col=0)
        temp_transposed_table = exp_table.transpose()
        temp_transposed_table = temp_transposed_table.dropna()
        exp_table = temp_transposed_table.transpose()
        df_output_feature_table = decontam_identify(
            table=self.asv_table,
            metadata=self.metadata_input,
            method='frequency',
            freq_concentration_column='quant_reading')
        df_output_feature_table = df_output_feature_table.round(decimals=6)
        exp_table = exp_table.round(decimals=6)
        with tempfile.TemporaryDirectory() as temp_dir_name:
            test_biom_fp = os.path.join(temp_dir_name, 'test_output.tsv')
            expected_biom_fp = os.path.join(temp_dir_name,
                                            'expected_output.tsv')
            df_output_feature_table.to_csv(test_biom_fp, sep="\t")
            exp_table.to_csv(expected_biom_fp, sep="\t")
            with open(test_biom_fp) as fh:
                test_table = biom.Table.from_tsv(fh, None, None, None)
            with open(expected_biom_fp) as th:
                expecter_table = biom.Table.from_tsv(th, None, None, None)
            self.assertEqual(test_table, expecter_table)

    def test_combined(self):
        exp_table = pd.read_csv(
            self.get_data_path('expected/combined-score-table.tsv'),
            sep='\t', index_col=0)
        df_output_feature_table = decontam_identify(
            table=self.asv_table,
            metadata=self.metadata_input,
            method='combined',
            prev_control_column='Sample_or_Control',
            prev_control_indicator='Control Sample',
            freq_concentration_column='quant_reading')
        df_output_feature_table = df_output_feature_table.fillna(0)
        df_output_feature_table = df_output_feature_table.round(decimals=6)
        exp_table = exp_table.fillna(0)
        exp_table = exp_table.round(decimals=6)
        with tempfile.TemporaryDirectory() as temp_dir_name:
            test_biom_fp = os.path.join(temp_dir_name, 'test_output.tsv')
            expected_biom_fp = os.path.join(temp_dir_name,
                                            'expected_output.tsv')
            df_output_feature_table.to_csv(test_biom_fp, sep="\t")
            exp_table.to_csv(expected_biom_fp, sep="\t")
            with open(test_biom_fp) as fh:
                test_table = biom.Table.from_tsv(fh, None, None, None)
            with open(expected_biom_fp) as th:
                expecter_table = biom.Table.from_tsv(th, None, None, None)
            self.assertEqual(test_table, expecter_table)

    def test_freq_prev_control_col_error_raised(self):
        with self.assertRaises(ValueError) as context:
            decontam_identify(
                table=self.asv_table,
                metadata=self.metadata_input,
                method='frequency',
                freq_concentration_column='quant_reading',
                prev_control_column='Sample_or_Control')
        self.assertIn("--p-prev-control-column given, but cannot be used",
                      str(context.exception))

    def test_freq_prev_control_col_indic_error_raised(self):
        with self.assertRaises(ValueError) as context:
            decontam_identify(
                table=self.asv_table,
                metadata=self.metadata_input,
                method='frequency',
                freq_concentration_column='quant_reading',
                prev_control_indicator='Control Sample')
        self.assertIn("--p-prev-control-indicator given, but cannot be used",
                      str(context.exception))

    def test_freq_params_error_raised(self):
        with self.assertRaises(ValueError) as context:
            decontam_identify(
                table=self.asv_table,
                metadata=self.metadata_input,
                method='frequency')
        self.assertIn("--p-freq-concentration-column is being utilized",
                      str(context.exception))

    def test_prev_no_indic_param_error_raised(self):
        with self.assertRaises(ValueError) as context:
            decontam_identify(
                table=self.asv_table,
                metadata=self.metadata_input,
                method='prevalence',
                prev_control_column='Sample_or_Control')
        self.assertIn("--p-prev-control-column and --p-prev-control-indicator",
                      str(context.exception))

    def test_prev_no_col_param_error_raised(self):
        with self.assertRaises(ValueError) as context:
            decontam_identify(
                table=self.asv_table,
                metadata=self.metadata_input,
                method='prevalence',
                prev_control_indicator='Control Sample')
        self.assertIn("--p-prev-control-column and --p-prev-control-indicator",
                      str(context.exception))

    def test_prev_freq_params_error_raised(self):
        with self.assertRaises(ValueError) as context:
            decontam_identify(
                table=self.asv_table,
                metadata=self.metadata_input,
                method='prevalence',
                prev_control_indicator='Control Sample',
                prev_control_column='Sample_or_Control',
                freq_concentration_column='quant_reading')
        self.assertIn("--p-freq-concentration-column given, but cannot be",
                      str(context.exception))


class TestIdentify_mixed_names(TestPluginBase):
    package = 'q2_quality_control.tests'

    def setUp(self):
        super().setUp()
        table = qiime2.Artifact.load(
            self.get_data_path('expected/decon_default_ASV_table.qza'))
        self.asv_table = table.view(qiime2.Metadata).to_dataframe()
        self.metadata_input = qiime2.Metadata.load(
            self.get_data_path(
                'expected/test_metadata_different_ordered_names.tsv'))

    def test_prevalence(self):
        exp_table = pd.read_csv(
            self.get_data_path('expected/prevalence-score-table.tsv'),
            sep='\t', index_col=0)
        temp_transposed_table = exp_table.transpose()
        temp_transposed_table = temp_transposed_table.dropna()
        exp_table = temp_transposed_table.transpose()
        df_output_feature_table = decontam_identify(
            table=self.asv_table,
            metadata=self.metadata_input,
            method='prevalence',
            prev_control_column='Sample_or_Control',
            prev_control_indicator='Control Sample')
        df_output_feature_table = df_output_feature_table.round(decimals=6)
        exp_table = exp_table.round(decimals=6)
        with tempfile.TemporaryDirectory() as temp_dir_name:
            test_biom_fp = os.path.join(temp_dir_name, 'test_output.tsv')
            expected_biom_fp = os.path.join(temp_dir_name,
                                            'expected_output.tsv')
            df_output_feature_table.to_csv(test_biom_fp, sep="\t")
            exp_table.to_csv(expected_biom_fp, sep="\t")
            with open(test_biom_fp) as fh:
                test_table = biom.Table.from_tsv(fh, None, None, None)
            with open(expected_biom_fp) as th:
                expecter_table = biom.Table.from_tsv(th, None, None, None)

            self.assertEqual(test_table, expecter_table)

    def test_frequency(self):
        exp_table = pd.read_csv(
            self.get_data_path('expected/frequency-score-table.tsv'),
            sep='\t', index_col=0)
        temp_transposed_table = exp_table.transpose()
        temp_transposed_table = temp_transposed_table.dropna()
        exp_table = temp_transposed_table.transpose()
        df_output_feature_table = decontam_identify(
            table=self.asv_table,
            metadata=self.metadata_input,
            method='frequency',
            freq_concentration_column='quant_reading')
        df_output_feature_table = df_output_feature_table.round(decimals=6)
        exp_table = exp_table.round(decimals=6)
        with tempfile.TemporaryDirectory() as temp_dir_name:
            test_biom_fp = os.path.join(temp_dir_name, 'test_output.tsv')
            expected_biom_fp = os.path.join(temp_dir_name,
                                            'expected_output.tsv')
            df_output_feature_table.to_csv(test_biom_fp, sep="\t")
            exp_table.to_csv(expected_biom_fp, sep="\t")
            with open(test_biom_fp) as fh:
                test_table = biom.Table.from_tsv(fh, None, None, None)
            with open(expected_biom_fp) as th:
                expecter_table = biom.Table.from_tsv(th, None, None, None)
            self.assertEqual(test_table, expecter_table)

    def test_combined(self):
        exp_table = pd.read_csv(
            self.get_data_path('expected/combined-score-table.tsv'),
            sep='\t', index_col=0)
        df_output_feature_table = decontam_identify(
            table=self.asv_table,
            metadata=self.metadata_input,
            method='combined',
            prev_control_column='Sample_or_Control',
            prev_control_indicator='Control Sample',
            freq_concentration_column='quant_reading')
        df_output_feature_table = df_output_feature_table.fillna(0)
        df_output_feature_table = df_output_feature_table.round(decimals=6)
        exp_table = exp_table.fillna(0)
        exp_table = exp_table.round(decimals=6)
        with tempfile.TemporaryDirectory() as temp_dir_name:
            test_biom_fp = os.path.join(temp_dir_name, 'test_output.tsv')
            expected_biom_fp = os.path.join(temp_dir_name,
                                            'expected_output.tsv')
            df_output_feature_table.to_csv(test_biom_fp, sep="\t")
            exp_table.to_csv(expected_biom_fp, sep="\t")
            with open(test_biom_fp) as fh:
                test_table = biom.Table.from_tsv(fh, None, None, None)
            with open(expected_biom_fp) as th:
                expecter_table = biom.Table.from_tsv(th, None, None, None)
            self.assertEqual(test_table, expecter_table)


class TestIdentify_more_names(TestPluginBase):
    package = 'q2_quality_control.tests'

    def setUp(self):
        super().setUp()
        table = qiime2.Artifact.load(
            self.get_data_path('expected/decon_default_ASV_table.qza'))
        self.asv_table = table.view(qiime2.Metadata).to_dataframe()
        self.metadata_input = qiime2.Metadata.load(
            self.get_data_path(
                'expected/test_metadata_more_sample_names.tsv'))

    def test_prevalence(self):
        exp_table = pd.read_csv(
            self.get_data_path('expected/prevalence-score-table.tsv'),
            sep='\t', index_col=0)
        temp_transposed_table = exp_table.transpose()
        temp_transposed_table = temp_transposed_table.dropna()
        exp_table = temp_transposed_table.transpose()
        df_output_feature_table = decontam_identify(
            table=self.asv_table,
            metadata=self.metadata_input,
            method='prevalence',
            prev_control_column='Sample_or_Control',
            prev_control_indicator='Control Sample')
        df_output_feature_table = df_output_feature_table.round(decimals=6)
        exp_table = exp_table.round(decimals=6)
        with tempfile.TemporaryDirectory() as temp_dir_name:
            test_biom_fp = os.path.join(temp_dir_name, 'test_output.tsv')
            expected_biom_fp = os.path.join(temp_dir_name,
                                            'expected_output.tsv')
            df_output_feature_table.to_csv(test_biom_fp, sep="\t")
            exp_table.to_csv(expected_biom_fp, sep="\t")
            with open(test_biom_fp) as fh:
                test_table = biom.Table.from_tsv(fh, None, None, None)
            with open(expected_biom_fp) as th:
                expecter_table = biom.Table.from_tsv(th, None, None, None)

            self.assertEqual(test_table, expecter_table)

    def test_frequency(self):
        exp_table = pd.read_csv(
            self.get_data_path('expected/frequency-score-table.tsv'),
            sep='\t', index_col=0)
        temp_transposed_table = exp_table.transpose()
        temp_transposed_table = temp_transposed_table.dropna()
        exp_table = temp_transposed_table.transpose()
        df_output_feature_table = decontam_identify(
            table=self.asv_table,
            metadata=self.metadata_input,
            method='frequency',
            freq_concentration_column='quant_reading')
        df_output_feature_table = df_output_feature_table.round(decimals=6)
        exp_table = exp_table.round(decimals=6)
        with tempfile.TemporaryDirectory() as temp_dir_name:
            test_biom_fp = os.path.join(temp_dir_name, 'test_output.tsv')
            expected_biom_fp = os.path.join(temp_dir_name,
                                            'expected_output.tsv')
            df_output_feature_table.to_csv(test_biom_fp, sep="\t")
            exp_table.to_csv(expected_biom_fp, sep="\t")
            with open(test_biom_fp) as fh:
                test_table = biom.Table.from_tsv(fh, None, None, None)
            with open(expected_biom_fp) as th:
                expecter_table = biom.Table.from_tsv(th, None, None, None)
            self.assertEqual(test_table, expecter_table)

    def test_combined(self):
        exp_table = pd.read_csv(
            self.get_data_path('expected/combined-score-table.tsv'),
            sep='\t', index_col=0)
        df_output_feature_table = decontam_identify(
            table=self.asv_table,
            metadata=self.metadata_input,
            method='combined',
            prev_control_column='Sample_or_Control',
            prev_control_indicator='Control Sample',
            freq_concentration_column='quant_reading')
        df_output_feature_table = df_output_feature_table.fillna(0)
        df_output_feature_table = df_output_feature_table.round(decimals=6)
        exp_table = exp_table.fillna(0)
        exp_table = exp_table.round(decimals=6)
        with tempfile.TemporaryDirectory() as temp_dir_name:
            test_biom_fp = os.path.join(temp_dir_name, 'test_output.tsv')
            expected_biom_fp = os.path.join(temp_dir_name,
                                            'expected_output.tsv')
            df_output_feature_table.to_csv(test_biom_fp, sep="\t")
            exp_table.to_csv(expected_biom_fp, sep="\t")
            with open(test_biom_fp) as fh:
                test_table = biom.Table.from_tsv(fh, None, None, None)
            with open(expected_biom_fp) as th:
                expecter_table = biom.Table.from_tsv(th, None, None, None)
            self.assertEqual(test_table, expecter_table)


class TestVizualization(TestPluginBase):
    package = 'q2_quality_control.tests'

    def setUp(self):
        super().setUp()
        self.input_table = {'test_dict': pd.DataFrame(
            [[1, 2, 3, 4, 5], [9, 10, 11, 12, 13]],
            columns=['abc', 'def', 'jkl', 'mno', 'pqr'],
            index=['sample-1', 'sample-2'])}
        temp_list = []
        seqs = ['ACGT', 'TTTT', 'AAAA', 'CCCC', 'GGG']
        indexes = ['abc', 'def', 'jkl', 'mno', 'pqr']
        for i in range(len(seqs)):
            temp_list.append(DNA(seqs[i], metadata={'id': indexes[i]}))
        self.input_seqs = temp_list
        self.input_scores = {'test_dict': pd.DataFrame(
            [[13.0, 0.969179],
             [16.0, 0.566067],
             [25.0, 0.019475],
             [10.0, 0.383949],
             [13.0, 0.969179]],
            index=['abc', 'def', 'jkl', 'mno', 'pqr'],
            columns=['prev', 'p'])}
        self.output_dir_obj = tempfile.TemporaryDirectory(
                prefix='q2-quality-control-decontam-test-temp-')
        self.output_dir = self.output_dir_obj.name

    def tearDown(self):
        self.output_dir_obj.cleanup()

    def assertScore_Viz_Basics(self, viz_dir):
        index_fp = os.path.join(viz_dir, 'index.html')
        self.assertTrue(os.path.exists(index_fp))
        with open(index_fp, 'r') as fh:
            index_contents = fh.read()
        self.assertIn('<td>1</td>\n                    '
                      '<td>4</td>\n                    '
                      '<td>20.00</td>\n', index_contents)
        self.assertTrue(os.path.exists(
            os.path.join(viz_dir, 'test_dict-identify-table-histogram.png')))

    def test_defaults(self):
        decontam_score_viz(output_dir=self.output_dir,
                           table=self.input_table,
                           decontam_scores=self.input_scores,
                           threshold=0.1,
                           rep_seqs=self.input_seqs, weighted=False)
        self.assertScore_Viz_Basics(self.output_dir)

    def test_defaults_no_reps(self):
        decontam_score_viz(output_dir=self.output_dir,
                           table=self.input_table,
                           decontam_scores=self.input_scores,
                           threshold=0.1,
                           rep_seqs=None, weighted=False)
        self.assertScore_Viz_Basics(self.output_dir)


if __name__ == '__main__':
    unittest.main()
