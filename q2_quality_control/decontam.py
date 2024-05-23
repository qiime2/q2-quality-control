# ----------------------------------------------------------------------------
# Copyright (c) 2017-2022, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import pandas as pd
import os
import tempfile
import subprocess
import qiime2.util
from ._utilities import _run_command


def _check_column_inputs_helper(table, metadata,
                                prev_control_column,
                                prev_control_indicator):
    if prev_control_column not in metadata.columns:
        raise ValueError('Prevalence column not found, please '
                         'select from:\n'
                         + str(', '.join(metadata.columns)))
    else:
        if prev_control_indicator not in list(
                metadata[prev_control_column]):
            raise ValueError('No control values found, please select '
                             'from:\n' +
                             str(', '.join(
                                 metadata[prev_control_column]
                                 .unique())))
        else:
            prev_controls = metadata.loc[
                metadata[prev_control_column] == prev_control_indicator]
            prev_sample_names = prev_controls.index.values
            indic = 0
            for name in prev_sample_names:
                if name in table.index.values:
                    indic = indic + 1
            if indic < 5:
                '''raise ValueError('At least 5 Control Samples needed '
                                 + str(indic) + ' found')'''
                print("We recommend 5 control samples - " +
                      str(indic) + " found")
            else:
                print("All appropriate inputs are found")


def _check_column_inputs(table, metadata, method, freq_concentration_column,
                         prev_control_column, prev_control_indicator):
    meta_sample_names = metadata.index.values
    no_info = []
    for name in table.index.values:
        if name not in meta_sample_names:
            no_info.append(name)
    if len(no_info) > 0:
        raise ValueError('The following samples have no '
                         'metadata:\n'
                         + str(', '.join(no_info)))
    if method == 'prevalence':
        _check_column_inputs_helper(table, metadata,
                                    prev_control_column,
                                    prev_control_indicator)

    elif method == 'frequency':
        if freq_concentration_column not in metadata.columns:
            raise ValueError('Frequency column not found, please '
                             'select from:\n'
                             + str(', '.join(metadata.columns)))
        else:
            print("All appropriate inputs are found")
    else:
        if ((freq_concentration_column not in metadata.columns) or
                (prev_control_column not in metadata.columns)):
            raise ValueError('Column id input error, please '
                             'select from:\n'
                             + str(', '.join(metadata.columns)))
        else:
            _check_column_inputs_helper(table, metadata,
                                        prev_control_column,
                                        prev_control_indicator)


def _decontam_identify_helper(track_fp, method):
    df = pd.read_csv(track_fp, sep='\t', index_col=0)
    df.index.name = '#OTU ID'
    df = df.drop(df.columns[(len(df.columns)-1)], axis=1)
    if method == 'frequency':
        df = df.drop(df.columns[[1, 2, 3]], axis=1)
    elif method == 'prevalence':
        df = df.drop(df.columns[[0, 2, 3]], axis=1)
    else:
        print("We need all of these columns")
    return df


def decontam_identify(table: pd.DataFrame,
                      metadata: qiime2.Metadata,
                      method: str = 'prevalence',
                      freq_concentration_column: str = None,
                      prev_control_column: str = None,
                      prev_control_indicator: str = None
                      ) -> (pd.DataFrame):
    metadata = metadata.to_dataframe()
    _check_column_inputs(table, metadata, method, freq_concentration_column,
                         prev_control_column, prev_control_indicator)
    with tempfile.TemporaryDirectory() as temp_dir_name:
        track_fp = os.path.join(temp_dir_name, 'track.tsv')
        ASV_dest = os.path.join(temp_dir_name, 'temp_ASV_table.csv')
        table.to_csv(os.path.join(ASV_dest))
        meta_dest = os.path.join(temp_dir_name, 'temp_metadata.csv')
        metadata.to_csv(os.path.join(meta_dest))

        cmd = ['run_decontam.R',
               '--asv_table_path', str(ASV_dest),
               '--threshold', str(0.1),
               '--decon_method', method,
               '--output_track', track_fp,
               '--meta_table_path', str(meta_dest),
               '--freq_con_column', str(freq_concentration_column),
               '--prev_control_or_exp_sample_column',
               str(prev_control_column),
               '--prev_control_sample_indicator',
               str(prev_control_indicator)]
        try:
            _run_command(cmd)
        except subprocess.CalledProcessError as e:
            if e.returncode == 2:
                raise ValueError("There was an issue running "
                                 "run_decontam.R please check your inputs")
            else:
                raise Exception("An error was encountered "
                                "while running Decontam in R "
                                "(return code %d), please inspect stdout"
                                " and stderr to learn more." % e.returncode)
        return _decontam_identify_helper(track_fp, method)


def decontam_remove(decontam_scores: pd.DataFrame,
                    table: pd.DataFrame,
                    rep_seqs: pd.Series,
                    threshold: float = 0.1
                    ) -> (pd.DataFrame, pd.Series):
    _check_inputs(**locals())
    decontam_scores['contaminant_seq'] = \
        decontam_scores['p'].astype(float) <= threshold

    decontam_scores = decontam_scores[decontam_scores['contaminant_seq']]
    table.drop(decontam_scores.index, axis=1, inplace=True)
    rep_seqs.drop(decontam_scores.index, inplace=True)

    return table, rep_seqs
