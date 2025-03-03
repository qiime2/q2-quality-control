# ----------------------------------------------------------------------------
# Copyright (c) 2017-2025, QIIME 2 development team.
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


def _check_prev_inputs(table, metadata, prev_control_column,
                       prev_control_indicator, method):
    if prev_control_column is None or prev_control_indicator is None:
        raise ValueError('For the ' + str(method) + ' method '
                         'please check that input parameters'
                         '--p-prev-control-column and '
                         '--p-prev-control-indicator are being utilized')
    if prev_control_column not in metadata.columns:
        raise ValueError('Prevalence column not found, please '
                         'select from:\n'
                         + str(', '.join(metadata.columns)))
    if prev_control_indicator not in list(metadata[prev_control_column]):
        raise ValueError('No control values found, please select '
                         'from:\n' +
                         str(', '.join(metadata[prev_control_column]
                                       .unique())))
    prev_controls = metadata.loc[
        metadata[prev_control_column] == prev_control_indicator]
    indic = prev_controls.index.intersection(table.index).size
    if indic < 5:
        print("We recommend 5 control samples - " +
              str(indic) + " found")


def _check_freq_inputs(metadata, freq_concentration_column, method):
    if freq_concentration_column is None:
        raise ValueError('For the ' + str(method) + ' method please check '
                         'that input parameter'
                         ' --p-freq-concentration-column is '
                         'being utilized')
    if freq_concentration_column not in metadata.columns:
        raise ValueError('Frequency column not found, please '
                         'select from:\n'
                         + str(', '.join(metadata.columns)))


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
        if freq_concentration_column is not None:
            raise ValueError('--p-freq-concentration-column given, but'
                             ' cannot be used with the prevalence method')
        _check_prev_inputs(table, metadata, prev_control_column,
                           prev_control_indicator, method)
    elif method == 'frequency':
        if prev_control_column is not None:
            raise ValueError('--p-prev-control-column given, but'
                             ' cannot be used with the frequency method')
        if prev_control_indicator is not None:
            raise ValueError('--p-prev-control-indicator given, but'
                             ' cannot be used with the frequency method')
        _check_freq_inputs(metadata, freq_concentration_column, method)
    else:
        _check_prev_inputs(table, metadata, prev_control_column,
                           prev_control_indicator, method)
        _check_freq_inputs(metadata, freq_concentration_column, method)


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
