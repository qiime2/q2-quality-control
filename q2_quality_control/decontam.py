# ----------------------------------------------------------------------------
# Copyright (c) 2017-2022, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import biom
import pandas as pd
import os
import tempfile
import subprocess
import qiime2.util
from ._utilities import _run_command


_WHOLE_NUM = (lambda x: x >= 0, 'non-negative')
_PER_NUM = (lambda x: 1 >= x >= 0, 'between 0 and 1')
_DECON_METHOD_STR = (lambda x: x in {'frequency', 'prevalence', 'combined'},
                     'frequency, prevalence, combined')
_BOOLEAN = (lambda x: type(x) is bool, 'True or False')
# Better to choose to skip, than to implicitly ignore things that KeyError
_SKIP = (lambda x: True, '')
_valid_inputs = {
    'table': _SKIP,
    'metadata': _SKIP,
    'threshold': _PER_NUM,
    'method': _DECON_METHOD_STR,
    'freq_concentration_column': _SKIP,
    'prev_control_column': _SKIP,
    'prev_control_indicator': _SKIP,
    'decontam_scores': _SKIP,
}


def _check_inputs(**kwargs):
    for param, arg in kwargs.items():
        check_is_valid, explanation = _valid_inputs[param]
        if not check_is_valid(arg):
            raise ValueError('Argument to %r was %r, should be %s.'
                             % (param, arg, explanation))


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
                raise ValueError('At least 5 Control Samples needed '
                                 + str(indic) + ' found')
            else:
                print("All appropriate inputs are found")


def _check_column_inputs(table, metadata, method, freq_concentration_column,
                         prev_control_column, prev_control_indicator):
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
    if method == 'combined':
        df = df.fillna(0)
    df = df.dropna(axis='columns')
    return df


def decontam_identify(table: pd.DataFrame,
                      metadata: qiime2.Metadata,
                      method: str = 'prevalence',
                      freq_concentration_column: str = None,
                      prev_control_column: str = None,
                      prev_control_indicator: str = None
                      ) -> (pd.DataFrame):
    _check_inputs(**locals())
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


def decontam_remove(decontam_scores: qiime2.Metadata,
                    table: pd.DataFrame,
                    threshold: float = 0.1
                    ) -> (biom.Table):
    _check_inputs(**locals())
    with tempfile.TemporaryDirectory() as temp_dir_name:
        df = decontam_scores.to_dataframe()
        df.loc[(df['p'].astype(float) <= threshold),
               'contaminant_seq'] = 'True'
        df.loc[(df['p'].astype(float) > threshold),
               'contaminant_seq'] = 'False'
        df = df[df.contaminant_seq == 'True']
        table = table.drop(df.index, axis=1)
        output = os.path.join(temp_dir_name, 'temp.tsv.biom')
        temp_transposed_table = table.transpose()
        temp_transposed_table.to_csv(output, sep="\t")
        with open(output) as fh:
            no_contam_table = biom.Table.from_tsv(fh, None, None, None)
        return no_contam_table
