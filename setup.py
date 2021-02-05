# ----------------------------------------------------------------------------
# Copyright (c) 2017-2021, QIIME 2 development team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

from setuptools import find_packages, setup

import versioneer


setup(
    name='q2-quality-control',
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
    license='BSD-3-Clause',
    packages=find_packages(),
    author="Nicholas Bokulich",
    author_email="nbokulich@gmail.com",
    description="Quality control methods for feature and sequence data.",
    url="https://github.com/qiime2/q2-quality-control",
    entry_points={
        'qiime2.plugins':
        ['q2-quality-control=q2_quality_control.plugin_setup:plugin']
    },
    package_data={
        'q2_quality_control': ['citations.bib', 'assets/*'],
        'q2_quality_control.tests': ['data/*'],
    },
    zip_safe=False,
)
