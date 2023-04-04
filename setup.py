# CONTAINS TECHNICAL DATA/COMPUTER SOFTWARE DELIVERED TO THE U.S. GOVERNMENT WITH UNLIMITED RIGHTS
#
# Contract No.: CA 80MSFC17M0022
# Contractor Name: Universities Space Research Association
# Contractor Address: 7178 Columbia Gateway Drive, Columbia, MD 21046
#
# Copyright 2017-2022 by Universities Space Research Association (USRA). All rights reserved.
#
# Developed by: William Cleveland and Adam Goldstein
#               Universities Space Research Association
#               Science and Technology Institute
#               https://sti.usra.edu
#
# Developed by: Daniel Kocevski
#               National Aeronautics and Space Administration (NASA)
#               Marshall Space Flight Center
#               Astrophysics Branch (ST-12)
#
# Licensed under the Apache License, Version 2.0 (the "License"); you may not use this file except
# in compliance with the License. You may obtain a copy of the License at
#
#    http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software distributed under the License
# is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or
# implied. See the License for the specific language governing permissions and limitations under the
# License.
#
import sys
from setuptools import setup

sys.path.append('src')

if __name__ == '__main__':
<<<<<<< HEAD
    setup(
        name='gdt-core',
        version=gdt_c.__version__,
        scripts=[
            "scripts/gdt-data"
        ],
        packages=find_packages(where='src'),
        package_dir={'': 'src'},
        url='github.com/USRA-STI/gdt',
        license='Apache 2.0',
        author='Cleveland, Goldstein, Kocevski',
        description='The Gamma-ray Data Tools (Core functions)',
        python_requires='>=3.8',
        install_requires=[
            'pyproj>=1.9.6',
            'numpy>=1.17.3',
            'scipy>=1.1.0',
            'matplotlib',
            'astropy>=3.1',
            'healpy>=1.12.4',
            'cartopy',
            'rich'
        ],
        extras_require={
            'docs': [
                'Sphinx==4.5.0', # for some reason newer versions break API docs
                'astropy_sphinx_theme',
                'nbsphinx',
                'ipython',
                'sphinx_automodapi',
                'notebook'
            ],
            'test': [
                'pytest'
            ]
        }
    )
    # Create the GDT data directory, if it doesn't exist
    gdt_c.data_path.mkdir(parents=True, exist_ok=True)
=======
    import gdt.core as core
>>>>>>> ec114d3 (moved the package details to pyproject.toml)

    # For backward compatibility
    setup()

    # create library data directory
    core.data_path.mkdir(parents=True, exist_ok=True)
