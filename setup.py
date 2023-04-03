from setuptools import setup, find_packages
import src.gdt.core as gdt_c

if __name__ == '__main__':
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
                'Sphinx',
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

