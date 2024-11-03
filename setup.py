from setuptools import setup, find_packages
from io import open
import os

version_module = {}
dir_name = os.path.dirname(__file__)
with open(os.path.join(dir_name, "src", "babachi", "version.py")) as fp:
    exec(fp.read(), version_module)
    __version__ = version_module['__version__']

with open(os.path.join(dir_name, "README.md"), encoding='utf8') as fh:
    long_description = fh.read()

setup(
    name='babachi',
    version=__version__,
    packages=find_packages('src'),
    include_package_data=True,
    package_data={'babachi': ['tests/*.bed']},
    long_description=long_description,
    entry_points={
        'console_scripts': [
            'babachi = babachi.bad_estimation:segmentation_start',
        ],
    },
    author="Sergey Abramov, Alexandr Boytsov",
    author_email='sabramov@altius.org',
    package_dir={'': 'src'},
    install_requires=[
        'docopt>=0.6.2',
        'numpy>=1.19.5,<2.0',
        'schema>=0.7.2',
        'contextlib2>=0.5.5',
        'pandas>=1.0.4',
        'matplotlib>=3.2.1',
        'seaborn>=0.10.1',
        'numba>=0.53.1',
        'pysam>=0.13.4',
        'scipy>=1.5.1',
        'validators>=0.18.2'
    ],
    extras_require={
        'dev': ['wheel', 'twine', 'setuptools_scm'],
        'vcf': ['pysam']
    },
    python_requires='>=3.6',
    url="https://github.com/autosome-ru/BABACHI",
)
