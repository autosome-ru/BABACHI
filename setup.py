from setuptools import setup, find_packages
from src.babachi import __version__
from os.path import join, dirname

setup(
    name='babachi',
    version=__version__,
    packages=find_packages(),
    package_data={'babachi': ['tests/*.tsv']},
    long_description=open(join(dirname(__file__), 'README.md')).read(),
    entry_points={
        'console_scripts': [
            'babachi = babachi.BADEstimation:segmentation_start',
        ],
    },
    author="Sergey Abramov, Alexandr Boytsov",
    author_email='aswq22013@gmail.com',
    package_dir={'': 'src'},
    install_requires=[
        'docopt>=0.6.2',
        'numpy>=1.18.0',
        'schema>=0.7.2',
        'contextlib2>=0.5.5',
        'pandas>=1.0.4',
        'matplotlib>=3.2.1',
        'seaborn>=0.10.1'
    ],
    python_requires='>=3.6',
    url="https://github.com/autosome-ru/BABACHI",
)
