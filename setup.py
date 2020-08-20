from setuptools import setup, find_packages
import os

version_module = {}
dir_name = os.path.dirname(__file__)
with open(os.path.join(dir_name, "src/babachi/version.py")) as fp:
    exec(fp.read(), version_module)
    __version__ = version_module['__version__']

setup(
    name='babachi',
    version=__version__,
    packages=find_packages(),
    package_data={'babachi': ['tests/*.tsv']},
    long_description=open(os.path.join(dir_name, 'README.md')).read(),
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
