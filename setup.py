from setuptools import setup, find_packages
from package import __version__
from os.path import join, dirname

setup(
    name='BAD-segmentation',
    version=__version__,
    packages=find_packages(),
    long_description=open(join(dirname(__file__), 'README.txt')).read(),
    entry_points={
        'console_scripts': [
            'segmentation = package.BADEstimation:segmentation_start',
        ],
    },
    install_requires=[
        'docopt>=0.6.2',
        'numpy>=1.18.0',
        'schema>=0.7.2',
        'contextlib2>=0.5.5'
    ],
    python_requires='>=3.6',
    url="https://github.com/wish1/BAD_segmentation",
)
