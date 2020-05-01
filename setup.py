from setuptools import setup, find_packages
from os.path import join, dirname

setup(
    name='BAD-segmentation',
    version='0.1',
    packages=find_packages(),
    long_description=open(join(dirname(__file__), 'README.txt')).read(),
    entry_points={
        'console_scripts': [
            'segmentation = scripts.PloidyEstimation:segmentation_start',
        ],
    },
    install_requires=[
        'docopt>=0.6.2',
        'numpy>=1.18.0',
        'schema>=0.7.2',
        'contextlib2>=0.5.5'
    ]
)
