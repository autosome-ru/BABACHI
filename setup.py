from setuptools import setup, find_packages
from os.path import join, dirname

setup(
    name='BAD_segmentation',
    version='1.0',
    packages=find_packages(),
    long_description=open(join(dirname(__file__), 'README.txt')).read(),
    entry_points={
        'console_scripts': [
            'segmentation = scripts.PloidyEstimation:test',
        ],
    },
    install_requires=[
        'docopt>=0.6.2',
        'numpy>=1.18.0',
    ]
)
