# setup.py

from setuptools import setup, find_packages

setup(
    name='vanheusden',
    version='0.1',
    packages=find_packages(),
    entry_points={
        'console_scripts': [
            'vanheusden=vanheusden.kmer_analysis:main',
        ],
    },
)
