from setuptools import setup, find_packages

setup(
    name='vanheusden',
    version='1.0',
    description='Genome Size Estimator',
    author='Jonathan Lam, Richard Xu, Edward Tang',
    author_email='email@ucsd.edu',
    packages=find_packages(),
    entry_points={
        "console_scripts": [
            "vanheusden=vanheusden.vanheusden:main"
        ],
    }
)