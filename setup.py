from setuptools import setup, find_packages
import os

setup(
    name='smolyak',
    version='0.0.1',
    author='EconForge',
    author_email='slyon@stern.nyu.edu',
    packages=find_packages(),
    url='https://github.com/EconForge/Smolyak.jl',
    description='Smolyak Grid interpolation.',
    long_description=open('README.md').read()
)
