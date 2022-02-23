import os

from setuptools import setup

base_dir = os.path.dirname(__file__)
install_requires = [line.rstrip() for line in open(os.path.join(base_dir, 'requirements.txt'))]

setup(
    name='lattice-tools',
    version='0.1.0',
    packages=['lattice'],
    url='https://github.com/Lattice-Data/lattice-tools/',
    license='MIT',
    author='Lattice Data Coordination',
    author_email='',
    description='External scripts used to interact with the Lattice Database',
    install_requires=install_requires
)
