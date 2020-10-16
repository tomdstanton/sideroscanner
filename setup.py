#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from setuptools import setup, find_packages

with open('requirements.txt') as f:
    requirements = f.read().splitlines()

setup(
    name='SideroScanner',
    version='0.0.2',
    packages=find_packages(),
    include_package_data = True,
    package_data={'sideroscanner': ['sideroscanner/data/irompdb/iromps.csv',
                                    'sideroscanner/data/irompdb/iromps.hmm',
                                    'sideroscanner/data/furdb/fur.meme']},
    scripts=['sideroscanner/sideroscanner',
             'sideroscanner/sideroscanner-buildhmms',
             'sideroscanner/sideroscanner-builddbs'],
    url='https://github.com/tomdstanton/sideroscanner',
    license='gpl-3.0',
    author='Tom Stanton',
    author_email='T.D.Stanton@sms.ed.ac.uk',
    keywords='microbial genomics amr virulence',
    description='A tool for annotating IROMPs in bacteria',
    long_description='Accuracy annotates Iron Regulated Outer Membrane '
                     'Proteins in bacteria and aims to provide as much '
                     'biologically relevent information as possible for '
                     'downstream wet/dry-lab validation.',
    install_requires=requirements,
    zip_safe=False,
    test_suite='nose.collector',
    tests_require=['nose'],
)
