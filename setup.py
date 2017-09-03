#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""The setup script."""

from setuptools import setup, find_packages

with open('README.md') as readme_file:
    readme = readme_file.read()

with open('INSTALL.md') as install_file:
    install = install_file.read()

requirements = [
    'numpy>=1.13.0'
    # TODO: put package requirements here
]

setup_requirements = [
    'pytest-runner',
    # TODO(andfranklin): put setup requirements (distutils extensions, etc.) here
]

test_requirements = [
    'pytest',
    # TODO: put package test requirements here
]

setup(
    name='sem_python',
    version='0.1.0',
    description="A 1D implementation of the 7-equation model in python.",
    long_description=readme + "\n\n" + install,
    author="Joshua Hansel",
    author_email='joshua.hansel@inl.gov',
    url='https://github.com/jhansel/sem_python',
    packages=find_packages(include=['sem_python']),
    entry_points={
        'console_scripts': [
            'sem=sem:main'
        ]
    },
    include_package_data=True,
    install_requires=requirements,
    license="MIT license",
    zip_safe=False,
    keywords='sem_python',
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: MIT License',
        'Natural Language :: English',
        "Programming Language :: Python :: 2",
        'Programming Language :: Python :: 2.6',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.3',
        'Programming Language :: Python :: 3.4',
        'Programming Language :: Python :: 3.5',
    ],
    test_suite='testing',
    # tests_require=test_requirements,
    # setup_requires=setup_requirements,
)
