"""
Setup file for package publishing
Created: 31/03/2026
Python 3.9.7
Carolin Sauer
"""
import setuptools

from copybara import helper

with open("README.md", encoding="utf-8") as fh:
    long_description = fh.read()

setuptools.setup(
    name="copybara-cf",
    version=f"{helper.__version__}",
    author="Carolin Sauer",
    author_email="csauer@ebi.ac.uk",
    url="https://github.com/cortes-ciriano-lab/copybara-cf",
    description="COPYBARA - copy number analysis for cfDNA sequencing data",
    long_description=long_description,
    long_description_content_type="text/markdown",
    packages=setuptools.find_packages(),
    entry_points={
        'console_scripts': [
            'copybara=copybara.copybara:main',
        ]
    },
    include_package_data=True,
    classifiers=(
        "Programming Language :: Python :: 3.9",
        "Operating System :: Unix",
        "Development Status :: 5 - Production/Stable",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics"
    ),
)
