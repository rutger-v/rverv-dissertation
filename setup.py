from setuptools import setup, find_packages

VERSION = "0.0.1"
DESCRIPTION = "River segment bearing analysis package"
LONG_DESCRIPTION = "Python package containing scripts for analysis of segment bearings of river networks"

setup(
    name = "rivernetworkpy",
    version = VERSION,
    author = "Rutger Vervoordeldonk",
    author_email = "rvervoordeldonk@outlook.com",
    description = DESCRIPTION,
    long_description = LONG_DESCRIPTION,
    packages = find_packages(),
    install_requires = [],
    keywords = ["River network", "Segment bearing", "Geomorphology"],
    classifiers = [
        "Development Status :: 1 - Planning",
        "Intended Audience :: Science/Research"
    ]
)