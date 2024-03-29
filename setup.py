from setuptools import setup, find_packages

setup(
    name='resistnet',
    version='1.2.0',
    # Specify where to find the packages
    package_dir={'': 'src'},
    packages=find_packages(where='src'),
    scripts=[
        'scripts/runResistnet.py',
        'scripts/simResistnet.py',
        'scripts/ensembleResistnet.py'
    ],
    install_requires=[
        "numpy",
        "pandas",
        "networkx",
        "pytest",
        "rpy2",
        "seaborn",
        "matplotlib",
        "pyogrio",
        "momepy",
        "geopy",
        "deap",
        "scipy"
    ],
    include_package_data=True,
    author='Tyler Chafin',
    author_email='tyler.chafin@bioss.ac.uk',
    description='A package for modelling environmental resistance networks',
    keywords='population genetics, environmental resistance, river networks',
    url='https://github.com/tkchafin/resistnet',
)