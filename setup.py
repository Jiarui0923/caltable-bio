from setuptools import setup, find_packages

LONG_DESCRIPTION = '''
The biology computational extension for `CalTable` package.
It provides data engine and computaton blocks for proteins data.

If there is any issue, please put up with an issue or contact Jiarui Li (jli78@tulane.edu)
'''
VERSION = '0.0.8'
NAME = 'CalTable-Bio'

install_requires = [
    "caltable @ git+https://github.com/Jiarui0923/CalTable@2.0.8",
    "py3Dmol",
    "plotly",
    "matplotlib"
]


setup(
    name=NAME,
    description="Calculate Table (CalTable) Biology Extention.",
    long_description=LONG_DESCRIPTION,
    long_description_content_type="text/markdown",
    version=VERSION,
    packages=find_packages(),
    install_requires=install_requires,
    url="https://git.tulane.edu/apl/caltable",
    author='Jiarui Li, Marco K. Carbullido, Jai Bansal, Samuel J. Landry, Ramgopal R. Mettu',
    author_email=('jli78@tulane.edu'),
    maintainer=("Jiarui Li"),
    maintainer_email="jli78@tulane.edu",
    license='GPLv3',
    classifiers=[
        "Intended Audience :: Developers",
        "Intended Audience :: Science/Research",
        "Operating System :: MacOS :: MacOS X",
        "Operating System :: Microsoft :: Windows",
        "Operating System :: POSIX",
        "Operating System :: OS Independent",
        'License :: OSI Approved :: GPLv3 License',
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.12",
        "Topic :: Software Development :: Build Tools",
    ],
    zip_safe=True,
    platforms=["any"],
    entry_points={
        'caltable.extensions': [
            'caltable_bio = caltable_bio'
        ],
    },
)
