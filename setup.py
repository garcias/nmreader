from setuptools import setup, find_packages

setup(
    name = 'nmreader',
    version = 0.1,
    description = 'Parse NMReady JCAMP-DX files for NMR data',
    url = 'https://github.com/garcias/nmreadyreader',
    author = 'Simon Garcia',
    author_email = 'garcias@kenyon.edu',
    license='MIT',
    packages=find_packages(),
    install_requires=['numpy', 'scipy', 'altair', 'pandas', 'dataclasses', ],
    python_requires='>=3.6',
    keywords='nmr parsing jcamp dx',
)

