from setuptools import setup, find_packages

setup(
    name = 'nmreader',
    version = 0.1,
    description = 'Parse NMReady JCAMP-DX files for NMR data',
    url = 'https://github.com/garcias/nmreader',
    author = 'Simon Garcia',
    author_email = 'garcias@kenyon.edu',
    license='MIT',
    keywords='nmr parsing jcamp dx',
    packages=find_packages(),
    install_requires=['numpy', 'scipy', 'pandas', 'dataclasses', ],
    python_requires='>=3.6',
    include_package_data=True,
)

