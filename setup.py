from setuptools import setup, find_packages

setup(
    name='ngs-pipeline-launcher',
    version='1.0.0',
    packages=find_packages(exclude=['tests*']),
    description='A launcher tool for sorting a sequencing run by pathogen and automatic initialization of their downstream pipelines.',
    url='https://github.com/provlab-bioinfo/ngs-pipeline-launcher',
    author='Andrew Lindsay',
    author_email='andrew.lindsay@albertaprecisionlabs.ca',
    include_package_data=True,
    keywords=[],
    zip_safe=False,
    install_requires=[
        'matplotlib',
        'numpy',
        'pandas',
        'openpyxl',
        'pip',
        'build',
        'search-tools @ git+https://github.com/provlab-bioinfo/search-tools'
    ],
    python_requires='>=3.10, <4'
)