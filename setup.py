from setuptools import setup, find_packages

with open("README.md", "r") as fh:
    long_description = fh.read()

classifiers = [
    'Development Status :: 4 - Beta',
    'Environment :: Console',
    'Intended Audience :: Science/Research',
    'License :: OSI Approved :: MIT License',
    'Programming Language :: Python :: 3',
    'Programming Language :: Python :: 3.2',
    'Programming Language :: Python :: 3.3',
    'Programming Language :: Python :: 3.4',
    'Programming Language :: Python :: 3.5',
    'Programming Language :: Python :: 3.6',
    'Programming Language :: Python :: 3.7',
    'Programming Language :: Python :: 3.8',
    'Topic :: Scientific/Engineering :: Bio-Informatics',
    'Operating System :: POSIX',
    'Operating System :: Unix',
    'Operating System :: MacOS'
]

setup(name='gwaspy',
      version='0.1.0',
      author='Lindokuhle Nkambule',
      author_email='lnkambul@broadinstitute.org',
      url='https://gwaspy.readthedocs.io/',
      project_urls={"GitHub": "https://github.com/atgu/GWASpy"},
      description='GWASpy: A Python package for performing GWAS QC, PCA, phasing, and genotype imputation.',
      long_description=long_description,
      long_description_content_type="text/markdown",
      license='MIT',
      packages=find_packages(),
      entry_points={
          'console_scripts': [
              'preimp_qc = gwaspy.preimp_qc.preimp_qc:main',
              'pca = gwaspy.pca.pca:main',
              'imputation = gwaspy.imputation.imputation:main',
              'phasing = gwaspy.phasing.phasing:main'
          ]
      },
      classifiers=classifiers,
      keywords='',
      install_requires=['hail', 'matplotlib', 'numpy', 'pandas', 'pylatex', 'plotly'],
      zip_safe=False
      )