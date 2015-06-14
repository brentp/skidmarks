from setuptools import setup

version="0.0.5"

setup(name='skidmarks',
      version=version,
      py_modules=['skidmarks'],
      description="find runs (non-randomness) in sequences",
      url="http://github.com/brentp/skidmarks/",
      long_description=open('README.md').read(),
      keywords='bioinformatics sequence randomness test',
      author='brentp',
      author_email='bpederse@gmail.com',
      license='MIT',
      test_suite='skidmarks.test_suite',
      zip_safe=False,
      install_requires=[
          # -*- Extra requirements: -*-
      ],
      scripts=[],
      entry_points={
      },
      classifiers=[
          "Programming Language :: Python :: 3",
          "Programming Language :: Python :: 2",
          "Topic :: Scientific/Engineering :: Bio-Informatics",
          ],
)


