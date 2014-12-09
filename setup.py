from setuptools import setup

version="0.0.4"

setup(name='skidmarks',
      version=version,
      py_modules=['skidmarks'],
      description="find runs (non-randomness) in sequences",
      url="http://github.com/brentp/skidmarks/",
      long_description=open('README.md').read(),
      classifiers=["Topic :: Scientific/Engineering :: Bio-Informatics"],
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
)


