import pathlib
from setuptools import setup

# The directory containing this file
HERE = pathlib.Path(__file__).parent

# The text of the README file
README = "README.MD"


setup(
  name = 'LammpsFileManipulation',
  packages = ['LammpsFileManipulation'],
  version = '0.13',
  license='MIT',
  description = 'This is a package designed to help streamline the process of preprocessing LAMMPS output files for scientific calculations/manipulations in Python. The class structures are built using pandas DataFrames making it easy to manipulate.',
  long_discription = README,
  long_discription_content_type = "text/markdown",
  author = 'Aaron Schwan',
  author_email = 'schwanaaron@gmail.com',
  url = 'https://github.com/AaronSchwan/LammpsFileManipulation',
  download_url = 'https://github.com/AaronSchwan/LammpsFileManipulation/releases/tag/v_0.13-alpha',
  keywords = ['LAMMPS','atomistic',"dump file"],
  install_requires=[
          'pandas',
          'numpy',
          'sys',
          'os',
          'gc',
          'ntpath',
          'pickle',
          'time'],
  classifiers=[
    'Development Status :: 3 - Alpha',
    'Intended Audience :: Developers',
    'Topic :: Software Development :: Build Tools',
    'License :: OSI Approved :: MIT License',
    'Programming Language :: Python :: 3',
    'Programming Language :: Python :: 3.6',
    'Programming Language :: Python :: 3.7',
    'Programming Language :: Python :: 3.8',
  ],
)
