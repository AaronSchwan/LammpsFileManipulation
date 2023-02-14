import pathlib
from setuptools import setup

# The directory containing this file
HERE = pathlib.Path(__file__).parent

# The text of the README file
README = "README.md"
with open(README,'r') as fh:
    long_description = fh.read()

setup(
  name = 'LammpsFileManipulation',
  packages = ['LammpsFileManipulation'],
  version = '1.0.1',
  license='MIT',
  description = 'This is a package designed to help streamline the process of preprocessing LAMMPS output files for scientific calculations/manipulations in Python. The class structures are built using pandas DataFrames making it easy to manipulate.',
  long_description = long_description,
  long_description_content_type = "text/markdown",
  author = 'Aaron Schwan',
  author_email = 'schwanaaron@gmail.com',
  url = 'https://github.com/AaronSchwan/LammpsFileManipulation',
  download_url = 'https://github.com/AaronSchwan/LammpsFileManipulation/releases/tag/v_0.17-alpha',
  keywords = ['LAMMPS','atomistic',"dump file"],
  install_requires=[
          'pandas',
          'numpy'],
  classifiers=[
    'Development Status :: 3 - Alpha',
    'Intended Audience :: Developers',
    'Topic :: Software Development :: Build Tools',
    'License :: OSI Approved :: MIT License',
    'Programming Language :: Python :: 3'],
)
