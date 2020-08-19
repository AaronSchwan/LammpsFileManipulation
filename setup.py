from distutils.core import setup
setup(
  name = 'LammpsFileManipulation',
  packages = ['LammpsFileManipulation'],
  version = '0.11',
  license='MIT',
  description = 'This is a package designed to help streamline the process of preprocessing LAMMPS output files for scientific calculations/manipulations in Python. The class structures are built using pandas DataFrames making it easy to manipulate.',
  author = 'Aaron Schwan',
  author_email = 'schwanaaron@gmail.com',
  url = 'https://github.com/AaronSchwan/LammpsFileManipulation',
  download_url = 'https://github.com/AaronSchwan/LammpsFileManipulation/releases/tag/v_0.1-alpha',
  keywords = ['LAMMPS','atomistic',"dump file"],
  install_requires=[
          'pandas',
          'numpy'],
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
