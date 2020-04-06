from setuptools import setup
import os
import re


def read_version():
    path = os.path.join(os.path.abspath(os.path.dirname(__file__)), 'sciutil/__init__.py')
    with open(path, 'r') as fh:
        return re.search(r'__version__\s?=\s?[\'"](.+)[\'"]', fh.read()).group(1)


def readme():
    with open('README.md') as f:
        return f.read()


setup(name='sciutil',
      version=read_version(),
      description='Utility functions for the sci* packages.',
      long_description=readme(),
      long_description_content_type='text/markdown',
      author='Ariane Mora',
      author_email='ariane.mora@uq.net.au',
      url='https://github.com/ArianeMora/sciutil',
      license='GPL3',
      project_urls={
          "Bug Tracker": "https://github.com/ArianeMora/sciutil/issues",
          "Documentation": "https://github.com/ArianeMora/sciutil",
          "Source Code": "https://github.com/ArianeMora/sciutil",
      },
      classifiers=[
          'Development Status :: 5 - Production/Stable',
          'Intended Audience :: Science/Research',
          'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
          'Natural Language :: English',
          'Operating System :: OS Independent',
          'Programming Language :: Python :: 3.6',
          'Programming Language :: Python :: 3.7',
          'Programming Language :: Python :: 3.8',
          'Topic :: Scientific/Engineering :: Bio-Informatics',
      ],
      keywords='util',
      packages=['sciutil'],
      entry_points={
          'console_scripts': [
              'sciutil = sciutil.__main__:main'
          ]
      },
      install_requires=[],
      python_requires='>=3.6',
      data_files=[("", ["LICENSE"])]
      )