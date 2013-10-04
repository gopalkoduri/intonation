from setuptools import setup

def readme():
    with open('README.rst') as f:
        return f.read()

setup(name='pypeaks',
      version='0.1',
      description='Python module to analyze, study and research the intonation, requires=[', requires=['scipy', 'numpy']numpy']
      of various intervals in different music traditions around the world',
      long_description=readme(),
      classifiers=[
        'Development Status :: 4 - Beta',
        'License :: OSI Approved :: GNU Affero General Public License v3',
        'Programming Language :: Python :: 2.7',
        'Topic :: Scientific/Engineering :: Information Analysis',
        'Topic :: Software Development :: Libraries :: Python Modules',
      ],
      keywords='python intonation research music intervals tuning',
      url='https://github.com/gopalkoduri/intonation',
      author='Gopala Krishna Koduri',
      author_email='gopala.koduri@gmail.com',
      license='GNU Affero GPL v3',
      packages=['intonation'],
      #data_files=[('examples', ['howto.ipynb', 'howto.py', 'ji-intervals.pickle' 'sample-histogram.pickle'])],
      install_requires=[
          'numpy',
          'matplotlib',
      ],
      zip_safe=False)
