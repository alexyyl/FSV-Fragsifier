from setuptools import setup, find_packages
import os
 
__version__ = '1.0.1'

# Read requirements from requirements.txt
with open(os.getcwd()+'/fragsifier/requirements.txt') as f:
    requirements = f.read().splitlines()
    
setup(
  name = 'fragsifier',     
  packages = find_packages(),  
  version = __version__,      
  license='BSD 3-Clause',    
  description = 'forensic STR sequence extraction tool',
  author = 'Alexander YY Liu',
  author_email = 'yliu575@aucklanduni.ac.nz',  
  url = 'https://github.com/alexyyl/FSV-Fragsifier', 
  download_url = 'https://github.com/alexyyl/FSV-Fragsifier/releases/tag/v_%s.tar.gz' % __version__, 
  keywords = ['forensic STR', 'sequence extraction'], 
  python_requires='>=3.6',
  install_requires=requirements,
  classifiers=[
    'Development Status :: 4 - Beta',
    'Intended Audience :: Science/Research',
    'Topic :: Scientific/Engineering :: Bio-Informatics',
    'License :: OSI Approved :: BSD License',  
    'Programming Language :: Python :: 3.6',
  ],
  entry_points={
    'console_scripts': ["fragsifier=fragsifier.fragsifier:main"]
  },
  package_data = {'': ['', '*', '*/*']},
)
