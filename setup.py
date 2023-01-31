from setuptools import setup, find_packages
from distutils.util import convert_path

setup(name = "gdr3binaryorbits",
    version = 0.1,
    description = "Visualize RV orbits from Gaia DR3",
    author = "Tharindu Jayasinghe",
    author_email = "",
    url = "https://github.com/tjayasinghe/gdr3binaryorbits",
    packages = find_packages(include=['gdr3binaryorbits', 'gdr3binaryorbits.*']),
    classifiers=[
      'Intended Audience :: Science/Research',
      'Operating System :: OS Independent',
      'Programming Language :: Python :: 3',
      'License :: OSI Approved :: MIT License',
      'Topic :: Scientific/Engineering :: Astronomy'
      ],
    zip_safe=False,
    python_requires=">=3.6",
    install_requires=["numpy", "scipy", "matplotlib", "pandas", "astropy"]
)