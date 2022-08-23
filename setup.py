import os
import sys
import setuptools
import subprocess
from map2loop import __version__

head, tail = os.path.split(sys.argv[0])

# try:
#     cmd = 'bash .git.sh'
#     subprocess.run(
#         cmd.split())
#     cmd = 'git tag -a {0} -m {0}'.format(
#         __version__)
#     subprocess.run(
#         cmd.split())
# except Exception as e:
#     print(e)

# def get_description():
#     long_description = ""
#     readme_file = os.path.join(head, "README.md")
#     with open(readme_file, "r") as fh:
#         long_description = fh.read()
#     return long_description


setuptools.setup(
    name="map2loop",
    version=__version__,
    author="The Loop Organisation",
    author_email="roy.thomson@monash.edu",
    description="Generate 3D model data from 2D maps.",
    # long_description=get_description(),
    long_description_content_type="text/markdown",
    # dependency_links deprecated and not recommended
    # dependency_links=["https://www.lfd.uci.edu/~gohlke/pythonlibs/rasterio‑1.2.10‑cp39‑cp39‑win_amd64.whl; platform_system=='Windows'"],
    # install_requires=[
        # "numpy",
        # "gdal==3.5.1",
        # "pandas",
        # "geopandas",
        # "pyproj",
        # "shapely",
        # "rasterio==1.2.10",
        # "descartes",
        # "matplotlib",
        # "networkx",
        # "scipy",
        # "owslib",
        # "tqdm",
        # "mplstereonet",
        # "hjson",
        # "loopprojectfile",
        # "beartype",
        # "cftime>=1.4.0",
        # "map2model-loop3d"
    # ],
    url="https://github.com/Loop3D/map2loop-2",
    packages=setuptools.find_packages(),
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "Natural Language :: English",
        "Programming Language :: Python",
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Topic :: Scientific/Engineering",
        "Topic :: Scientific/Engineering :: GIS",
    ],
    python_requires=">=3.7",
)
