import setuptools
from setuptools.command.install import install
import sys
from distutils.command.sdist import sdist as sdist_orig
from distutils.errors import DistutilsExecError

class CondaDependencies(sdist_orig):
    def run(self):
        try:
            print("Installing dependencies with conda...")
            deplist_path = "./dependencies.txt"
            deps = []
            with open(deplist_path, 'r') as f:
                for line in f:
                    deps.append(line.strip())

            command = 'conda install -c loop3d -c conda-forge -y python=3.7'.split() + deps
            self.spawn(command)
        except DistutilsExecError:
            self.warn('WARNING: Could not install dependencies using conda!')
        super().run()


long_description = ""
with open("README.md", "r") as fh:
    long_description = fh.read()


setuptools.setup(
    name="map2loop",
    version="1.0.8",
    author="Yohan de Rose",
    author_email="contact@loop3d.org",
    description="Generate 3D model data from 2D maps.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/Loop3D/map2loop-2",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    cmdclass={
        'sdist': CondaDependencies,
    },
    # install_requires=[
    #     'map2model-loop3d',
    #     'numpy',
    #     'pandas',
    #     'geopandas',
    #     'matplotlib',
    #     'hjson',
    #     'networkx',
    #     'shapely',
    #     'rasterio',
    #     'scipy'
    # ],
    python_requires='>=3.6',
)
