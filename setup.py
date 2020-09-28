import setuptools
from setuptools.command.develop import develop
import subprocess


class CondaDependencies(develop):
    def run(self):
        try:
            deplist_path = "./dependencies.txt"
            deps = []
            with open(deplist_path, 'r') as f:
                for line in f:
                    deps.append(line.strip())

            command = 'conda install -c defaults -c conda-forge -y python=3.7'.split() + \
                deps
            print(command)
            subprocess.call(command, shell=True)
            command = 'conda install -c loop3d -y mplstereonet'.split()
            print(command)
            subprocess.call(command, shell=True)
            command = 'conda install -c loop3d -y hjson'.split()
            print(command)
            subprocess.call(command, shell=True)
            command = 'conda install -c loop3d -y map2model'.split()
            print(command)
            subprocess.call(command, shell=True)
        except Exception as e:
            self.error('Could not install dependencies using conda!')

        develop.run(self)


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
        'develop': CondaDependencies,
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
