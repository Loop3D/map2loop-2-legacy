import sys
import setuptools
from setuptools.command.install import install


class CondaDependencies(install):
    user_options = install.user_options + [
        ('conda', None, None),  # a 'flag' option
        # ('someval=', None, None) # an option that takes a value
    ]

    def initialize_options(self):
        install.initialize_options(self)
        self.conda = None
        #self.someval = None

    def finalize_options(self):
        print("Using conda for installing dependencies =", self.conda)
        if self.conda:
            import subprocess

            deps = []
            with open('./dependencies.txt', 'r') as f:
                for line in f:
                    deps.append(line.strip())
            
            command = 'conda install -c loop3d -y python=3.7'.split() + deps
            subprocess.run(command)

        install.finalize_options(self)

    def run(self):
        install.run(self)


long_description = ""
with open("README.md", "r") as fh:
    long_description = fh.read()


setuptools.setup(
    name="map2loop",
    version="1.0.8",
    author="Yohan de Rose",
    author_email="contact@loop3d.org",
    description="Generate 3D model data using 2D maps.",
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
        'install': CondaDependencies,
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
