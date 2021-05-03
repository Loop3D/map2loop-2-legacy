import os
import sys
import setuptools
from setuptools.command.develop import develop
import subprocess
import platform

head, tail = os.path.split(sys.argv[0])


class CondaDependencies(develop):
    def run(self):
        try:
            deplist_path = os.path.join(head, "dependencies.txt")
            deps = []
            with open(deplist_path, 'r') as f:
                for line in f:
                    deps.append(line.strip())

            _platform = platform.platform()

            if _platform.startswith("Windows"):
                _shell = True
            else:  # Linux or Mac
                _shell = False

            command = 'conda install -c conda-forge -c loop3d -y python=3.7'.split() + \
                deps
            print(command)
            subprocess.run(command, shell=_shell)
        except Exception as e:
            self.error('Could not install dependencies using conda!')

        develop.run(self)


long_description = ""
readme_file = os.path.join(head, "README.md")
with open(readme_file, "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="map2loop",
    version="1.1.5",
    author="The Loop Organisation",
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
    python_requires='>=3.6',
)
