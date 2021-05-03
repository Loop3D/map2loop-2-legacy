import os
import sys
import time
import setuptools
import importlib
from setuptools.command.develop import develop
import subprocess
import platform
from map2loop import __version__

head, tail = os.path.split(sys.argv[0])

try:
    # proc = subprocess.Popen(
    # 'git tag'.split(), stdin=subprocess.PIPE, stdout=subprocess.PIPE)
    # existing_tags = proc.communicate()[0]
    # tag_clean = 'git tag -d {}'.format(existing_tags.decode('ascii'))
    # subprocess.run(
    # tag_clean.split())
    tag_create = 'git tag -a {0} -m {0}'.format(__version__)
    subprocess.run(
        tag_create.split())
except Exception as e:
    print(e)


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


def get_description():
    long_description = ""
    readme_file = os.path.join(head, "README.md")
    with open(readme_file, "r") as fh:
        long_description = fh.read()
    return long_description


setuptools.setup(
    name="map2loop",
    version=__version__,
    author="The Loop Organisation",
    author_email="yohan.derose@monash.edu",
    description="Generate 3D model data from 2D maps.",
    long_description=get_description(),
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
    install_requires=[
        'setuptools_scm',
    ],
    include_package_data=True,
)
