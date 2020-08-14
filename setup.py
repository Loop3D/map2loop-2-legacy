import setuptools

long_description = ""
with open("README.md", "r") as fh:
    long_description = fh.read()


setuptools.setup(
    name="map2loop",
    version="1.0.6",
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
    install_requires=[
        'map2model-loop3d',
        'numpy',
        'pandas',
        'geopandas',
        'matplotlib',
        'hjson',
        'networkx',
        'shapely',
        'rasterio',
        'scipy'
    ],
    python_requires='>=3.6',
)
