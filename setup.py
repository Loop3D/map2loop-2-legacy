import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="map2loop-2",
    version="0.0.1",
    author="Yohan de Rose",
    author_email="yohan.derose@monash.edu",
    description="An example package hosting for m2l.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/Loop3D/map2loop-2",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    # install_requires=[
    #     'GDAL',
    #     'numpy',
    #     'pandas',
    #     'geopandas',
    #     'shapely'
    # ],
    python_requires='>=3.6',
)
