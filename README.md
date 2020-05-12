# README

Generate model input data from geological maps. Revision of objectives in [https://github.com/Loop3D/map2loop](https://github.com/Loop3D/map2loop)

## Dependencies

If you wish to, create your own python virtual environment with the following modules to run the map2loop examples.

- python=3.7
- numpy
- pandas
- geopandas
- matplotlib
- rasterio
- networkx
- owslib
- pyamg
- descartes
- mplstereonet

## Build with Docker

Download and install the docker containerisation framework and CLI [https://docs.docker.com/get-docker/](https://docs.docker.com/get-docker/)

Clone this repo and navigate inside. Run the following to build the 'm2l' container.

```bash
docker build -t m2l .
```

Then execute the following  to run the example notebooks.

```bash
docker run -it -p 8888:8888 m2l
```

## Install via PyPi

Still to come...

### Bugs

### References