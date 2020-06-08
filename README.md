# Map2Loop 2.0

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

Download and install the docker containerisation software and CLI [https://docs.docker.com/get-docker/](https://docs.docker.com/get-docker/)

1. Clone this repo and navigate inside using the recurse submodules flag to fetch the example data. 

    ```bash
    git clone --recurse-submodules https://github.com/Loop3D/map2loop-2
    ```

2. Run the following command and click on the link Jupyter outputs to access the original [map2loop](https://github.com/Loop3D/map2loop) notebooks.

    ```bash
    docker-compose --compatibility up
    ```

3. To jump into a bash shell in the container itself, open a new terminal and issue the following command. 

    ```bash
    docker exec -it map2loop-2_dev_1 bash
    ```

## Install via PyPi

Still to come...

### Known Issues

- Developing using docker on a Windows host will mean you will not have GPU passthrough funtionality - so can't use a discrete graphics card in the container even if you have one.

### References