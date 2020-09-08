# Map2Loop 2.0

Generate model input data from geological maps. High-level implementation and expansion of [https://github.com/Loop3D/map2loop](https://github.com/Loop3D/map2loop)

## Build with Conda

### Usage

If you just want to use map2loop, issue the following command.

```bash
conda install -c loop3d map2loop -y
```

### Development

If you want to develop or contribute, clone this repo and navigate into it. Issue the following command and all dependencies, as well as map2loop, will be installed in the currently active conda environment.

```bash
python dependencies.py && pip install -e .
```

## Build with Docker

Download and install the docker containerisation software and CLI [https://docs.docker.com/get-docker/](https://docs.docker.com/get-docker/)

### Usage

If you do not wish to develop, you can pull the complete docker image with map2loop and LoopStructural preinstalled, and then run the container to interact with the notebooks.

```bash
docker pull yohanderose/map2loop-2_dev:0.0.2
docker run -it -p 8888:8888 yohanderose/map2loop-2_dev:0.0.2
```

### Development

1. Clone this repo and navigate inside using the recurse submodules flag to fetch the example data. 

    ```bash
    git clone --recurse-submodules https://github.com/Loop3D/map2loop-2
    ```

2. Run the following command and click on the link Jupyter outputs to access the original [map2loop](https://github.com/Loop3D/map2loop) notebooks.

    ```bash
    docker-compose up
    ```

3. To jump into a bash shell in the container itself, open a new terminal and issue the following command. 

    ```bash
    docker exec -it map2loop-2_dev_1 bash
    ```

### Known Issues

- Developing using docker on a Windows host will mean you will not have GPU passthrough and can't use a discrete graphics card in the container even if you have one.
- If Jupyter pops up and requires a token or password, it likely means port 8888 is already in use. To fix, either make docker map to another port on the host ie -p 8889:8888 or stop any other instances on 8888.

### References