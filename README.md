# Map2Loop 2.0

Generate 3D geological model inputs from geographical maps — a high-level implementation and extension of [https://github.com/Loop3D/map2loop](https://github.com/Loop3D/map2loop) by Prof. Mark Jessell at UWA. To see an example interactive model build using map2loop and LoopStructural, follow this link: <a href="http://geo.loop-gis.org/models/vtkleaflet_2021-03-19-11-40.html">3D Model from the Hamersley region, Western Australia</a>)

## Install

You will need some flavour of conda (a python package manager, [see here](https://docs.anaconda.com/anaconda/install/index.html)), as well as Python ≥ 3.6

### Run 🐍

To just use map2loop, issue the following

```bash
conda install -c conda-forge -c loop3d map2loop -y
```

### Development 🔧

If you want to tinker yourself/contribute, clone the source code with

```bash
git clone https://github.com/Loop3D/map2loop-2.git
```

Or get the source + example notebooks with

```bash
git clone https://github.com/Loop3D/map2loop-2.git
git clone --single-branch --branch yohan https://github.com/Loop3D/map2loop2-notebooks
```

Navigate into map2loop-2, and issue the following to install map2loop and its dependencies. _Note_: The 'develop' flag makes your source changes take effect on saving, so you only need to run this once

```bash
python setup.py develop
```

---

## Building with Docker

Fair warning, we recommend conda to almost everyone. With great software development power comes great environment setup inconvenience. You'll need to download and install the [docker containerisation software](https://docs.docker.com/get-docker/), and the docker and docker-compose CLI.

### Development ️️️️🐋

1. Clone this repo and navigate inside as per above
2. Run the following and click on the Jupyter server forwarded link to access and edit the notebooks

   ```bash
   docker-compose up --build
   ```

3. To hop into a bash shell in a running container, open a terminal and issue

   ```bash
   docker ps
   ```

   Find the container name or ID and then run

   ```bash
   docker exec -it <container_NAMEorID> bash
   # Probably -> docker exec -it  map2loop-2_dev_1 bash
   ```

---

## Usage

Our notebooks cover use cases in more detail, but here is an example of processing Loop's South Australia remote geospatial data in just 20 lines of Python.

First, lets import map2loop and define a bounding box. You can use GIS software to find one or use [Loop's Graphical User Interface](https://loop3d.github.io/downloads.html) for the best experience and complete toolset. Remember what projection your coordinates are in!

```python
from map2loop.project import Project

bbox_3d = {
    'minx': 250805.1529856466,
    'miny': 6405084.328058686,
    'maxx': 336682.921539395,
    'maxy': 6458336.085975628,
    'base': -3200,
    'top': 1200
}
```

![sa example](docs/Untitled.png?raw=true)

Then, specify: the state, directory for the output, the bounding box and projection from above - and hit go! That's it.

```python
proj = Project(loopdata_state="SA")

proj.update_config(out_dir='sa-remote',
                   overwrite="true",
                   bbox_3d=bbox_3d,
                   proj_crs={'init': 'EPSG:28354'},
                #    drift_prefix=['T', 'Q', 'water', 'void'],
                #    quiet='no-figures',
                   )
proj.run()
```

This is a minimal example and a small part of Loop.

Our _documentation and other resources outline how to extend map2loop and port to the LoopStructural modelling engine. We are working to incorporate geophysical tools and best provide visualisation and workflow consolidation in the GUI._

_Loop is led by Laurent Ailleres (Monash University) with a team of Work Package leaders from:_

- _Monash University: Roy Thomson, Lachlan Grose and Robin Armit_
- _University of Western Australia: Mark Jessell, Jeremie Giraud, Mark Lindsay and Guillaume Pirot_
- _Geological Survey of Canada: Boyan Brodaric and Eric de Kemp_

---

### Known Issues and FAQs

- Developing with docker on Windows means you won't have GPU passthrough and can’t use a discrete graphics card in the container even if you have one.
- If Jupyter links require a token or password, it may mean port 8888 is already in use. To fix, either make docker map to another port on the host ie -p 8889:8888 or stop any other instances on 8888.
- Sometimes the submodules misbehave. Ensure you specify the recurse-submodules flag with an 's' as opposed to recurse-submodule, and double-check your .git directory is clean, see [https://github.com/Loop3D/map2loop-2/issues/41](https://github.com/Loop3D/map2loop-2/issues/41)

### Links

[https://loop3d.github.io/](https://loop3d.github.io/)

[https://github.com/Loop3D/LoopStructural](https://github.com/Loop3D/LoopStructural)
