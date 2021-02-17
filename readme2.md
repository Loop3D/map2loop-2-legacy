# Map2Loop 2.0

From an Anaconda prompt, execute the following commands:

`git clone https://github.com/Loop3D/LoopStructural`

`git clone -b mj2021 https://github.com/Loop3D/map2loop-2`


If {map2loop-2} repository already exists,

`cd map2loop-2` 

and execute `git clone -b mj2021 https://github.com/Loop3D/map2loop-2`

Next:

`git clone https://github.com/Loop3D/loopprojectfile`

`git clone https://github.com/Loop3D/map2loop2-notebooks`

`conda create --name loop2`

`conda activate loop2`

`cd LoopStructural`

`python setup.py install build_ext --inplace`

`cd ../loopprojectfile`

`python setup.py install`

`cd ../map2loop-2`

`python setup.py develop`

`conda install meshio pyamg -c conda-forge -y`

`pip install lavavu`

you may need to install cython `pip install cython`

may need  `conda install jupyter`

may need `conda install ipywidgets`

may need `conda install -c conda-forge gdal`


if maps don't appear in notebook 2 uncomment first three lines of first cell and rerun and then restart jupyter server