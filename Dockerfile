FROM continuumio/miniconda3

COPY . map2loop-2/

# Tools for developing in the contiainer
RUN apt update && apt install -y build-essential \
	mlocate git cmake vim python3 python3-dev python3-vtk7 python3-numpy \
	xvfb libglew-dev zlib1g zlib1g-dev libgl1-mesa-dev libx11-dev libtiff-dev libtiff-tools
RUN updatedb && ldconfig

## > M2l
RUN conda install python=3.7 setuptools jupyter pybind11 pytest -y
RUN cd map2loop-2 && python setup.py develop

## > Structural
RUN git clone https://github.com/Loop3D/LoopStructural
RUN pip install -r /LoopStructural/requirements.txt
RUN pip install lavavu
RUN pip install -e /LoopStructural
# # > Gempy
# RUN pip install -r /map2loop-2/engines/gempy-requirements.txt
# RUN pip install gempy
# # > PyNoddy
# RUN git clone https://github.com/cgre-aachen/pynoddy
# RUN cp /map2loop-2/engines/noddy /usr/bin
# RUN pip install pynoddy

# > Ensemble Generator
# RUN git clone --branch docker --single-branch https://Loop3D:2fae576a7fb2b205dddfc7a8004a694812623142@github.com/Loop3D/ensemble_generator
# RUN pip install -e /ensemble_generator

# > Tini
ENV TINI_VERSION v0.18.0
ADD https://github.com/krallin/tini/releases/download/${TINI_VERSION}/tini /tini
RUN chmod +x /tini
RUN printf "#\041/bin/sh \n rm -f /tmp/.X99-lock && xvfb-run -s '-screen 0 1600x1200x16' \$@" >> /usr/local/bin/xvfbrun.sh && \
	chmod +x /usr/local/bin/xvfbrun.sh
#ENTRYPOINT ["/tini", "--", "/usr/local/bin/xvfbrun.sh"]

# > Leaflet Map Notebooks
RUN conda install -c conda-forge ipywidgets
RUN conda install -c conda-forge ipyleaflet
RUN conda install -c conda-forge folium
RUN pip install ipyfilechooser
RUN jupyter nbextension enable --py --sys-prefix ipyleaflet

# > Demo Notebooks
RUN git clone --single-branch --branch yohan https://github.com/Loop3D/map2loop2-notebooks
RUN conda update --all -y

EXPOSE 8888
EXPOSE 9999
