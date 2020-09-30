FROM continuumio/miniconda3

# Synced development folder
# COPY . /map2loop-2

# Install deps for compiling m2m
RUN apt update && apt install -y build-essential mlocate git cmake vim python3-dev python3-vtk7 xvfb 
RUN updatedb

RUN git clone https://gist.github.com/yohanderose/083a04767328de71128b542d300e75dc vimstuff
RUN cp vimstuff/.vimrc /etc/vim/vimrc

RUN git clone https://github.com/Loop3D/map2loop-2

# Create m2l conda environment:
RUN conda env create -f /map2loop-2/environment.yml
ENV PATH /opt/conda/envs/m2l/bin:$PATH
ENV CONDA_DEFAULT_ENV m2l
RUN /bin/bash -c "source activate m2l"

# Install new package
RUN pip install pybind11 pytest jupyter
RUN pip install -e map2loop-2

# # Fetch and install model engines
# > Structural
RUN pip install lavavu
RUN git clone https://github.com/Loop3D/LoopStructural
RUN pip install -r /LoopStructural/requirements.txt
RUN pip install -e /LoopStructural 
# # > Gempy
# RUN pip install -r /map2loop-2/engines/gempy-requirements.txt
# RUN pip install gempy
# # > PyNoddy
# RUN git clone https://github.com/cgre-aachen/pynoddy
# RUN cp /map2loop-2/engines/noddy /usr/bin
# RUN pip install pynoddy

# # Fetch and install ensemble generator
# # RUN git clone --branch docker --single-branch https://Loop3D:2fae576a7fb2b205dddfc7a8004a694812623142@github.com/Loop3D/ensemble_generator
# # RUN pip install -e /ensemble_generator


# Add Tini
ENV TINI_VERSION v0.18.0
ADD https://github.com/krallin/tini/releases/download/${TINI_VERSION}/tini /tini
RUN chmod +x /tini
# script for xvfb-run.  all docker commands will effectively run under this via the entrypoint
RUN printf "#\041/bin/sh \n rm -f /tmp/.X99-lock && xvfb-run -s '-screen 0 1600x1200x16' \$@" >> /usr/local/bin/xvfbrun.sh && \
    chmod +x /usr/local/bin/xvfbrun.sh
# note we use xvfb which to mimic the X display for lavavu
ENTRYPOINT ["/tini", "--", "/usr/local/bin/xvfbrun.sh"]

# Fetch, install and setup original repo
# RUN git clone https://github.com/Loop3D/map2loop
# # Build map2model from source
# RUN mkdir /map2loop/m2m_source/build 
# RUN cd /map2loop/m2m_source/build && cmake .. && make -j2 && cp map2model ../../m2m_cpp
# # Install deps and then the package itself
# RUN conda install -c conda-forge ipywidgets
# RUN conda install -c conda-forge ipyleaflet
# RUN conda install -c conda-forge folium
# RUN pip install ipyfilechooser
# RUN jupyter nbextension enable --py --sys-prefix ipyleaflet
# RUN pip install -e /map2loop

# Execute jupyter on run 
CMD ["jupyter","notebook","--ip=0.0.0.0", "--allow-root", "--no-browser", "--NotebookApp.password=''"]

EXPOSE 8888