FROM continuumio/miniconda3

# Synced development folder
COPY . /map2loop-2

# Install deps for compiling m2m
RUN apt update && apt install -y build-essential git cmake python3-vtk7

# Create m2l conda environment:
RUN conda env create -f /map2loop-2/environment.yml
ENV PATH /opt/conda/envs/m2l/bin:$PATH
RUN /bin/bash -c "source activate m2l"

# Fetch, install and setup original repo
RUN git clone https://github.com/Loop3D/map2loop
RUN pip install /map2loop

# Fetch and install model engines
# RUN git clone https://github.com/Loop3D/LoopStructural
# RUN pip install -r /LoopStructural/requirements.txt
# RUN pip install /LoopStructural
RUN pip install -r /map2loop-2/gempy-requirements.txt
RUN pip install gempy

# Build map2model from source
ADD maps/m2m /
RUN /bin/bash -c "chmod +x build-m2m.sh"
RUN /bin/bash -c "chmod -R 777 /m2m_source"
RUN ./build-m2m.sh

# Execute jupyter on run 
CMD ["jupyter","notebook","map2loop","--ip=0.0.0.0", "--allow-root", "--no-browser"]

EXPOSE 8888