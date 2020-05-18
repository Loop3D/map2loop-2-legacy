FROM continuumio/miniconda3

# Install deps for compiling m2m
RUN apt update && apt install -y build-essential git cmake

# Create m2l conda environment:
COPY environment.yml .
RUN conda env create -f environment.yml
ENV PATH /opt/conda/envs/m2l/bin:$PATH
RUN /bin/bash -c "source activate m2l"

# Fetch, install and setup original repo
RUN git clone https://github.com/Loop3D/map2loop
RUN pip install /map2loop

# Build map2model from source
ADD m2m /
RUN /bin/bash -c "chmod +x build-m2m.sh"
RUN /bin/bash -c "chmod -R 777 /m2m_source"
RUN ./build-m2m.sh


CMD ["jupyter","notebook","map2loop","--ip=0.0.0.0", "--allow-root", "--no-browser"]