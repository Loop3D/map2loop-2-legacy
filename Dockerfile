FROM continuumio/miniconda3

# Install deps for compiling m2m
RUN apt update && apt install -y build-essential git cmake

# Create m2l the environment:
COPY environment.yml .
RUN conda env create -f environment.yml
ENV PATH /opt/conda/envs/m2l/bin:$PATH
RUN /bin/bash -c "source activate m2l"

COPY m2m .
RUN /bin/bash -c "chmod -R 777 /m2m_source"
# Execute m2m build script
RUN ./build-m2m.sh
