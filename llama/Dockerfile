FROM ubuntu:20.04

ENV LANG=C.UTF-8 LC_ALL=C.UTF-8
ENV PATH /opt/conda/bin:$PATH

RUN apt-get update --fix-missing && \
    apt-get install -y wget bzip2 ca-certificates curl git unzip && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

RUN wget --quiet https://repo.anaconda.com/miniconda/Miniconda2-latest-Linux-x86_64.sh -O ~/miniconda.sh && \
    /bin/bash ~/miniconda.sh -b -p /opt/conda && \
    rm ~/miniconda.sh && \
    /opt/conda/bin/conda clean -tipsy && \
    ln -s /opt/conda/etc/profile.d/conda.sh /etc/profile.d/conda.sh && \
    echo ". /opt/conda/etc/profile.d/conda.sh" >> ~/.bashrc && \
    echo "conda activate base" >> ~/.bashrc

ENV TINI_VERSION v0.16.1
ADD https://github.com/krallin/tini/releases/download/${TINI_VERSION}/tini /usr/bin/tini
RUN chmod +x /usr/bin/tini


RUN mkdir /tools && \
    cd /tools && \
    wget --quiet https://github.com/cov-lineages/llama/archive/master.zip && \
    unzip master.zip && \
    cd llama-master/ && \
    conda env create -f environment.yml

RUN conda init bash && \
    . ~/.bashrc && \
    conda create --name llama && \
    conda activate llama && \
    cd /tools/llama-master/ && \
    python setup.py install
   
#ENTRYPOINT [ "/usr/bin/tini", "--" ]
#CMD [ "/bin/bash" ]
#ENTRYPOINT ["conda activate llama"]

RUN echo "source activate llama" > ~/.bashrc
ENV PATH /opt/conda/envs/env/bin:$PATH
