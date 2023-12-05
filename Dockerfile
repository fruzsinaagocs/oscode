FROM ubuntu:18.04

RUN apt-get update && apt-get install -y \
    git python3-pip
RUN pip3 install --no-cache --upgrade pip && \
    pip3 install --no-cache notebook && \
    pip3 install scipy numpy matplotlib ipywidgets jupyter_contrib_nbextensions \
    hide_code cite2c RISE && \
    jupyter contrib nbextension install --user && \
    jupyter nbextension enable hide_code --user --py && \
    python3 -m cite2c.install

#RUN git clone https://github.com/eigenteam/eigen-git-mirror
#ENV CPLUS_INCLUDE_PATH=${PWD}/eigen-git-mirror/:${CPLUS_INCLUDE_PATH}
#RUN echo $CPLUS_INCLUDE_PATH
RUN mkdir oscode
RUN git clone --single-branch --recursive --branch master https://github.com/fruzsinaagocs/oscode oscode/

RUN cd oscode && \
    python3 setup.py install

# create user with a home directory
ARG NB_USER
ARG NB_UID
ENV USER ${NB_USER}
ENV HOME /home/${NB_USER}

RUN adduser --disabled-password \
    --gecos "Default user" \
    --uid ${NB_UID} \
    ${NB_USER}
WORKDIR ${HOME}
USER ${USER}

# Make sure the contents of our repo are in ${HOME}
COPY . ${HOME}
USER root
RUN chown -R ${NB_UID} ${HOME}
USER ${NB_USER}
