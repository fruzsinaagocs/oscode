FROM python:3.7-slim
# install the notebook package
RUN pip install --no-cache --upgrade pip && \
    pip install --no-cache notebook

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


RUN apt-get update && apt-get install -y \
    git \
RUN git clone https://github.com/eigenteam/eigen-git-mirror 
RUN export CPLUS_INCLUDE_PATH="$PWD/eigen-git-mirror/:$CPLUS_INCLUDE_PATH"
RUN pip install -r requirements.txt
RUN pip install .

