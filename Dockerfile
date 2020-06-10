FROM ubuntu

RUN apt-get update && apt-get install -y \
    git python3-pip
RUN pip install --no-cache --upgrade pip && \
    pip install --no-cache notebook && \
    pip install scipy numpy matplotlib

RUN git clone https://github.com/eigenteam/eigen-git-mirror 
RUN export CPLUS_INCLUDE_PATH="$PWD/eigen-git-mirror/:$CPLUS_INCLUDE_PATH"
RUN mkdir oscode
RUN git clone https://github.com/fruzsinaagocs/oscode ${HOME}/oscode/

RUN cd oscode && \
    pip install .

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


