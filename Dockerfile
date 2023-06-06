FROM debian:buster-slim
SHELL ["/bin/bash", "-c"]

RUN mkdir -p /usr/share/man/man1 && \
    apt-get -qq update && \
    apt-get -qq -y install --no-install-recommends \
        build-essential \
        gnupg \
        libfftw3-dev \
        default-jdk \
        curl \
        python3 \
        python3-dev \
        python3-pip \
        zlib1g-dev \
        git \
        wget \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

RUN echo "deb [signed-by=/usr/share/keyrings/cloud.google.gpg] http://packages.cloud.google.com/apt cloud-sdk main" | tee -a /etc/apt/sources.list.d/google-cloud-sdk.list && \
    curl https://packages.cloud.google.com/apt/doc/apt-key.gpg | apt-key --keyring /usr/share/keyrings/cloud.google.gpg  add - && \
    apt-get update -y && apt-get install google-cloud-cli -y
      
RUN ln -s /usr/bin/python3 /usr/bin/python

RUN python -m pip install --upgrade pip --no-cache-dir && \
    python -m pip install scanpy pandas numpy scipy dill pyranges scikit-learn --no-cache-dir