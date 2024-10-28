# python 3.7.12 full image
FROM python:3.11.8

# Set the working directory in the container
WORKDIR /app

# Install system dependencies
RUN apt-get update && apt-get install -y \
    git \
    gcc \
    g++ \
    make \
    libpython3-dev \
    && rm -rf /var/lib/apt/lists/*

# pip install package dependencies
RUN pip install scikit-learn==1.3.0 biopython==1.83 pandas

# pip install autopeptideml
RUN pip install autopeptideml

# cleanup
RUN apt-get remove -y git gcc g++ make libpython3-dev && apt-get autoremove -y

# environment path
ENV PATH="/app/autopeptideml/scripts:${PATH}"
