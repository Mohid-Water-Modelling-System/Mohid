# MOHID Docker Repository

This directory contains Dockerfiles and build instructions focused on compiling and running MOHID tools.


## List of public Docker Images
Public Docker images are pre-built versions of MOHID tools and models, hosted on a Docker registry. These images can be pulled directly and are ready to use without requiring local compilation or additional dependencies.

 ### Models
- **MOHID Water**
 > docker pull SOON

- **MOHID Land**
 > docker pull SOON
 
 ### Tools
 - **Compare2HDFfiles**
 > docker pull SOON
- **Convert2HDF5**
 > docker pull SOON
- **Convert2NC**
 > docker pull SOON


## How to generate your own Images
If you prefer to customize the Docker images or if the public images are not yet available, you can generate your own images using the provided Dockerfiles. This process involves downloading the source code and dependencies, and then building the images locally. Generating your own images may also be useful for testing or modifying the tools to suit specific needs.

### Models
- **MOHID Water**
 > docker build -f "Solutions\dockers\MohidWater\dockerfile" -t mohidwater:latest "."

- **MOHID Land**
 > docker build -f "Solutions\dockers\MohidLand\dockerfile" -t mohidland:latest "."
 
 
 ### Tools
 - **Compare2HDFfiles**
 > docker build -f "Solutions\dockers\Compare2HDFfiles\dockerfile" -t compare2hdffiles:latest "."
- **Convert2HDF5**
 > docker build -f "Solutions\dockers\Convert2HDF5\dockerfile" -t convert2hdf5:latest "."
- **Convert2NC**
 > docker build -f "Solutions\dockers\Convert2NC\dockerfile" -t convert2nc:latest "."

### Key Differences Between Hosted and Locally Built Images

- **Hosted Images**:
  - Pre-built and optimized for general use.
  - Can be pulled directly from a Docker registry (e.g., Docker Hub).
  - Simplifies deployment without requiring build tools or source code.

- **Locally Built Images**:
  - Requires downloading source code and dependencies.
  - Enables customization of tools or configurations before building.
  - Useful for testing new features or working on development versions.

By understanding these differences, you can choose the approach that best fits your workflow and requirements.

We will try to update public hosted images yearly.