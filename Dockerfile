FROM ubuntu:22.04
#ARG DEBIAN_FRONTEND=noninteractive
ENV DEBIAN_FRONTEND=noninteractive

# Update and install necessary dependencies
RUN apt-get -y update && apt-get upgrade -y
RUN apt-get -y install curl build-essential bison flex subversion git mercurial wget \
    libgsl0-dev vim flex bison cmake zlib1g-dev libboost-system-dev libboost-thread-dev \
    libopenmpi-dev openmpi-bin gnuplot sudo

# Install OpenFOAM
RUN curl -s https://dl.openfoam.com/add-debian-repo.sh | bash && \
    wget -q -O - https://dl.openfoam.com/add-debian-repo.sh | bash && \
    apt-get -y update && apt-get install -y openfoam2312-default

# Install additional numerical libraries
RUN apt-get -y install libblas-dev liblapack-dev libopenblas-dev liblapacke-dev 

# Set up user
WORKDIR "/home/"
RUN adduser user && usermod -aG sudo user && apt install sudo -y
RUN echo "user:0000" | chpasswd
USER user

# Ensure OpenFOAM environment loads correctly
RUN echo ". /usr/lib/openfoam/openfoam2312/etc/bashrc" >> ~/.bashrc
SHELL [ "/bin/bash", "-c" ]

# Clone the repository
WORKDIR "/home/user"
RUN git clone https://github.com/ManuelaRC/TestingOF2312.git
RUN echo "Cloned respository ManuelaRC/TestingOF2312"
WORKDIR "/home/user/TestingOF2312"

# Set environment variables for user OpenFOAM directory
ARG FOAM_USER_LIBBIN="/home/user/OpenFOAM/user-v2312/platforms/linux64GccDPInt32Opt/lib/"
ENV FOAM_USER_LIBBIN=$FOAM_USER_LIBBIN
ARG WM_PROJECT_USER_DIR="/home/user/OpenFOAM/user-v2312/"
ENV WM_PROJECT_USER_DIR=$WM_PROJECT_USER_DIR

USER root
# Grant user permissions to OpenFOAM root directory
RUN sudo chmod -R 777 /usr/lib/openfoam/openfoam2312/


#### SETTING NEW TURBULENCE MODEL
# Copy kEpsilonFp to the RAS turbulence models directory
RUN cp -r kEpsilonFp /usr/lib/openfoam/openfoam2312/src/TurbulenceModels/turbulenceModels/RAS/
# Modify turbulentTransportModels.C
RUN sed -i '/makeRASModel(kEpsilon);/a #include "kEpsilonFp.H"\nmakeRASModel(kEpsilonFp);' \
    /usr/lib/openfoam/openfoam2312/src/TurbulenceModels/incompressible/turbulentTransportModels/turbulentTransportModels.C
# Compile the turbulence models
WORKDIR "/usr/lib/openfoam/openfoam2312/src/TurbulenceModels"
RUN source /usr/lib/openfoam/openfoam2312/etc/bashrc && ./Allwclean && ./Allwmake
# Confirmation message
RUN echo "kEpsilonFp installed"


#### SETTING NEW FVOPTION
WORKDIR "/home/user/TestingOF2312"
RUN ls -la

# Copy calibratedActuatorDisk to the fvOption folder
RUN cp -r calibratedActuatorDisk /usr/lib/openfoam/openfoam2312/src/fvOptions/sources/derived/
# Modify turbulentTransportModels.C
RUN sed -i '\|$(derivedSources)/actuationDiskSource/actuationDiskSource.C|a $(derivedSources)/calibratedActuatorDisk/calibratedActuatorDisk.C' \
    /usr/lib/openfoam/openfoam2312/src/fvOptions/Make/files
RUN sed -i '\|$(derivedSources)/actuationDiskSource/actuationDiskSource.C|a $(derivedSources)/calibratedActuatorDisk/myLookupTable4Turbines/myLookupTable4Turbines.C' \
    /usr/lib/openfoam/openfoam2312/src/fvOptions/Make/files
WORKDIR "/usr/lib/openfoam/openfoam2312/src/fvOptions"
RUN source /usr/lib/openfoam/openfoam2312/etc/bashrc && wclean .
RUN source /usr/lib/openfoam/openfoam2312/etc/bashrc && wmake
# Confirmation message
RUN echo "fvOption::calibratedActuatorDisk installed"


#### COPYING TUTORIALS
WORKDIR "/home/user/TestingOF2312"
RUN cp -r tutorials /usr/lib/openfoam/openfoam2312/tutorials/incompressible/simpleFoam/
# Confirmation message
RUN echo "tutorials of calibratedActuatorDisk copied to tutorials/incompressible/simpleFoam/"

USER user

# Set working directory back to user's home
WORKDIR "/home/user"

