# docker ps                 : list containers
# docker images             : list the images
# docker build ...          : build a image
# docker run -it ...        : run an image into a container
# docker volume create

# See https://docs.docker.com/reference/dockerfile/
# Dockerfile commands to build an image:
# FROM                      : insert an image
# COPY                      : copy files into the image
# ADD                       : more complex than COPY (URL?)
# RUN                       : run a command during image build
# CMD/ENTRYPOINT            : run a command when the container is launched (only 1 CMD/ENTRYPOINT per Dockerfile !!!)
ARG MY_OS
FROM ${MY_OS}
RUN if [ "${MY_OS}" = "fedora:41" ]; then \
dnf -y install gcc g++ gfortran; \
dnf -y install perl make patch file bc diffutils util-linux tcl tk hostname libX11-devel unzip; \
dnf -y install openmpi openmpi-devel; \
dnf -y install mesa-libGLU-devel texlive-scheme-full; \
else \
dnf -y install gcc gcc-c++ gcc-gfortran; \
dnf -y install python3 perl make patch file bc diffutils util-linux tcl tk hostname libX11-devel unzip; \
fi 
ENV USER=root
# Build during run:
#ENTRYPOINT [ "/bin/bash" , "-c", "cd /tmp;./configure -force_provided_openmpi -without-pdflatex && make" ]

