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
FROM fedora:latest
RUN dnf -y install patch g++ gfortran
RUN dnf -y install libX11-devel mesa-libGLU-devel texlive-scheme-full
RUN dnf -y install file bc diffutils util-linux tcl tk hostname libX11-devel unzip
RUN dnf -y install openmpi openmpi-devel
ENV USER=root
# Build during run:
ENTRYPOINT [ "/bin/sh" , "-c", "cd /tmp;./configure -force_system_mpi -without-pdflatex && make" ]

