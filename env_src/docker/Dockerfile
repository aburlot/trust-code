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
# RUN                       : command during build
# CMD/ENTRYPOINT 	    : command during run with parameters (CMD can be overidden in the CLI)
ARG MY_OS
FROM ${MY_OS}
RUN useradd -ms /bin/bash trust 
WORKDIR /home/trust
RUN if [ "${MY_OS}" = "fedora:41" ]; then \
dnf -y install gcc g++ gfortran; \
dnf -y install wget python3 perl make patch file bc diffutils util-linux tcl tk hostname libX11-devel unzip; \
dnf -y install openmpi openmpi-devel; \
dnf -y install mesa-libGLU-devel texlive-scheme-full; \
else \
dnf -y install gcc gcc-c++ gcc-gfortran; \
dnf -y install wget python3 perl make patch file bc diffutils util-linux tcl tk hostname libX11-devel unzip; \
dnf -y install openmpi openmpi-devel; \
fi 
#ENV USER=root
# Build during run:
CMD [ "/bin/bash" , "-c", "cd /home/trust;./configure -force_system_mpi -without-pdflatex && make optim" ]
#ENTRYPOINT [ "/bin/bash" , "-c", "echo Container OK" ]
