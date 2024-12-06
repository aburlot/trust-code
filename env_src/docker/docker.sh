# Pour eviter de creer des images multiples: --rm
MY_OS=rockylinux:8
#MY_OS=redhat/ubi8 # N'a pas openmpi...
#MY_OS=fedora:41

# Nom de l'image
image=pledac/`echo trust-$MY_OS | awk '{gsub("\\\/","-",$0);gsub(":","-",$0);print $0}'`
TRUST_ROOT=`pwd`/../..
exec=$TRUST_ROOT/exec/TRUST_mpi_opt

# Build the container to build TRUST on a specific OS
docker build --rm --build-arg MY_OS=$MY_OS -t $image .

# Run the container 
docker run -v $TRUST_ROOT:/home/trust $image `[ -f $exec ] && echo "Binary $exec ready."`

# Build the container with the binary for the specific OS
cp $exec .
docker build --rm --build-arg MY_OS=$MY_OS -f Dockerfile2 -t $image .
rm -f ./TRUST_mpi_opt

# Test the container locally with a data file:
echo "system \"echo Docker works\"" > dumb.data
docker run -v `pwd`:/tmp $image dumb . || exit -1
rm -f dumb.data


# Container pushed on hub.docker.com: Le tag pourrait etre plus tot...
docker push $image	|| exit -1

# Docker image in a compressed format for CCRT to use by pcocc:
save=$TRUST_ROOT/`basename $image.tar.gz`
docker save $image | gzip > $save
echo "$save created."


