# Pour eviter de creer des images multiples: --rm
MY_OS=rockylinux:8
#MY_OS=redhat/ubi8 # N'a pas openmpi...
#MY_OS=fedora:41
TRUST_ROOT=`pwd`/../../
exec=$TRUST_ROOT/exec/TRUST_mpi_opt
# Build the container to build TRUST on a specific OS
docker build --rm --build-arg MY_OS=$MY_OS -t my_$MY_OS .

docker run -v $TRUST_ROOT:/home/trust my_$MY_OS `[ -f $exec ] && echo "Binary $exec ready."`

# Build the container with the binary for the specific OS
image=trust_$MY_OS
cp $exec .
docker build --rm --build-arg MY_OS=$MY_OS -f Dockerfile2 -t $image .
rm -f ./TRUST_mpi_opt
touch toto.data
docker run $image toto . || exit -1

# Changement de / et : en _ pour le nom des fichiers:
name_image=`echo $image | awk '{gsub("\\\/","_",$0);gsub(":","_",$0);print $0}'`

# Container pushed on hub.docker.com:
tag=pledac/$name_image
docker tag $image $tag 	|| exit -1
docker push $tag	|| exit -1

# Docker image in a compressed format:
save=$TRUST_ROOT/$name_image.tar.gz
docker save $image | gzip > $save
echo "$save created."


