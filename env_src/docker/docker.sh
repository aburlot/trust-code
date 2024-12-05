# Pour eviter de creer des images multiples: --rm
MY_OS=redhat/ubi8
#MY_OS=fedora:41
TRUST_ROOT=`pwd`/../../
# Build the container to build TRUST on a specific OS
docker build --rm --build-arg MY_OS=$MY_OS -t my_$MY_OS .
docker run -v $TRUST_ROOT:/home/trust my_$MY_OS

# Build the container with the binary for the specific OS
cp ../../exec/TRUST_mpi_opt .
docker build --rm --build-arg MY_OS=$MY_OS -f Dockerfile2 -t trust_$MY_OS .
rm -f ./TRUST_mpi_opt
docker run trust_$MY_OS .

# Changement de / et : en _ pour le nom des fichiers:
MY_OS=`echo $MY_OS | awk '{gsub("\\\/","_",$0);gsub(":","_",$0);print $0}'`

# Container pushed on hub.docker.com:
tag="pledac/"trust_$MY_OS
docker tag trust_$MY_OS $tag
docker push $tag

# Docker image in a compressed format:
save=$TRUST_ROOT/trust_"$MY_OS".tar.gz
docker save trust_$MY_OS | gzip > $save
echo "$save created."


