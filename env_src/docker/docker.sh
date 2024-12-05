# Pour eviter de creer des images multiples: --rm
MY_OS=redhat/ubi8
#MY_OS=fedora:41
docker build --rm --build-arg MY_OS=$MY_OS -t my_$MY_OS . && docker run -v /home/formation/s-sac-dm2s-train9/fedora/trust-code:/tmp my_$MY_OS
