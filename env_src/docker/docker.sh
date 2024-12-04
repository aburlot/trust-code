# Pour eviter de creer des images multiples: --rm
docker build --rm -t my_fedora . && docker run -v /home/formation/s-sac-dm2s-train9/fedora/trust-code:/tmp my_fedora
