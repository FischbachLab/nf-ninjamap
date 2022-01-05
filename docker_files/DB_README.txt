# Docker
docker container run -it --rm -v ${PWD}:/data biocontainers/grinder:v0.5.4-1-deb_cv1

# Grinder Command
grinder \
    -cf 100 \
    -rd 151 \
    -id 800 \
    -mo FR \
    -dc '-~*NX' \
    -md poly4 3e-3 3.3e-8 \
    -rs 1712 \
    -am exponential \
    -ql 30 10 \
    -fq 1 \
    -od output_dir \
    -bn b_vulgatus \
    -rf