set -ex

cp  ../VERSION ./VERSION
cp ../environment.yml ./environment.yml

function cleanup {
    rm  ./VERSION
    rm  ./environment.yml
}

trap cleanup EXIT


# bop it
docker build -t open2c/distiller_env:latest . 
docker run -it open2c/distiller_env:latest apt list | sed 's/\x1b\[[0-9;]*m//g' > ./apt.list
docker run -it open2c/distiller_env:latest conda list > ./conda.list
docker images
