export IMAGE="dip-dsa_r-shiny_app:latest"
export FILE="Dockerfile"
DOCKER_BUILDKIT=1 docker build --no-cache -f $FILE -t $IMAGE .
docker run -p 8161:8161 $IMAGE