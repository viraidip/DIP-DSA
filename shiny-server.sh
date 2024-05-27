export IMAGE="dip-dsa_r-shiny_app:latest"
export FILE="Dockerfile"

# stop old container
docker rm $(docker stop $(docker ps -a --filter ancestor="dip-dsa_r-shiny_app" --format="{{.ID}}"))
docker rmi -f $IMAGE

DOCKER_BUILDKIT=1 docker build --no-cache -f $FILE -t $IMAGE .
docker run -p 8161:8161 -d --name "dip-dsa_r-shiny_app" --restart always $IMAGE