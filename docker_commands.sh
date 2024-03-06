docker rmi dip-dsa:latest
docker build . --no-cache -t dip-dsa

docker run -p 8161:8161 dip-dsa