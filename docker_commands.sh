docker rmi dip-dsa_r-shiny_app:latest
docker build . --no-cache -t dip-dsa_r-shiny_app

docker run -p 8161:8161 dip-dsa_r-shiny_app