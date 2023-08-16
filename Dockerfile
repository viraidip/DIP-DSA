FROM rocker/shiny:latest

COPY ./app/ /srv/DIP-DSA/
COPY ./data/ /srv/DIP-DSA/data/

RUN apt-get update && apt-get install libcurl4-openssl-dev libssl-dev -y
RUN install2.r shiny shinydashboard ggplot2 DT comprehenr hash plotly plyr jsonlite tools reticulate shinyvalidate ggvenn
RUN R -e "install.packages(c('BiocManager'), dependencies=TRUE, repos='https://cloud.r-project.org')"
RUN R -e "BiocManager::install('Biostrings')"

RUN chmod -R 777 /srv/DIP-DSA/
EXPOSE 3838
WORKDIR /srv/DIP-DSA

CMD ["R", "-e", "shiny::runApp(port=3838, host=127.0.0.1"]
