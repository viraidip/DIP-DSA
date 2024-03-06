FROM rocker/shiny:latest

COPY ./app/ /DIP-DSA/

RUN apt-get update && apt-get install libcurl4-openssl-dev libssl-dev -y
RUN install2.r shiny shinydashboard ggplot2 DT stringr dplyr \
        hash plotly plyr jsonlite shinyvalidate
RUN R -e "install.packages(c('BiocManager'), dependencies=TRUE, repos='https://cloud.r-project.org')"
RUN R -e "BiocManager::install('Biostrings')"
RUN R -e "BiocManager::install('ComplexHeatmap')"

RUN chmod -R 777 /DIP-DSA/
WORKDIR /DIP-DSA/

CMD ["R", "--no-save", "-e", "shiny::runApp('app.R', port=8161, host='0.0.0.0')"]

EXPOSE 8161