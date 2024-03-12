FROM rocker/shiny:latest

RUN apt-get update && apt-get install libcurl4-openssl-dev libssl-dev -y
RUN install2.r shiny shinydashboard ggplot2 DT stringr dplyr \
        hash plotly plyr jsonlite shinyvalidate
RUN R -e "install.packages(c('BiocManager'), dependencies=TRUE, repos='https://cloud.r-project.org')"
RUN R -e "BiocManager::install('Biostrings')"
RUN R -e "BiocManager::install('ComplexHeatmap')"

COPY ./app/ /srv/shiny-server/dipdsa
COPY ./shiny-server.conf /etc/shiny-server/shiny-server.conf
RUN chmod -R 777 /srv/shiny-server/

RUN sed -i -e 's/\blisten 3838\b/listen 8161/g' /etc/shiny-server/shiny-server.conf

USER shiny

EXPOSE 8161
CMD ["/usr/bin/shiny-server"]
