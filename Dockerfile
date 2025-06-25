FROM rocker/shiny:latest

RUN apt-get update
RUN apt-get dist-upgrade -y && \
    apt-get install libcurl4-openssl-dev libssl-dev xz-utils -y

RUN wget https://nodejs.org/dist/v20.19.3/node-v20.19.3-linux-x64.tar.xz \
    && tar -xf node-v20.19.3-linux-x64.tar.xz \
    && cp node-v20.19.3-linux-x64/bin/node /opt/shiny-server/bin/node \
    && cp node-v20.19.3-linux-x64/bin/node /opt/shiny-server/ext/node/bin/node \
    && rm -rf node-v20.19.3-linux-x64*

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