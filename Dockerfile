FROM rocker/shiny:latest

RUN apt-get update && \
    apt-get upgrade -y && \
    apt-get purge -y dirmngr && \
    apt-get install -y dirmngr=2.4.4-2ubuntu17.4 gnupg2 && \
    apt-get install -y libcurl4-openssl-dev libssl-dev xz-utils wget npm && \
    rm -rf /var/lib/apt/lists/*

RUN wget https://nodejs.org/dist/v25.8.2/node-v25.8.2-linux-x64.tar.xz \
    && tar -xf node-v25.8.2-linux-x64.tar.xz \
    && cp node-v25.8.2-linux-x64/bin/node /opt/shiny-server/bin/node \
    && cp node-v25.8.2-linux-x64/bin/node /opt/shiny-server/ext/node/bin/node \
    && rm -rf node-v25.8.2-linux-x64*

WORKDIR /opt/shiny-server
RUN npm install \
    handlebars@4.7.9 \
    flatted@3.4.2 \
    underscore@1.13.8 \
    tar@7.5.11 \
    qs@6.14.2 \
    minimatch@10.2.1 \
    glob@11.1.0 \
    cross-spawn@7.0.5 \
    brace-expansion@5.0.5 \
    tar@7.5.11 \
    qs@6.14.2 \
    minimatch@10.2.3


RUN install2.r --error --skipinstalled \
    shiny shinydashboard ggplot2 DT stringr dplyr \
    hash plotly plyr jsonlite shinyvalidate BiocManager \
    && R -e "BiocManager::install(c('Biostrings', 'ComplexHeatmap'), ask=FALSE)" \
    && rm -rf /tmp/downloaded_packages

COPY ./app/ /srv/shiny-server/dipdsa
COPY ./shiny-server.conf /etc/shiny-server/shiny-server.conf

RUN chmod -R 777 /srv/shiny-server/ \
    && sed -i -e 's/\blisten 3838\b/listen 8161/g' /etc/shiny-server/shiny-server.conf

USER shiny

EXPOSE 8161
CMD ["/usr/bin/shiny-server"]
