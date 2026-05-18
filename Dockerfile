FROM rocker/shiny:latest

RUN apt-get update && \
    apt-get upgrade -y && \
    apt-get install -y --no-install-recommends \
    dirmngr gnupg2 libcurl4-openssl-dev libssl-dev \
    xz-utils wget ca-certificates build-essential python3 \
    && rm -rf /var/lib/apt/lists/*

RUN wget https://nodejs.org/dist/v22.22.0/node-v22.22.0-linux-x64.tar.xz \
    && tar -xf node-v22.22.0-linux-x64.tar.xz \
    && cp -R node-v22.22.0-linux-x64/bin/* /usr/local/bin/ \
    && cp -R node-v22.22.0-linux-x64/lib/* /usr/local/lib/ \
    && ln -sf /usr/local/bin/node /opt/shiny-server/bin/node \
    && rm -rf node-v22.22.0-linux-x64*

WORKDIR /tmp/patch
COPY package.json /tmp/patch/package.json
RUN npm install -g nan && \
    npm install --omit=dev && \
    cp -R node_modules/* /opt/shiny-server/node_modules/ && \
    rm -rf /tmp/patch

RUN install2.r --error --skipinstalled \
    shiny shinydashboard ggplot2 DT stringr dplyr \
    hash plotly plyr jsonlite shinyvalidate BiocManager \
    && R -e "BiocManager::install(c('Biostrings', 'ComplexHeatmap'), ask=FALSE)" \
    && rm -rf /tmp/downloaded_packages

WORKDIR /srv/shiny-server
COPY ./app/ /srv/shiny-server/dipdsa
COPY ./shiny-server.conf /etc/shiny-server/shiny-server.conf

RUN chmod -R 777 /srv/shiny-server/ \
    && sed -i -e 's/\blisten 3838\b/listen 8161/g' /etc/shiny-server/shiny-server.conf

USER shiny

EXPOSE 8161
CMD ["/usr/bin/shiny-server"]
