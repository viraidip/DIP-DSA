FROM rocker/shiny:latest

RUN apt-get update && \
    apt-get upgrade -y && \
    apt-get install -y --no-install-recommends \
    dirmngr gnupg2 libcurl4-openssl-dev libssl-dev \
    xz-utils wget ca-certificates build-essential python3 curl && \
    mkdir -p /etc/apt/keyrings && \
    curl -fsSL https://deb.nodesource.com/gpgkey/nodesource-repo.gpg.key | gpg --dearmor -o /etc/apt/keyrings/nodesource.gpg && \
    echo "deb [signed-by=/etc/apt/keyrings/nodesource.gpg] https://deb.nodesource.com/node_22.x nodistro main" | tee /etc/apt/sources.list.min.d/nodesource.list || \
    echo "deb [signed-by=/etc/apt/keyrings/nodesource.gpg] https://deb.nodesource.com/node_22.x nodistro main" | tee /etc/apt/sources.list.d/nodesource.list && \
    apt-get update && \
    apt-get install -y nodejs && \
    rm -rf /var/lib/apt/lists/*

RUN ln -sf /usr/bin/node /opt/shiny-server/bin/node

WORKDIR /tmp/patch
COPY package.json /tmp/patch/package.json
RUN npm install -g nan && \
    npm install --omit=dev --no-audit --no-fund && \
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

RUN sed -i -e 's/\blisten 3838\b/listen 8161/g' /etc/shiny-server/shiny-server.conf && \
    chown -R shiny:shiny /srv/shiny-server/

RUN mkdir -p /srv/shiny-server/dipdsa/tmp && \
    chown -R shiny:shiny /srv/shiny-server/dipdsa/tmp

USER shiny

EXPOSE 8161
CMD ["/usr/bin/shiny-server"]
