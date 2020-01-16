FROM rocker/shiny-verse:3.6.1

MAINTAINER juozas@vugene.eu


# System requirements
RUN apt-get update && apt-get install -y \
  build-essential \
  curl \
  parallel \
  python3-pip \
  libmagick++-dev \
  && rm -rf /var/lib/apt/lists/*

# R packages
RUN install2.r --error --dep TRUE \
  kableExtra \
  DT

# Configuration
COPY shiny-server.conf /etc/shiny-server/shiny-server.conf

# Start the server
EXPOSE 3838
CMD ["/usr/bin/shiny-server.sh"]
