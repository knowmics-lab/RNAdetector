FROM ubuntu:latest

RUN DEBIAN_FRONTEND=noninteractive /usr/bin/apt update && \
    DEBIAN_FRONTEND=noninteractive /usr/bin/apt install -y perl dialog software-properties-common &&  \
    apt-key adv --keyserver keyserver.ubuntu.com --recv-keys 51716619E084DAB9 && \
    add-apt-repository 'deb https://cloud.r-project.org/bin/linux/ubuntu bionic-cran35/' &&  \
    DEBIAN_FRONTEND=noninteractive /usr/bin/apt update && \
    DEBIAN_FRONTEND=noninteractive /usr/bin/apt install -y init && \
    DEBIAN_FRONTEND=noninteractive /usr/bin/apt dist-upgrade -y && \
    DEBIAN_FRONTEND=noninteractive apt-get clean

RUN DEBIAN_FRONTEND=noninteractive /usr/bin/apt install -y htop nginx php-fpm php-json php-sqlite3 php-intl php-mbstring php-xml
RUN DEBIAN_FRONTEND=noninteractive /usr/bin/apt install -y supervisor git composer python-setuptools python-dev curl
RUN DEBIAN_FRONTEND=noninteractive /usr/bin/apt install -y tophat bwa samtools python-htseq subread cutadapt fastqc salmon
RUN DEBIAN_FRONTEND=noninteractive /usr/bin/apt install -y r-base-dev libxml2-dev libcurl4-openssl-dev
RUN systemctl enable nginx && systemctl enable php7.2-fpm && systemctl enable supervisor

COPY bootstrap.sh /usr/local/bin/
COPY installPackages.R /usr/local/bin/
COPY setup.sh /usr/local/bin/
COPY repo.tar.gz /repo.tar.gz

#RUN /usr/bin/Rscript /usr/local/bin/installPackages.R

RUN /bin/bash /usr/local/bin/setup.sh

RUN ln -s /usr/local/bin/bootstrap.sh /

EXPOSE 80

ENTRYPOINT ["/bootstrap.sh"]
CMD ["/sbin/init"]
