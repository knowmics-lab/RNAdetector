FROM ubuntu:latest

RUN DEBIAN_FRONTEND=noninteractive /usr/bin/apt update && \
    DEBIAN_FRONTEND=noninteractive /usr/bin/apt install -y perl dialog software-properties-common &&  \
    apt-key adv --keyserver keyserver.ubuntu.com --recv-keys 51716619E084DAB9 && \
    add-apt-repository 'deb https://cloud.r-project.org/bin/linux/ubuntu bionic-cran35/' &&  \
    add-apt-repository ppa:ondrej/php -y && \
    DEBIAN_FRONTEND=noninteractive /usr/bin/apt update && \
    DEBIAN_FRONTEND=noninteractive /usr/bin/apt install -y init && \
    DEBIAN_FRONTEND=noninteractive /usr/bin/apt dist-upgrade -y && \
    DEBIAN_FRONTEND=noninteractive apt-get clean

RUN DEBIAN_FRONTEND=noninteractive /usr/bin/apt install -y nginx php7.3-cli php7.3-fpm php7.3-json php7.3-sqlite3 \ 
                                                           php7.3-intl php7.3-mbstring php7.3-xml supervisor \ 
                                                           git composer python-setuptools python3-setuptools \ 
                                                           python-dev curl tophat bwa samtools python-htseq \ 
                                                           subread cutadapt fastqc salmon r-base-dev \ 
                                                           libxml2-dev libcurl4-openssl-dev
RUN systemctl enable nginx && systemctl enable php7.3-fpm && systemctl enable supervisor

COPY bootstrap.sh        /usr/local/bin/
COPY installPackages.R   /usr/local/bin/
COPY setup.sh            /usr/local/bin/
COPY CIRI.pl             /usr/local/bin/CIRI.pl
COPY repo.tar.gz         /repo.tar.gz
COPY nginx.conf          /nginx.conf
COPY worker.conf         /worker.conf
COPY genkey.sh           /genkey.sh
COPY import_reference.sh /import_reference.sh

RUN /usr/bin/Rscript /usr/local/bin/installPackages.R

RUN /bin/bash /usr/local/bin/setup.sh

RUN ln -s /usr/local/bin/bootstrap.sh /

RUN DEBIAN_FRONTEND=noninteractive /usr/bin/apt autoclean && \ 
    DEBIAN_FRONTEND=noninteractive /usr/bin/apt clean


EXPOSE 80
VOLUME [ "/rnadetector/ws/storage/app/" ]

ENTRYPOINT ["/bootstrap.sh"]
CMD ["/sbin/init"]
