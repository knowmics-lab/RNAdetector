FROM ubuntu:focal

ENV DOCKER_USER_ID 501
ENV DOCKER_USER_GID 20

ENV BOOT2DOCKER_ID 1000
ENV BOOT2DOCKER_GID 50

ENV SUPERVISOR_VERSION 4.2.0

ENV DEBUG false

# Tweaks to give Apache/PHP write permissions to the app
RUN usermod -u ${BOOT2DOCKER_ID} www-data && \
    usermod -G staff www-data && \
    useradd -r mysql && \
    usermod -G staff mysql && \
    groupmod -g $(($BOOT2DOCKER_GID + 10000)) $(getent group $BOOT2DOCKER_GID | cut -d: -f1) && \
    groupmod -g ${BOOT2DOCKER_GID} staff

ENV DEBIAN_FRONTEND noninteractive

# Setup repositories
RUN /usr/bin/apt update --fix-missing && \
    /usr/bin/apt install -y perl dialog software-properties-common &&  \
    apt-key adv --keyserver keyserver.ubuntu.com --recv-keys 51716619E084DAB9 && \
    apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9 && \
    add-apt-repository 'deb https://cloud.r-project.org/bin/linux/ubuntu focal-cran40/' &&  \
    add-apt-repository ppa:ondrej/php -y && \
    add-apt-repository universe && \
    add-apt-repository multiverse && \
    /usr/bin/apt update && \
    /usr/bin/apt dist-upgrade -y && \
    apt-get clean

# Install common packages
RUN /usr/bin/apt install -y wget unzip curl tar openjdk-14-jre grep tzdata rename apt-utils nano beanstalkd apache2 \
                            jq pwgen mysql-server zip pigz git python-setuptools python3-setuptools pandoc python-dev \
                            libssl-dev r-base-dev cmake wget cython swig python3-pip libxml2-dev libcurl4-gnutls-dev \
                            libyaml-dev libgit2-dev php7.4-xdebug composer php7.4-cli libapache2-mod-php7.4 \
                            php7.4-mysql php7.4-apcu php7.4-gd php7.4-xml php7.4-mbstring php7.4-zip php7.4-curl \
                            php7.4-intl php7.4-json cufflinks tabix bwa samtools bedtools subread cutadapt fastqc \
                            rna-star hisat2 rna-star && \
    apt autoremove -y && \
    echo "ServerName localhost" >> /etc/apache2/apache2.conf

# Fix certificate issues
RUN apt-get update && \
    apt-get install ca-certificates-java && \
    apt-get clean && \
    update-ca-certificates -f

# Install supervisor 4
RUN curl -L https://pypi.io/packages/source/s/supervisor/supervisor-${SUPERVISOR_VERSION}.tar.gz | tar xvz && \
    cd supervisor-${SUPERVISOR_VERSION}/ && \
    python3 setup.py install && \
    cd .. && rm -r supervisor-${SUPERVISOR_VERSION}/

# Enable mod_rewrite in apache
RUN a2enmod rewrite

# Remove pre-installed database
RUN rm -rf /var/log/mysql && rm -rf /var/lib/mysql

# Create RNAdetector and temp dir
RUN mkdir -p /rnadetector/tmp/

# Install trim_galore
RUN PW=$(pwd) && cd /rnadetector/tmp/ && \
    curl -fsSL https://github.com/FelixKrueger/TrimGalore/archive/0.6.5.tar.gz -o trim_galore.tar.gz && \
    tar -zxvf trim_galore.tar.gz && \
    cp TrimGalore-0.6.5/trim_galore /usr/local/bin/ && cd "$PW"

# Install pandoc 2.9.2.1
RUN PW=$(pwd) && cd /rnadetector/tmp/ && \
    apt remove -y pandoc && \
    curl -fsSL https://github.com/jgm/pandoc/releases/download/2.9.2.1/pandoc-2.9.2.1-1-amd64.deb -o pandoc.deb && \
    dpkg -i pandoc.deb && cd "$PW"

# Install salmon
RUN PW=$(pwd) && cd /rnadetector/tmp/ && \
    curl -fsSL https://github.com/COMBINE-lab/salmon/releases/download/v1.4.0/salmon-1.4.0_linux_x86_64.tar.gz -o salmon.tar.gz && \
    tar -zxvf salmon.tar.gz && \
    mv salmon-latest_linux_x86_64/ /opt/salmon/ && \
    ln -s /opt/salmon/bin/salmon /usr/bin/salmon && cd "$PW"

# Install bbmap for re-pair utility
RUN PW=$(pwd) && cd /rnadetector/tmp/ && \
    curl -fsSL https://sourceforge.net/projects/bbmap/files/latest/download -o bbmap.tar.gz && \
    tar -zxvf bbmap.tar.gz --directory=/opt/ && \
    chmod 755 /opt/bbmap/* && \
    rm bbmap.tar.gz && cd "$PW"

# Install stringtie from source
RUN PW=$(pwd) && cd /rnadetector/tmp/ && \
    curl -fsSL http://ccb.jhu.edu/software/stringtie/dl/stringtie-2.0.6.tar.gz -o stringtie.tar.gz && \
    tar -xzvf stringtie.tar.gz && \
    cd /rnadetector/tmp/stringtie-2.0.6 && \
    make release && \
    cp stringtie /usr/local/bin/stringtie && \
    cp prepDE.py /usr/local/bin/prepDE.py && \
    cd /rnadetector/tmp && \
    rm -rf stringtie-2.0.5.Linux_x86_64/  && \
    rm stringtie.tar.gz  && cd "$PW"

# Install latest version of CIRIquant
RUN PW=$(pwd) && cd /rnadetector/tmp/ && \
    curl -fsSL https://downloads.sourceforge.net/project/ciri/CIRIquant/CIRIquant_v1.0.tar.gz -o CIRIquant.tar.gz && \
    tar -zxvf CIRIquant.tar.gz && cd /rnadetector/tmp/CIRIquant && \
    python setup.py install && cd "$PW"

# Install latest version of fastq-pair
RUN PW=$(pwd) && cd /rnadetector/tmp/ && \
    git clone https://github.com/linsalrob/fastq-pair.git && \
    cd fastq-pair/ && mkdir build && cd build && cmake .. && make && make install && cd "$PW"

# Install latest version of htseq-count
RUN PW=$(pwd) && cd /rnadetector/tmp/ && \
    curl -fsSL https://github.com/simon-anders/htseq/archive/release_0.11.1.tar.gz -o htseq.tar.gz && \
    tar -xzvf htseq.tar.gz && \
    cd /rnadetector/tmp/htseq-release_0.11.1 && \
    python setup.py build && \
    python setup.py install && \
    cd /rnadetector/tmp && \
    rm -rf htseq-release_0.11.1/ && \
    rm htseq.tar.gz && cd "$PW"

# Install deeptools
RUN pip3 install deeptools

# Remove tempo directory
RUN rm -r /rnadetector/tmp/

# Add all support files
ADD  scripts/bootstrap.sh                            /usr/local/bin/
ADD  scripts/installPackages.R                       /usr/local/bin/
ADD  scripts/setup.sh                                /usr/local/bin/
ADD  scripts/create_mysql_users.sh                   /usr/local/bin/
ADD  scripts/create_mysql_users_fast.sh              /usr/local/bin/
ADD  scripts/CIRI1.pl                                /usr/local/bin/CIRI1.pl
ADD  scripts/CIRI2.pl                                /usr/local/bin/CIRI2.pl
ADD  scripts/genkey.sh                               /genkey.sh
ADD  scripts/import_reference.sh                     /import_reference.sh
ADD  scripts/update_check.sh                         /update_check.sh
ADD  scripts/update_run.sh                           /update_run.sh
ADD  scripts/start-apache2.sh                        /usr/local/bin/start-apache2.sh
ADD  scripts/start-mysqld.sh                         /usr/local/bin/start-mysqld.sh
ADD  conf/apache_default.conf                        /etc/apache2/sites-available/000-default.conf
ADD  conf/supervisord-apache2.conf                   /etc/supervisor/conf.d/supervisord-apache2.conf
ADD  conf/supervisord-mysqld.conf                    /etc/supervisor/conf.d/supervisord-mysqld.conf
ADD  conf/supervisord-beanstalkd.conf                /etc/supervisor/conf.d/supervisord-beanstalkd.conf
ADD  conf/supervisord-worker.conf                    /etc/supervisor/conf.d/supervisord-worker.conf
ADD  conf/supervisord.conf                           /etc/supervisor/supervisord.conf
ADD  conf/mysqld_innodb.cnf                          /etc/mysql/conf.d/mysqld_innodb.cnf
COPY database.tar.bz2                                /opt/database.tar.bz2
COPY repo.tar.gz                                     /repo.tar.gz

RUN chmod 755 /usr/local/bin/*
RUN chmod 755 /genkey.sh
RUN chmod 755 /import_reference.sh

# Install R packages
RUN /usr/bin/Rscript /usr/local/bin/installPackages.R

# Complete the setup
RUN /bin/bash /usr/local/bin/setup.sh

# Link bootstrap script
RUN ln -s /usr/local/bin/bootstrap.sh /

# Clean up the container
RUN /usr/bin/apt autoclean && \
    /usr/bin/apt clean && \
    rm -rf /var/lib/apt/lists/*

EXPOSE 80
VOLUME [ "/rnadetector/ws/storage/app/" ]

CMD ["/bootstrap.sh"]
