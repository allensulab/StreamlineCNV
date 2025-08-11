FROM nfcore/tools:latest
LABEL description="StreamlineCNV"

RUN apt-get update -y && apt-get upgrade -y && apt-get install -y build-essential cmake r-base libhts3 libhts-dev fastqc samtools bcftools gnuplot autoconf automake libtool libgtextutils0v5 libgtextutils-dev libgd-perl libperlio-gzip-perl bwa python3-seaborn python3-sklearn python3-sklearn-lib python-sklearn-doc

RUN mkdir -p /tmp && cd /tmp &&\
    git clone https://github.com/shahcompbio/hmmcopy_utils &&\ 
    cd hmmcopy_utils &&\
    cmake . &&\
    make &&\
    cp bin/* /usr/local/bin

RUN mkdir -p /tmp && cd /tmp &&\
    git clone https://github.com/Debian/fastx-toolkit &&\
    cd fastx-toolkit &&\
    sed -i '105s/usage();/usage();break;/' src/fasta_formatter/fasta_formatter.cpp &&\
    ./reconf &&\
    ./configure &&\
    make &&\
    make install

RUN R -e 'if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager", lib="/usr/lib/R/library/")'
RUN R -e 'BiocManager::install("HMMcopy", lib="/usr/lib/R/library/")'
RUN R -e 'BiocManager::install("IRanges", lib="/usr/lib/R/library/")'

RUN mkdir -p /opt/scripts/
COPY ./scripts/ /opt/scripts/
