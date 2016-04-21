FROM centos:centos6

MAINTAINER Artem Tarasov <artem.tarasov@embl.de>
WORKDIR /root
RUN rm /bin/sh && ln -s /bin/bash /bin/sh

RUN yum -y update && yum -y install zlib-static libxml2-static libxml2-devel make cmake git wget autoconf automake unzip && yum clean all

# install g++ 4.9
RUN yum -y localinstall https://www.softwarecollections.org/en/scls/rhscl/devtoolset-3/epel-6-x86_64/download/rhscl-devtoolset-3-epel-6-x86_64.noarch.rpm
RUN yum -y install devtoolset-3-gcc devtoolset-3-binutils devtoolset-3-gcc-c++ && yum clean all
RUN /usr/bin/scl enable devtoolset-3 true

ENV CC=/opt/rh/devtoolset-3/root/usr/bin/gcc
ENV CXX=/opt/rh/devtoolset-3/root/usr/bin/g++

# let's link to liblzma statically
RUN curl -O http://tukaani.org/xz/xz-5.2.2.tar.gz &&\
    tar xzf xz-5.2.2.tar.gz &&\
    cd xz-5.2.2 &&\
    ./configure CFLAGS='-fPIC -O2 -mtune=generic' && make && make install

# install conda
RUN wget -q http://repo.continuum.io/miniconda/Miniconda-latest-Linux-x86_64.sh &&\
    bash Miniconda-latest-Linux-x86_64.sh -b -p /miniconda &&\
    rm Miniconda-latest-Linux-x86_64.sh
ENV PATH /miniconda/bin:$PATH

RUN conda create -y -n py3 python=3; source activate py3; pip install auditwheel pypatchelf
RUN conda create -y -n py2 python=2 cffi