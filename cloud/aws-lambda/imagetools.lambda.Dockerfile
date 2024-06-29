FROM fedora:38

RUN dnf install -y wget python3-pip unzip git bc blosc java-21-openjdk python3-lxml libtiff-tools && dnf clean all

RUN wget https://downloads.openmicroscopy.org/bio-formats/7.3.0/artifacts/bftools.zip && \
    unzip bftools.zip -d /usr/local/share/applications && \
    wget https://github.com/glencoesoftware/bioformats2raw/releases/download/v0.9.2/bioformats2raw-0.9.2.zip && \
    unzip bioformats2raw-0.9.2.zip -d /usr/local/share/applications && \
    rm -rf bioformats2raw.zip && \
    rm -rf bftools.zip

WORKDIR /usr/local/bin

RUN ln -s /usr/local/share/applications/bftools/tiffcomment tiffcomment && \
    ln -s /usr/local/share/applications/bftools/showinf showinf && \
    ln -s /usr/local/share/applications/bioformats2raw-0.9.2/bin/bioformats2raw  bioformats2raw

COPY ./imagetools /imagetools

WORKDIR /imagetools
RUN pip3 install . && pip3 install awslambdaric && pip3 cache purge

# Copy credentials - AWS provisioned containers run with the permissions of a default user named “sbx_user1051”
RUN useradd -ms /bin/bash sbx_user1051 
COPY ./imagetools/.deriva /home/sbx_user1051/.deriva
# Copy lambda handler
COPY ./imagetools/cloud/aws-lambda/lambda /lambda

# Passing relative Path of handler to awslambdaric causes path resolution issues, so we set PYTHONPATH to the location of the handler
ENV PYTHONPATH /lambda
ENV JAVA_HOME /etc/alternatives/jre

# Since AWS Lambda has permission to only write in /tmp, we execute extract_scenes in /tmp
WORKDIR /tmp
USER sbx_user1051

ENTRYPOINT [ "/usr/bin/python3", "-m", "awslambdaric" ]
CMD ["lambda.handler"]