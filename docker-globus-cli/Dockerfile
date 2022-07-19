# Set the base image to Ubuntu:xenial
FROM ubuntu:xenial


# Update the repository sources list
# Install base packages: pip, git, vim
RUN apt-get update && DEBIAN_FRONTEND=noninteractive apt-get install -y \
    python3-pip \
    git \
    vim

#Obtain and run get-pip.py necessary for Globus CLI install
WORKDIR /opt/
ADD https://bootstrap.pypa.io/get-pip.py .
RUN python3 get-pip.py

#Create environment settings necessary for Globus CLI install
ENV LC_ALL=C.UTF-8
ENV LANG=C.UTF-8
ENV GLOBUS_CLI_INSTALL_DIR="$(python -c 'import site; print(site.USER_BASE)')/bin"
RUN echo "GLOBUS_CLI_INSTALL_DIR=$GLOBUS_CLI_INSTALL_DIR"
ENV PATH="$GLOBUS_CLI_INSTALL_DIR:$PATH"
RUN echo 'export PATH="'"$GLOBUS_CLI_INSTALL_DIR"':$PATH"' >> "$HOME/.bashrc"

#Install Globus CLI
RUN pip install --upgrade globus-cli

#Install virtualenv
RUN pip install virtualenv

#Pull in necessary automation examples
RUN git clone https://github.com/globus/automation-examples
WORKDIR automation-examples
RUN chmod +x cleanup_cache.py cli-sync.sh globus_folder_sync.py share-data.sh share_data.py

#Run and install requirements for virtualenv
RUN virtualenv venv
RUN . venv/bin/activate
RUN pip install -r requirements.txt
