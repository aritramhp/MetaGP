# For more information, please refer to https://aka.ms/vscode-docker-python
FROM python:3.9.19-slim-bullseye

# Keeps Python from generating .pyc files in the container
ENV PYTHONDONTWRITEBYTECODE=1

# Turns off buffering for easier container logging
ENV PYTHONUNBUFFERED=1

RUN mkdir /app
WORKDIR /app
COPY . /app


RUN apt update
RUN apt upgrade
RUN pip install --upgrade pip
RUN apt-get -y install build-essential
# Install pip requirements
RUN pip install -r /app/requirements.txt
# Install java
RUN apt-get update && \
    apt-get install -y openjdk-11-jre-headless && \
    apt-get clean;

# Creates a non-root user with an explicit UID and adds permission to access the /app folder
# For more info, please refer to https://aka.ms/vscode-docker-python-configure-containers
#RUN adduser -u 5678 --disabled-password --gecos "" appuser && chown -R appuser /app
#USER appuser

# During debugging, this entry point will be overridden. For more information, please refer to https://aka.ms/vscode-docker-python-debug
CMD ["bash", "MetaGP.sh"]