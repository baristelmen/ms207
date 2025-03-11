<!-- PROJECT LOGO -->
<br />
<div align="center">
  
  <h3 align="center">MS207 - Travaux dirigés (TD) #5</h3>

  <p align="center">
    Finite element simulation with CAST3M and Fenicsx/Dolfinx
    <br />
    <a href="https://github.com/baristelmen/MS207/archive/refs/heads/main.zip">Download</a>
    ·
    <a href="https://github.com/baristelmen/MS207/issues">Report Bug</a>
    ·
    <a href="https://github.com/baristelmen/MS207/pulls">Request Feature</a>
  </p>
</div>

## Table of Contents

* [About this repository](#about-this-repository)
  * [Overview](#overview)
  * [Features](#features)
* [Getting started](#getting-started)
  * [Setting up the environment](#setting-up-the-environment)
  * [Using this repo](#using-this-repo)
* [Getting Help](#getting-help)
  * [Issues tab](#1-issues-tab)
  * [Contact maintainer](#contact-maintainer)
* [Why containers?](#why-containers)
* [Contributing](#contributing)
* [License](#license)


## About this repository

### Overview

This repository provides a Docker Compose setup for running a JupyterLab environment with pre-installed CAST3M and jupyterlab. This setup allows you to perform finite element simulations using the CAST3M and Fenics/dolfinx. Essentialy, the compose file will be fetching image specially created for **MEC_4MS07_TA** course which is named,  ``` MS207 ```. 

### Features

- **Jupyter Lab**: Jupyter Lab is a powerful web-based interactive development environment for Python and other programming languages. It comes pre-configured with essential libraries and tools for data analysis, machine learning, and more.

- **KASM Workspace**: KASM allows the creation of custom container images with VNC support for remote desktop access. To build a custom image, use a base image like kasmweb/ubuntu-focal-core and install a window manager (e.g., Openbox) along with x11vnc and novnc. Configure x11vnc to start on container launch and expose the required ports. A proper kasm_user_setup script ensures user settings persist. That kind of image will allow you to use a fully fledged desktop without using VirtualBox variants.

- **Docker Compose**: Docker Compose simplifies the management of multi-container applications. With just a few commands, you can spin up your Jupyter Lab and MFront containers and start working on your projects.


## Getting Started

### Setting up the environment

This repo uses Docker to run its services contrary to classical VMs. Ensure your system meets the specified requirements before proceeding.

#### System Requirements

**Linux:**
- A compatible CPU with hardware virtualization support (for running Docker containers).
- At least 2GB of RAM.
- A network connection

**Windows:**
- Windows 10 64-bit: Home, Pro, or Enterprise Edition.
- Virtualization must be enabled in BIOS/UEFI settings.
- Hyper-V must be disabled if enabled.
- Docker Desktop for Windows requires Microsoft Hyper-V.
- A network connection.

**macOS:**
- A 2010 or newer Mac with macOS 10.13 or later.
- macOS Sierra or later.
- At least 4GB of RAM.
- A network connection.

Docker supports multiple Linux distros at its current stage and provides .deb and .rpm packages. To proceed with the installation process, follow the link below for official instructions.

https://docs.docker.com/engine/install/

This repo uses a docker-compose plugin. Normally, it should be installed via standard docker installation. If not, refer to Docker documentation for proper installation instructions.

## Using this repo

To use this repo, one can either download or clone. We do not provide versions or tags on this repo.

To clone this repository,

```console
git clone https://github.com/baristelmen/MS207.git
```

To download this repository
```console
https://github.com/baristelmen/MS207/archive/refs/heads/main.zip
```
Once you have the source code, follow the steps below (Identical for each OS):

1. Unzip the folder on your computer 
2. Open you favorite terminal on your preferred OS and navigate to the unzipped folder
3. An environment file (.env) is provided in this repository. This file can be either modified or left as it is.
4. At the source code folder, execute the following command
   ```shell
   docker compose up -d
5. This command will download the appropriate image, based on your CPU architecture, on your local machine. The download could take some time. So make sure that you have enough coffee.
6. To access the interface, open up your favorite browser and tap https://localhost:6901 or https://127.0.0.1:6901. This will allow you to access the desktop environment
7. The username is **kasm_user** and the password is defined in docker-compose.yaml (you can change it). By default, the password is **ms207ms207**.
8. The Jupyter lab is installed and mildly configured on your computer. To access your Jupyter lab interface, you have to be in the workspace image (https://localhost:6901 or https://127.0.0.1:6901) and you can launch it via the browser inside the image.
   ```shell
   http://localhost:8888/lab
9. There is no password associated with Jupyter Lab because it runs inside the container. 
10. The data stored in the docker images can be accessed from your computer via the data folder given in this repository. The interactive sharing folder definition is pre-configured in **docker-compose.yaml** file. You can always access to your data from this folder including your simulation results or your custom Python codes. 

## Getting Help

If you encounter any issues or have questions while using this environment, here are a few ways you can seek assistance:

#### **Issues Tab**

For bugs, feature requests, or general questions, please open an issue in the [Issues tab](#). Be sure to provide as much detail as possible, including:
- Steps to reproduce the issue.
- Any error messages or logs.
- Your environment details (e.g., Docker version, system configuration).
 
#### **Contact Maintainer**

If your issue is specific to this repository, feel free to contact the maintainer directly by creating an issue or posting a comment. 

## Why containers?

Containers have become a popular solution for developing, deploying, and managing applications because they offer numerous advantages over traditional approaches, such as virtual machines or bare-metal setups. Here’s why using containers is a great choice for this project and in general:

#### **Consistency Across Environments**
One of the key benefits of containers is that they ensure consistency across different environments. Whether you're running your application locally, in a staging environment, or in production, containers ensure that:
- The application behaves the same way regardless of the underlying operating system or infrastructure.
- Dependencies, libraries, and system settings are included within the container, minimizing "it works on my machine" issues.

In this project, containers ensure that the TFEL/MFront installation and JupyterLab environment behave identically on any machine where Docker is available.

#### **Simplified Setup and Deployment**
Containers make it easy to package an entire application, including its dependencies, into a single portable unit. This simplifies the setup process by allowing users to:
- Avoid the hassle of installing and configuring TFEL/MFront or other dependencies manually.
- Launch the entire environment with a single `docker-compose up` command.

This project leverages Docker Compose to automatically set up the environment, meaning users don't have to worry about complex installation instructions.

#### **Isolation and Security**
Containers provide a high level of isolation between applications running on the same system. Each container has its own filesystem, network interfaces, and resource limits, reducing the risk of conflicts between applications. Benefits include:
- **Resource Isolation**: You can limit the CPU and memory usage of containers, ensuring that this environment does not interfere with other services running on the same machine.
- **Security**: Containers are isolated from the host system, and potential vulnerabilities within the application are less likely to affect the broader system.

In this project, the JupyterLab environment with TFEL/MFront is completely isolated from the host, ensuring that any changes or experiments conducted inside the container do not affect your local setup.

#### **Portability**
Containers make your application highly portable. Since Docker containers can run anywhere Docker is installed, you can:
- Move the application seamlessly between local development environments, cloud services, or even between different operating systems.
- Share the exact same environment with collaborators without worrying about platform differences or missing dependencies.

For this project, this means you can share your simulation environment, ensuring others can run the same TFEL/MFront simulations without additional setup.

#### **Efficient Resource Usage**
Unlike traditional virtual machines, containers are lightweight. They share the host system's kernel, which significantly reduces overhead compared to running a full-fledged VM. As a result:
- Containers start much faster than virtual machines.
- They use less disk space and memory, making them more efficient in handling large-scale applications or microservices.

This project benefits from the lightweight nature of containers, allowing users to run simulations without heavy resource overhead.

#### **Version Control and Reproducibility**
With containers, you can specify exactly which versions of software dependencies to use in the `Dockerfile` or `docker-compose.yml`. This leads to:
- **Reproducibility**: Anyone using the same container image will get the same results, which is particularly important for scientific computing and simulation environments.
- **Version Control**: You can version your container images, ensuring that older simulations can be rerun in the same environment, even if the software stack has since been updated.

#### **Scalability**
Containers are designed to be stateless and lightweight, making them ideal for scaling applications across multiple servers. Using orchestration tools like Kubernetes or Docker Swarm, you can:
- Easily scale up or down based on demand.
- Deploy containerized applications in a cluster, handling high workloads with efficiency.

Although this project may not require large-scale deployment, the containerized setup ensures it can scale if needed, such as running multiple simulations in parallel.

#### **Collaboration and Sharing**
Containers are highly portable and easy to distribute. You can:
- Share your work with others by distributing the Docker image, enabling your colleagues or community to run the exact same environment on their machines.
- Use container registries like Docker Hub to publish pre-built images, allowing anyone to pull and run the container without building it themselves.

## Contributing
Contributions to this repository are welcome! If you encounter issues or have improvements to suggest, please feel free to open an issue or submit a pull request.

## License
This repository is open-source and distributed under the MIT License. You are free to use, modify, and distribute it according to the terms of the license.
