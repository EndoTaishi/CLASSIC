Singularity Containers {#Containers}
============


**This is from Ed Wisernig. It needs editing before being added to a repo**

# What are containers?

A recent journal paper describes containers as "a software-based packaging and distribution tool that collects all elements and dependencies of a Linux-based application. Containers store the runtime environment together with one or more applications for ease of transportation, installation, and execution across a variety of operating systems (Linux, Mac, Windows)." (for more information, please see the References section)
Docker is the container type we will be using. Singularity is the platform that will allow us to build, run or shell into one of these Docker containers.

For example, in order to run CLASSIC model, we need several specific tools such as compilers (GNU, Intel etc.), libraries (MPI, NetCDF etc.) that have to be present on our machine to be able to run the model.
If, however, we use containers instead, then we simply access the model through a container, a process that eliminates the cumbersome and time consuming process of installing all of these compilers and libraries locally.

Another significant advantage is that the version of these libraries is "frozen" in the container. For example, we could have a container with the library versions at the time of a major model release version. So, even several years later, we can run the model with the original intended library versions.

# How to use Singularity containers?
In order to use Singularity containers, one must first make certain that a local installation of Singularity is available.

On a Linux machine (Ubuntu in our particular case) one may use the following command to install singularity.

`sudo apt install singularity-container`

Please see this link for more detailed instructions on how to set up Singularity on Mac and Windows:

[http://singularity.lbl.gov/singularity-tutorial](http://singularity.lbl.gov/singularity-tutorial)

# Hello World
Once singularity has been installed, let's try out a simple hello world example:

`singularity run shub://vsoch/hello-world`

As we can see in the above example,

# Creating Singularity containers

On Linux, the easiest way to install singularity is by using the standard repo:

`sudo apt-get install -y singularity-container`

Then, download our container: **FLAG need to make our own**

`singularity pull shub://eduardwisernig/testSingularity`

And shell into it:

`singularity shell eduardwisernig-testSingularity-master.simg`

One can successfully set up singularity on a macOS platform using vagrant using the method described below:

`http://singularity.lbl.gov/install-mac#option-1-singularityware-vagrant-box`

For more information about Virtualbox, please visit:

`https://www.virtualbox.org/wiki/VirtualBox`

For more information about Vagrant, please visit:

`https://www.vagrantup.com/intro/index.html`

Regardless of wether you have your own immediate environment set up or if you are shelled into our container, now you can get ready to start using the model.


## Github
## Singularity-hub


# Disclaimer
Some of the text used in this file is copied and pasted in from journal papers or other sources. While not referenced individually, all the original sources can be found in the *References* section.

# References
`http://journals.ametsoc.org/doi/pdf/10.1175/BAMS-D-15-00255.1`
