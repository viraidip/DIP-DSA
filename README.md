[![DOI](https://zenodo.org/badge/528365333.svg)](https://zenodo.org/doi/10.5281/zenodo.10455047)

# Defective Interfering Particle Deletion Site Analyzer (DIP-DSA)

This application allows the user to investigate datasets that include Defective
Interfering Particles. It gives different analyses about the DI RNAs of these
particles. They contain a long deletion site internally.

## How to set up

As a first step the right R version (3.6.3) and some additional packages need
to be installed. These additional libraries are listed here:

- shiny 1.7.2
- shinydashboard 0.7.2
- ggplot2 3.3.6
- DT 0.24
- plotly 4.10.0
- Biostrings 2.54.0
- hash 2.2.6.2
- jsonlite 1.8.2
- comprehenr 0.6.10
- plyr 1.8.8

The libraries can be installed by running the following command from the root
folder:

```
bash setup_script.sh
```

For further support on the installation of specific libraries please use the
original documentation of the corresponding library.

### Running the app

To run the application locally type the following command into the commandline.

```
bash run_DIP-DSA_local.sh
```

It will automatically open the application on http://localhost:3838.

#### Building a docker container

Alternatively a docker container can be built to deploy the application.
The commands for this are given in docker_commands.sh. It builds the container
and runs it locally.

```
bash docker_commands.sh
```

It will automatically open the application on http://localhost:3838.

#### Pull docker image from docker hub

A prebuild docker image can be pulled form docker hub. The image can be found on
[docker hub](https://hub.docker.com/repository/docker/lohmannjens/dip-dsa/general).
To run it the following steps need to be performed

```
docker pull lohmannjens/dip-dsa
docker run -p 3838:3838 lohmannjens/dip-dsa
```

After this the appication can be accessed on http://localhost:3838.

### Example dataset

When uploading a custom data set the user has to provide the data in a specific
format. It is a csv file with four columns. The order correct ordering of the
columns is crucial! For more information about the content see the following
table:

|                 | Column1 | Column2  | Column3 | Column4      |
|-----------------|---------|----------|---------|--------------|
| header names    | Segment                                              | Start                               |End                                | NGS_read_count                                                 |
| column data type| character                                            | integer                             | integer                           | integer                                                        |
| description     | Name of the segment that the DI RNA is comming from  | Start position of the deletion site | End position of the deletion site | Number of counts in the NGS data of the specific DI RNA sample |

