# Defective Interfering Particle Deletion Site Analyzer (DIP-DSA)

This application allows the user to investigate datasets that include Defective
Interfering Particles. 

## How to set up

As a first step the right R version (3.6.3) and some additional packages need to be 
installed. These libraries are listed here:

- shiny 1.7.2
- shinydashboard 0.7.2
- ggplot2 3.3.6
- DT 0.24
- plotly 4.10.0
- Biostrings 2.54.0
- hash 2.2.6.2
- jsonlite 1.8.2
- comprehenr 0.6.10

For support on the installation of specific libraries please use the original
documentation of the corresponding library.

### Adding data

To have access to the predefined datasets the user has to download them
[here](https://ecosia.de "data source").
After the download has finished move the ```data``` folder to the root folder of the
repository. The root folder should look like this afterwards:

```
+DIP-DSA
|
+--README.md
+--run_DIP-DSA_local.sh
+--setup_script.sh
+--data/
+--src/
```

### Running the app

To run the application locally type the following command into the commandline.

```
bash run_DIP-DSA_local.sh
```

It will automatically open the application on http://localhost:1234
