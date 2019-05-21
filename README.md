# Description
Quantative Omics Discovery Base is a fast, efficient web-based tool for exploring functional enrichment of expression sets.

ui.R contains the structure of visual components
server.R contains core functions that communicate with ui.R
global.R contains additional core functions and is the link between the server functions and the global environment in R

Additional scripts are found in /bin and each represents a major section of the tool.

# Installation

In order to run this tool, first install

devtools::install_github('martinry/qob')

and download required annotation data from

https://www.dropbox.com/s/gs6l1alc37i1cx2/Large_files.zip?dl=0

Store this folder in your home directory.


Finally, install the QODB in your home directory

git clone https://github.com/martinry/qodb-shiny.git


The software has been tested in RStudio. Open one of the following files in RStudio:
- ui.R
- server.R
- global.R

and click the options for "Run App". Select "Run External". This way the software will run in your browser.

Click Run App to launch.
