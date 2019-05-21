Quantative Omics Discovery Base is a fast, efficient web-based tool for exploring functional enrichment of expression sets.

ui.R contains the structure of visual components
server.R contains core functions that communicate with ui.R
global.R contains additional core functions and is the link between the server functions and the global environment in R

Additional scripts are found in /bin and each represents a major section of the tool.

In order to run this tool, first install

devtools::install_github('martinry/qob')


