# Description
ProteoMill is a fast, efficient web-based tool for exploring functional enrichment of expression sets. It is freely available at https://proteomill.com/

# Directory structure
ui.R contains the structure of visual components
server.R contains core functions that communicate with ui.R
global.R contains additional core functions and is the link between the server functions and the global environment in R Shiny

bin/ contains additional scripts, each of which represents some major functionality of the tool.
contact/ is the output directory for the in-app contact form.
data/ contains all public and private (for testing) datasets.
lib/ contains all annotation data from STRING and Reactome
reports/ is a template document used for the "Generate report" feature.
www/ contains various web elements: CSS, javascript, icons, fonts, and a sub-directory, doc/, which contains all in-app modal
displays such as News.

