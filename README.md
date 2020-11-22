# Description
ProteoMill is a fast, efficient web-based tool for exploring functional enrichment of expression sets.

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

# Run locally
The software has been tested in RStudio. Open one of the following files in RStudio:
- ui.R
- server.R
- global.R

and click the options for "Run App". Select "Run External". This way the software will run in your browser.

Click Run App to launch.
