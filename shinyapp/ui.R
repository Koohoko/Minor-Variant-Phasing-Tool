#Copyright (c) 2017 Haogao Gu. All rights reserved.
library(shiny)
#runApp("shinyapp")

# Define UI for miles per gallon application
shinyUI(pageWithSidebar(
    
    # Application title
    headerPanel("Minor Variant Phasing Tool"),
    
    sidebarPanel(
        fileInput('bamfile', 'Choose bam File'),
        fileInput('reffile', 'Choose Reference File'),
        fileInput('snpfile', 'Choose SNPs File'),
        uiOutput("selector"),
        submitButton("Update View"),
        hr(),
        a(href="https://github.com/Koohoko", "Copyright (c) 2017 Haogao Gu. All rights reserved.")
    ),
    
    
    mainPanel(imageOutput("myImage"),
              textOutput('text'))
))