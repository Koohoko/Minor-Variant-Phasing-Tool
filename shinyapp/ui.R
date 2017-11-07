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
        uiOutput("selector"),
        wellPanel(
        checkboxInput("usesnp","Use Reference SNPs",F),
        fileInput('snpfile', 'Choose SNPs File')
        ),
        #actionButton("goButton", "Update View",icon = icon("refresh")),
        submitButton("Update View"),
        hr(),
        a(href="https://github.com/Koohoko/Minor-Variant-Phasing-Tool",
          "Copyright (c) 2017 Haogao Gu. All rights reserved.")
    ),
    
    
    mainPanel(imageOutput("myImage"),
              textOutput('text'))
))