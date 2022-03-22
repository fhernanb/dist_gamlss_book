#
# This is the server logic of a Shiny web application. You can run the
# application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(gamlss)
library(gamlss.dist)
library(dplyr)

source("aux_functions.R")


# Define server logic required to draw a histogram
shinyServer(function(input, output) {
  
  output$nopar <- reactive({
    the_dist <- identifica(input$distribution)
    the_dist$nopar
    })
  
  outputOptions(output, 'nopar', suspendWhenHidden = FALSE)
  
  # output$type <- reactive({
  #   the_dist <- identifica(input$distribution)
  #   the_dist$type
  # })
  # 
  # outputOptions(output, 'type', suspendWhenHidden = FALSE)
  # output$id<-reactive({identifica(input$distribution)$nopar})

  output$distPlot <- renderPlot({
    the_dist <- identifica(dist=input$distribution, input=input$mu)
    if (the_dist$type == "Discrete")
      plot_discrete(input)
    else
      plot_continuous(input)
    })

})
