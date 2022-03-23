#
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(gamlss)
library(gamlss.dist)


# Define UI for application that draws a histogram
shinyUI(fluidPage(

    # Application title
    titlePanel("All gamlss distributions"),

    # Sidebar with a slider input for number of bins
    sidebarLayout(
        sidebarPanel(
            
            # Para elegir la distribucion
            selectInput(inputId="distribution",
                        label="Select the gamlss distribution:",
                        choices=grep("^[A-Z ]+$", ls("package:gamlss.dist"), value = TRUE),
                        selected="EXP"),
            
            # Para elegir los valores de los parametros
            conditionalPanel(condition="output.nopar==1",
                             numericInput(inputId="mu",
                                          label=HTML("Select &mu;:"),
                                          min=-20,
                                          max=20,
                                          step=0.1,
                                          value=NA)),
            
            conditionalPanel(condition="output.nopar==2",
                             numericInput(inputId="mu",
                                          label=HTML("Select &mu;:"),
                                          min=-20,
                                          max=20,
                                          step=0.1,
                                          value=NA),
                             
                             numericInput(inputId="sigma",
                                          label=HTML("Select &sigma;:"),
                                          min=-20,
                                          max=20,
                                          step=0.1,
                                          value=NA)),
            
            conditionalPanel(condition="output.nopar==3",
                             numericInput(inputId="mu",
                                          label=HTML("Select &mu;:"),
                                          min=-20,
                                          max=20,
                                          step=0.1,
                                          value=NA),
                             
                             numericInput(inputId="sigma",
                                          label=HTML("Select &sigma;:"),
                                          min=-20,
                                          max=20,
                                          step=0.1,
                                          value=NA),
                             
                             numericInput(inputId="nu",
                                          label=HTML("Select &nu;:"),
                                          min=-20,
                                          max=20,
                                          step=0.1,
                                          value=NA)),
            
            conditionalPanel(condition="output.nopar==4",
                             numericInput(inputId="mu",
                                          label=HTML("Select &mu;:"),
                                          min=-20,
                                          max=20,
                                          step=0.1,
                                          value=NA),
                             
                             numericInput(inputId="sigma",
                                          label=HTML("Select &sigma;:"),
                                          min=-20,
                                          max=20,
                                          step=0.1,
                                          value=NA),
                             
                             numericInput(inputId="nu",
                                          label=HTML("Select &nu;:"),
                                          min=-20,
                                          max=20,
                                          step=0.1,
                                          value=NA),
                             
                             numericInput(inputId="tau",
                                          label=HTML("Select &tau;:"),
                                          min=-20,
                                          max=20,
                                          step=0.1,
                                          value=NA)),
            
            # Para ingresar el maximo y minimo
            sliderInput(inputId = "minimo",
                        label = HTML("Select the minimum value for the random variable:"),
                        min = -30,
                        max =  30,
                        value = 0.01,
                        step= 0.01,
                        animate = TRUE),
            sliderInput(inputId = "maximo",
                        label = HTML("Select the maximum value for the random variable:"),
                        min = -30,
                        max =  30,
                        value = 0.99,
                        step= 0.01,
                        animate = TRUE)
        ),

        # Show a plot of the generated distribution
        mainPanel(
            plotOutput("distPlot"),
            verbatimTextOutput("value"),
            verbatimTextOutput("summary")
        )
    )
))
