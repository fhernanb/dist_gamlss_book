identifica <- function(x) {
  d <- gamlss.family(x)
  which_param <- case_when(
    d$nopar == 1 ~ ",mu=input$mu",
    d$nopar == 1 ~ ",mu=input$mu,sigma=input$sigma",
    d$nopar == 1 ~ ",mu=input$mu,sigma=input$sigma,nu=input$nu",
    d$nopar == 1 ~ ",mu=input$mu,sigma=input$sigma,nu=input$nu,tau=input$tau"
  )
  
  list(family=d$family, parameters=d$parameters, nopar=d$nopar,
       disc_cont=d$type, which_param=which_param)
}

