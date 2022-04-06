# Robert Rigby and Mikis Stasinopoulos
# 28-2-18
# to DO 
#  i) binomial model???
# -------------------------------------------------------------------------
gamlssRobust1 <- function(object, 
                          bound = abs(qNO(1/(2*n))), 
                         n.iter = 20, 
                         c.crit = 0.001, 
                          trace = FALSE)
{
#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
# SUMMARY
# I) set m1 to object
# II) for i to n.iter
#  a) get deviance and residual from m1
#  b) using residual find bound 1 and 2 for y
#  c) bound y and refit m1 
#  d) break if |olddeviance-newdeviance| is small
# III) calculate the new residuals from the the fitted model but old y 
#  VI) output m1 with yRobust resRobust and new residuals
#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------  
# this is to replicate rqres within gamlssMLpred enviroment
# it is used as in gamlss()
  rqres <- function (pfun = "pNO", 
                     type = c("Continuous", "Discrete", "Mixed"),
                     censored = NULL,  
                     ymin = NULL, 
                     mass.p = NULL, 
                     prob.mp = NULL,
                     y = y,
                     ... )
  { }
  body(rqres) <-  eval(quote(body(rqres)), envir = getNamespace("gamlss"))
  ##---------------------------------------------------------------------------------------
  ##---------------------------------------------------------------------------------------    
# MAIN FUNCTION STARTS HERE
# get the family
   fname <- as.gamlss.family(object$family) 
    pfun <- paste("p", fname$family[[1]],sep="")
    qfun <- paste("q", fname$family[[1]],sep="")
    lpar <- length(fname$parameters)
   dtype <- fname$type  
   if (dtype!="Continuous") warning("the gamlssRobust is not tested for discrete distributions")
   if (fname$family[[1]] %in% .gamlss.bi.list) stop("No robust versions of binomial models exist")
# rename object  
      m1 <- object
      n <- length(object$y)
# get  y for new data  if binomial       
  if (fname$family[1] %in% .gamlss.bi.list) bd <- object$bd
       y <- m1$y
#---------------------------------------------
#---------------------------------------------       
# iterations  Q
for (i in 1:n.iter)# starts loop
  {
    olddev <- deviance(m1) # old deviance
    # abs(k) is the bound on the normalized quantile residuals
        k <- rep(-bound, length(y))
       pk <- pNO(k)
      res <- m1$residuals # old residuals
#---------------------------------------------
# finding bound 1 and two to bound the y variable      
  if(lpar==1) 
    {
      if (fname$family[[1]] %in% .gamlss.bi.list)
      {
        bound1 <-  eval(call(qfun, p=  pk, mu=fitted(m1, "mu"), bd=bd))
        bound2 <-  eval(call(qfun, p=1-pk, mu=fitted(m1, "mu"), bd=bd))
      } else
      {
        bound1 <-  eval(call(qfun, p=  pk, mu=fitted(m1, "mu")))
        bound2 <-  eval(call(qfun, p=1-pk, mu=fitted(m1, "mu")))
      } 
  }
    if(lpar==2) 
    {
      if (fname$family[[1]] %in% .gamlss.bi.list)
      {
        bound1 <-  eval(call(qfun, p=  pk, mu=fitted(m1, "mu"),sigma=fitted(m1, "sigma"), bd=bd))
        bound2 <-  eval(call(qfun, p=1-pk, mu=fitted(m1, "mu"),sigma=fitted(m1, "sigma"), bd=bd))
      } else
      {
        bound1 <-  eval(call(qfun, p=  pk, mu=fitted(m1, "mu"),sigma=fitted(m1, "sigma")))
        bound2 <-  eval(call(qfun, p=1-pk, mu=fitted(m1, "mu"),sigma=fitted(m1, "sigma")))
      } 
    }
    if(lpar==3) 
    {
      if (fname$family[[1]] %in% .gamlss.bi.list)
      {
        bound1 <-  eval(call(qfun, p=pk, mu=fitted(m1, "mu"), sigma=fitted(m1, "sigma"), nu=fitted(m1, "nu"), bd=bd))
        bound2 <-  eval(call(qfun, p=1-pk, mu=fitted(m1, "mu"), sigma=fitted(m1, "sigma"), nu=fitted(m1, "nu"), bd=bd))
      } else
      {
        bound1 <-  eval(call(qfun, p=pk, mu=fitted(m1, "mu"), sigma=fitted(m1, "sigma"), nu=fitted(m1, "nu")))
        bound2 <-  eval(call(qfun, p=1-pk, mu=fitted(m1, "mu"), sigma=fitted(m1, "sigma"), nu=fitted(m1, "nu")))
      }
    }
    if(lpar==4) 
    {
      if (fname$family[[1]] %in% .gamlss.bi.list)
      {
        bound1 <-  eval(call(qfun, p=pk,   mu=fitted(m1, "mu"), sigma=fitted(m1, "sigma"), nu=fitted(m1, "nu"), tau=fitted(m1, "tau"), bd=bd))
        bound2 <-  eval(call(qfun, p=1-pk, mu=fitted(m1, "mu"), sigma=fitted(m1, "sigma"), nu=fitted(m1, "nu"), tau=fitted(m1, "tau"), bd=bd))
      } else
      {
        bound1 <-  eval(call(qfun, p=pk,   mu=fitted(m1, "mu"), sigma=fitted(m1, "sigma"), nu=fitted(m1, "nu"), tau=fitted(m1, "tau")))
        bound2 <-  eval(call(qfun, p=1-pk, mu=fitted(m1, "mu"), sigma=fitted(m1, "sigma"), nu=fitted(m1, "nu"), tau=fitted(m1, "tau")))
      } 
    }
#---------------------------------------------    
# y is now bounded above and below
    y1 <- ifelse(res > k,     y,   bound1)
    y1 <- ifelse(res < (-k),  y1,  bound2)
    if (trace) cat("------", i, "\n")
    #
     assign("y1", y1, envir=.GlobalEnv)
     on.exit(rm(y1, envir=.GlobalEnv))
     m1 <- update(m1, formula=update(formula(m1), y1~.), start.from=m1, trace=trace)
    # m1$y <- y1
    # m1 <- update(m1,formula=update(formula(m1), m1$y~.), start.from=m1, trace=trace)
    newdev <- deviance(m1)
    if (abs(newdev-olddev)< c.crit)  break
} # finish loop
#-------------------------------------------------------------
#-------------------------------------------------------------       
# get the new residuals-------------------------------------------
# jump depending on the number of parameters 
if(lpar==1) 
{
  if (fname$family[[1]] %in% .gamlss.bi.list)
  {
    ures <-  call("rqres", pfun=pfun, type=dtype,
                  ymin=fname$rqres[[1]][["ymin"]], 
                  y=y, bd=bd, mu= fitted(m1))   
  } else
  {
    ures <- call("rqres", pfun=pfun, type=dtype, 
                 ymin=fname$rqres[[1]][["ymin"]], 
                 y=y, mu= fitted(m1))   
  } 
}
else if(lpar==2)
{
  if (fname$family[[1]] %in% .gamlss.bi.list)
  {
    ures <-  call("rqres", pfun=pfun, type=dtype, 
                  ymin=fname$rqres[[1]][["ymin"]], 
                  y=y, mu= fitted(m1), sigma= fitted(m1,"sigma"), bd=bd)
  } else
  {
    ures <- call("rqres", pfun=pfun,  type=dtype, 
                 ymin=fname$rqres[[1]][["ymin"]], 
                 y=y, mu= fitted(m1), sigma= fitted(m1,"sigma"))
  } 
}
else if(lpar==3)
{
  if (fname$family[[1]] %in% .gamlss.bi.list)
  {
    ures <- call("rqres", pfun=pfun,  type=dtype, 
                 ymin=fname$rqres[[1]][["ymin"]], 
                 y=y, mu= fitted(m1), sigma= fitted(m1,"sigma"),
                 nu =  fitted(m1, "nu"), bd=bd)
  } else
  {
    ures <- call("rqres", pfun=pfun,  type=dtype, 
                 ymin=fname$rqres[[1]][["ymin"]], 
                 y=y, mu= fitted(m1), sigma= fitted(m1,"sigma"),
                 nu =  fitted(m1, "nu"))
  } 
}
else 
{
  if (fname$family[[1]] %in% .gamlss.bi.list)
  {
    ures <- call("rqres", pfun=pfun,  type=dtype, 
                 ymin=fname$rqres[[1]][["ymin"]], 
                 y=y, mu= fitted(m1), sigma= fitted(m1,"sigma"),
                 nu =  fitted(m1, "nu"), tau = fitted(m1, "tau"), bd=bd)
  } else
  {
    ures <- call("rqres", pfun=pfun,  type=dtype, 
                 ymin=fname$rqres[[1]][["ymin"]], 
                 y=y, mu= fitted(m1), sigma= fitted(m1,"sigma"),
                 nu =  fitted(m1, "nu"), tau = fitted(m1, "tau"))
  } 
}
       # get output   
        m1$yRobust <- m1$y
              m1$y <- y       
m1$residualsRobust <- resid(m1)
      m1$residuals <- eval(ures)
          m1$bound <- bound
m1
}  
#--------------------------------------------------------------------------
#--------------------------------------------------------------------------
#--------------------------------------------------------------------------
#--------------------------------------------------------------------------
#--------------------------------------------------------------------------
# -------------------------------------------------------------------------
gamlssRobustW <- function(object, 
                         bound = abs(qNO(1/(2*n))),
                      CD.bound = abs(qNO(1/(5*n))), 
                         n.iter = 20, 
                         c.crit = 0.001, 
                         trace = FALSE)
{
  #--------------------------------------------------------------------------------------
  #--------------------------------------------------------------------------------------
  # SUMMARY
  # I) set m1 to object and get original weights
  # II) for j to n.iter loop 1 
  #     for i to n.iter loop 2 
  #  a) get deviance and residual from m1
  #  b) using residual find bound 1 and 2 for y
  #  c) bound y and refit m1 
  #  d) break loop 2 if |olddeviance-newdeviance| is small
  #  f) calculate the new residuals from the the fitted model but old y 
  #  g) set weights euql to zero if resiadual are big
  #  h) break loop 1  if |olddeviance-newdeviance| is small
  #  VI) output m1 
  #--------------------------------------------------------------------------------------
  #--------------------------------------------------------------------------------------  
  # this is to replicate rqres within gamlssMLpred enviroment
  # it is used as in gamlss()
  rqres <- function (pfun = "pNO", 
                     type = c("Continuous", "Discrete", "Mixed"),
                     censored = NULL,  
                     ymin = NULL, 
                     mass.p = NULL, 
                     prob.mp = NULL,
                     y = y,
                     ... )
  { }
  body(rqres) <-  eval(quote(body(rqres)), envir = getNamespace("gamlss"))
  ##---------------------------------------------------------------------------------------
  ##---------------------------------------------------------------------------------------    
  # MAIN FUNCTION STARTS HERE
  # get the family
  fname <- as.gamlss.family(object$family) 
   pfun <- paste("p", fname$family[[1]],sep="")
   qfun <- paste("q", fname$family[[1]],sep="")
   lpar <- length(fname$parameters)
  dtype <- fname$type  
  if (dtype!="Continuous") warning("the gamlssRobust is not tested for discrete distributions")
  if (fname$family[[1]] %in% .gamlss.bi.list) stop("No robust versions of binomial models exist")
  # rename object  
    m1 <- object
orig.weights <- m1$weights
   wei <- rep(1, length(m1$y))
   wei <- orig.weights * wei
     n <- length(object$y)
    k2 <- rep(-CD.bound, length(y))
  # get  y for new data  if binomial       
  if (fname$family[1] %in% .gamlss.bi.list) bd <- object$bd
     y <- m1$y
#---------------------------------------------
#---------------------------------------------   
loop2 <-  function(n.iter, c.crit , trace  ) 
     {
      olddev <- deviance(m1)
       for (i in 1:n.iter) # starts loop 2
       {
              k <- rep(-bound, length(y))
             pk <- pNO(k)
            res <- m1$residuals # old residuals
         #---------------------------------------------
         # finding bound 1 and two to bound the y variable      
         if(lpar==1) 
         {
           if (fname$family[[1]] %in% .gamlss.bi.list)
           {
             bound1 <-  eval(call(qfun, p=  pk, mu=fitted(m1, "mu"), bd=bd))
             bound2 <-  eval(call(qfun, p=1-pk, mu=fitted(m1, "mu"), bd=bd))
           } else
           {
             bound1 <-  eval(call(qfun, p=  pk, mu=fitted(m1, "mu")))
             bound2 <-  eval(call(qfun, p=1-pk, mu=fitted(m1, "mu")))
           } 
         }
         if(lpar==2) 
         {
           if (fname$family[[1]] %in% .gamlss.bi.list)
           {
             bound1 <-  eval(call(qfun, p=  pk, mu=fitted(m1, "mu"),sigma=fitted(m1, "sigma"), bd=bd))
             bound2 <-  eval(call(qfun, p=1-pk, mu=fitted(m1, "mu"),sigma=fitted(m1, "sigma"), bd=bd))
           } else
           {
             bound1 <-  eval(call(qfun, p=  pk, mu=fitted(m1, "mu"),sigma=fitted(m1, "sigma")))
             bound2 <-  eval(call(qfun, p=1-pk, mu=fitted(m1, "mu"),sigma=fitted(m1, "sigma")))
           } 
         }
         if(lpar==3) 
         {
           if (fname$family[[1]] %in% .gamlss.bi.list)
           {
             bound1 <-  eval(call(qfun, p=pk, mu=fitted(m1, "mu"), sigma=fitted(m1, "sigma"), nu=fitted(m1, "nu"), bd=bd))
             bound2 <-  eval(call(qfun, p=1-pk, mu=fitted(m1, "mu"), sigma=fitted(m1, "sigma"), nu=fitted(m1, "nu"), bd=bd))
           } else
           {
             bound1 <-  eval(call(qfun, p=pk, mu=fitted(m1, "mu"), sigma=fitted(m1, "sigma"), nu=fitted(m1, "nu")))
             bound2 <-  eval(call(qfun, p=1-pk, mu=fitted(m1, "mu"), sigma=fitted(m1, "sigma"), nu=fitted(m1, "nu")))
           }
         }
         if(lpar==4) 
         {
           if (fname$family[[1]] %in% .gamlss.bi.list)
           {
             bound1 <-  eval(call(qfun, p=pk,   mu=fitted(m1, "mu"), sigma=fitted(m1, "sigma"), nu=fitted(m1, "nu"), tau=fitted(m1, "tau"), bd=bd))
             bound2 <-  eval(call(qfun, p=1-pk, mu=fitted(m1, "mu"), sigma=fitted(m1, "sigma"), nu=fitted(m1, "nu"), tau=fitted(m1, "tau"), bd=bd))
           } else
           {
             bound1 <-  eval(call(qfun, p=pk,   mu=fitted(m1, "mu"), sigma=fitted(m1, "sigma"), nu=fitted(m1, "nu"), tau=fitted(m1, "tau")))
             bound2 <-  eval(call(qfun, p=1-pk, mu=fitted(m1, "mu"), sigma=fitted(m1, "sigma"), nu=fitted(m1, "nu"), tau=fitted(m1, "tau")))
           } 
         }
         #---------------------------------------------    
         # y is now bounded above and below
            y1 <- ifelse(res > k,     y,   bound1)
            y1 <- ifelse(res < (-k),  y1,  bound2)
         if (trace) cat("------", i, "\n")
         #
            assign("y1", y1, envir=.GlobalEnv)
            on.exit(rm(y1, envir=.GlobalEnv))
            mod <- update(m1, formula=update(formula(m1), y1~.), weights = wei, start.from=m1, trace=trace)
          #m1$y <- y1
          #m1 <- update(m1,formula=update(formula(m1), m1$y~.), start.from=m1, trace=trace)
            newdev <- deviance(mod)
         if (abs(newdev-olddev) < c.crit)  break
            olddev <-  newdev
       } # finish loop 1 here     
mod       
}
#---------------------------------------------
#---------------------------------------------     
# loop  1
for (j in 1:n.iter) # loop 1
{
    assign("wei", wei, envir=.GlobalEnv)
    on.exit(rm(wei, envir=.GlobalEnv))
  #-------------------------------------------------------------
  m1 <-  loop2(n.iter=n.iter, c.crit=c.crit , trace=trace )
  #-------------------------------------------------------------       
  # get the new residuals---------------------------------------
  # jump depending on the number of parameters 
  if(lpar==1) 
    {
      if (fname$family[[1]] %in% .gamlss.bi.list)
      {
        ures <-  call("rqres", pfun=pfun, type=dtype,
                      ymin=fname$rqres[[1]][["ymin"]], 
                      y=y, bd=bd, mu= fitted(m1))   
      } else
      {
        ures <- call("rqres", pfun=pfun, type=dtype, 
                     ymin=fname$rqres[[1]][["ymin"]], 
                     y=y, mu= fitted(m1))   
      } 
    }
    else if(lpar==2)
    {
      if (fname$family[[1]] %in% .gamlss.bi.list)
      {
        ures <-  call("rqres", pfun=pfun, type=dtype, 
                      ymin=fname$rqres[[1]][["ymin"]], 
                      y=y, mu= fitted(m1), sigma= fitted(m1,"sigma"), bd=bd)
      } else
      {
        ures <- call("rqres", pfun=pfun,  type=dtype, 
                     ymin=fname$rqres[[1]][["ymin"]], 
                     y=y, mu= fitted(m1), sigma= fitted(m1,"sigma"))
      } 
    }
    else if(lpar==3)
    {
      if (fname$family[[1]] %in% .gamlss.bi.list)
      {
        ures <- call("rqres", pfun=pfun,  type=dtype, 
                     ymin=fname$rqres[[1]][["ymin"]], 
                     y=y, mu= fitted(m1), sigma= fitted(m1,"sigma"),
                     nu =  fitted(m1, "nu"), bd=bd)
      } else
      {
        ures <- call("rqres", pfun=pfun,  type=dtype, 
                     ymin=fname$rqres[[1]][["ymin"]], 
                     y=y, mu= fitted(m1), sigma= fitted(m1,"sigma"),
                     nu =  fitted(m1, "nu"))
      } 
    }
    else if(lpar==4)
    {
      if (fname$family[[1]] %in% .gamlss.bi.list)
      {
        ures <- call("rqres", pfun=pfun,  type=dtype, 
                     ymin=fname$rqres[[1]][["ymin"]], 
                     y=y, mu= fitted(m1), sigma= fitted(m1,"sigma"),
                     nu =  fitted(m1, "nu"), tau = fitted(m1, "tau"), bd=bd)
      } else
      {
        ures <- call("rqres", pfun=pfun,  type=dtype, 
                     ymin=fname$rqres[[1]][["ymin"]], 
                     y=y, mu= fitted(m1), sigma= fitted(m1,"sigma"),
                     nu =  fitted(m1, "nu"), tau = fitted(m1, "tau"))
      } 
      
    } 
    m1$residuals <- eval(ures)
    wei <- ifelse((abs( m1$residuals) > abs(CD.bound)),  0, wei)  
} # finish loop 1
  m1$yRobust <- m1$y
  m1$y <- y       
  m1$residualsRobust <- resid(m1)
  m1$residuals <- eval(ures)
  m1$bound <- bound
  m1$CD.bound <- CD.bound  
  m1
}  


#--------------------------------------------------------------------------
#--------------------------------------------------------------------------
#--------------------------------------------------------------------------
#--------------------------------------------------------------------------
#--------------------------------------------------------------------------
gamlssRobust <- function(object, 
                         bound = abs(qNO(1/n)), 
                      CD.bound = abs(qNO(1/(3*n))), 
                        n.iter = 20, 
                        c.crit = 0.001, 
                         trace = FALSE,
                    weight.out = c("once", "iterative", "prior")) 
{ # 'iterative" the  weights are updated with the iterative ones. Which meens that
  # if an observation is weighted out in the begining remains weighted out).
  # This is more stable
  # "prior" the prior weights are used for udating very unstable
  #--------------------------------------------------------------------------------------
  #--------------------------------------------------------------------------------------
  # this is to replicate rqres within gamlssMLpred enviroment
  # it is used as in gamlss()
  rqres <- function (pfun = "pNO", 
                     type = c("Continuous", "Discrete", "Mixed"),
                     censored = NULL,  
                     ymin = NULL, 
                     mass.p = NULL, 
                     prob.mp = NULL,
                     y = y,
                     ... )
  { }
  body(rqres) <-  eval(quote(body(rqres)), envir = getNamespace("gamlss"))
  ##---------------------------------------------------------------------------------------
  ##---------------------------------------------------------------------------------------    
  # MAIN FUNCTION STARTS HERE
  # get the family
  weight.out <- match.arg(weight.out) 
  fname <- as.gamlss.family(object$family[[1]]) 
   pfun <- paste("p", fname$family[[1]],sep="")
   qfun <- paste("q", fname$family[[1]],sep="")
   lpar <- length(fname$parameters)
  dtype <- fname$type  
  if (dtype!="Continuous") warning("the gamlssRobust is not tested for discrete distributions")
  if (fname$family[[1]] %in% .gamlss.bi.list) stop("No robust versions of binomial models exist")
  # rename object  
     m1 <- object 
   # wei <- object$weights
      n <- length(object$y)
   # get  y for new data  if binomial       
  if (fname$family[1] %in% .gamlss.bi.list) bd <- object$bd
      y <- m1$y   
     k2 <- rep(-CD.bound, length(y))
    wei <- ifelse((abs(resid(object)) > abs(CD.bound)),  0, object$weights)  
  for (i in 1:n.iter)# starts loop
  {
      olddev <- deviance(m1)
    # abs(k) is the bound on the normalized quantile residuals
          k <- rep(-bound, length(y))
         pk <- pNO(k)
        res <- m1$residuals# is this correct?
    #     
    #---------------------------------------------
    if(lpar==1) 
    {
      if (fname$family[[1]] %in% .gamlss.bi.list)
      {
        bound1 <-  eval(call(qfun, p=  pk, mu=fitted(m1, "mu"), bd=bd))
        bound2 <-  eval(call(qfun, p=1-pk, mu=fitted(m1, "mu"), bd=bd))
      } else
      {
        bound1 <-  eval(call(qfun, p=  pk, mu=fitted(m1, "mu")))
        bound2 <-  eval(call(qfun, p=1-pk, mu=fitted(m1, "mu")))
      } 
    }
    if(lpar==2) 
    {
      if (fname$family[[1]] %in% .gamlss.bi.list)
      {
        bound1 <-  eval(call(qfun, p=  pk, mu=fitted(m1, "mu"),sigma=fitted(m1, "sigma"), bd=bd))
        bound2 <-  eval(call(qfun, p=1-pk, mu=fitted(m1, "mu"),sigma=fitted(m1, "sigma"), bd=bd))
      } else
      {
        bound1 <-  eval(call(qfun, p=  pk, mu=fitted(m1, "mu"),sigma=fitted(m1, "sigma")))
        bound2 <-  eval(call(qfun, p=1-pk, mu=fitted(m1, "mu"),sigma=fitted(m1, "sigma")))
      } 
    }
    if(lpar==3) 
    {
      if (fname$family[[1]] %in% .gamlss.bi.list)
      {
        bound1 <-  eval(call(qfun, p=pk, mu=fitted(m1, "mu"), sigma=fitted(m1, "sigma"), nu=fitted(m1, "nu"), bd=bd))
        bound2 <-  eval(call(qfun, p=1-pk, mu=fitted(m1, "mu"), sigma=fitted(m1, "sigma"), nu=fitted(m1, "nu"), bd=bd))
      } else
      {
        bound1 <-  eval(call(qfun, p=pk, mu=fitted(m1, "mu"), sigma=fitted(m1, "sigma"), nu=fitted(m1, "nu")))
        bound2 <-  eval(call(qfun, p=1-pk, mu=fitted(m1, "mu"), sigma=fitted(m1, "sigma"), nu=fitted(m1, "nu")))
      }
    }
    if(lpar==4) 
    {
      if (fname$family[[1]] %in% .gamlss.bi.list)
      {
        bound1 <-  eval(call(qfun, p=pk,   mu=fitted(m1, "mu"), sigma=fitted(m1, "sigma"), nu=fitted(m1, "nu"), tau=fitted(m1, "tau"), bd=bd))
        bound2 <-  eval(call(qfun, p=1-pk, mu=fitted(m1, "mu"), sigma=fitted(m1, "sigma"), nu=fitted(m1, "nu"), tau=fitted(m1, "tau"), bd=bd))
      } else
      {
        bound1 <-  eval(call(qfun, p=pk,   mu=fitted(m1, "mu"), sigma=fitted(m1, "sigma"), nu=fitted(m1, "nu"), tau=fitted(m1, "tau")))
        bound2 <-  eval(call(qfun, p=1-pk, mu=fitted(m1, "mu"), sigma=fitted(m1, "sigma"), nu=fitted(m1, "nu"), tau=fitted(m1, "tau")))
      } 
    }
    #---------------------------------------------    
    # y is now bounded above and below
    y1 <- ifelse(res > -bound,     y,   bound1)
    y1 <- ifelse(res < bound,  y1,  bound2)
  # wei <- ifelse(res > k2,     wei1, 0)
  # wei <- ifelse(res < (-k2),  wei1, 0)
   if (weight.out == "iterative") 
     wei <- ifelse((abs(res) > abs(CD.bound)),  0, wei)  
   if (weight.out == "prior")     
     wei <- ifelse((abs(res) > abs(CD.bound)),  0,  object$weights)  
   #cat(wei, "\n")
    if (trace) cat("------", i, "\n")
    #y1 <<- y 
    assign("y1", y1, envir=.GlobalEnv)
    assign("wei", wei, envir=.GlobalEnv)
    on.exit(rm(y1, wei, envir=.GlobalEnv))
    m1 <- update(m1, formula=update(formula(m1), y1~.), weights=wei, start.from=m1, trace=trace)
    # m1$y <- y1
    # m1 <- update(m1,formula=update(formula(m1), m1$y~.), start.from=m1, trace=trace)
    newdev <- deviance(m1)
    if (abs(newdev-olddev)< c.crit)  break
  } # finish loop
    
  # get the residuals--------------------------------------------
  # jump depending on the number of parameters 
  if(lpar==1) 
  {
    if (fname$family[[1]] %in% .gamlss.bi.list)
    {
      ures <-  call("rqres", pfun=pfun, type=dtype,
                    ymin=fname$rqres[[1]][["ymin"]], 
                    y=y, bd=bd, mu= fitted(m1))   
    } else
    {
      ures <- call("rqres", pfun=pfun, type=dtype, 
                   ymin=fname$rqres[[1]][["ymin"]], 
                   y=y, mu= fitted(m1))   
    } 
  }
  else if(lpar==2)
  {
    if (fname$family[[1]] %in% .gamlss.bi.list)
    {
      ures <-  call("rqres", pfun=pfun, type=dtype, 
                    ymin=fname$rqres[[1]][["ymin"]], 
                    y=y, mu= fitted(m1), sigma= fitted(m1,"sigma"), bd=bd)
    } else
    {
      ures <- call("rqres", pfun=pfun,  type=dtype, 
                   ymin=fname$rqres[[1]][["ymin"]], 
                   y=y, mu= fitted(m1), sigma= fitted(m1,"sigma"))
    } 
  }
  else if(lpar==3)
  {
    if (fname$family[[1]] %in% .gamlss.bi.list)
    {
      ures <- call("rqres", pfun=pfun,  type=dtype, 
                   ymin=fname$rqres[[1]][["ymin"]], 
                   y=y, mu= fitted(m1), sigma= fitted(m1,"sigma"),
                   nu =  fitted(m1, "nu"), bd=bd)
    } else
    {
      ures <- call("rqres", pfun=pfun,  type=dtype, 
                   ymin=fname$rqres[[1]][["ymin"]], 
                   y=y, mu= fitted(m1), sigma= fitted(m1,"sigma"),
                   nu =  fitted(m1, "nu"))
    } 
  }
  else 
  {
    if (fname$family[[1]] %in% .gamlss.bi.list)
    {
      ures <- call("rqres", pfun=pfun,  type=dtype, 
                   ymin=fname$rqres[[1]][["ymin"]], 
                   y=y, mu= fitted(m1), sigma= fitted(m1,"sigma"),
                   nu =  fitted(m1, "nu"), tau = fitted(m1, "tau"), bd=bd)
    } else
    {
      ures <- call("rqres", pfun=pfun,  type=dtype, 
                   ymin=fname$rqres[[1]][["ymin"]], 
                   y=y, mu= fitted(m1), sigma= fitted(m1,"sigma"),
                   nu =  fitted(m1, "nu"), tau = fitted(m1, "tau"))
    } 
  }
  # get output   
          m1$yRobust <- m1$y
                m1$y <- y       
  m1$residualsRobust <- resid(m1)
        m1$residuals <- eval(ures)
     m1$obs.weighted <- which(wei==0)
      m1$obs.bounded <- which(y1!=y)
            m1$bound <- bound
         m1$CD.bound <- CD.bound
  m1
}  
#---------------------------------------------------------

