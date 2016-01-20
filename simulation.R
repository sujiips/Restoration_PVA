require(hierfstat)

rm(list=ls())
Dir <- "/home/sujiips/teste"
setwd(Dir)

source('function_novo.R')
parms = list(
  Start <- 1
  , rep  <- 43
  , nloci <- 10
  , time <- 500
  , xDim <- 1667
  , yDim <- 30
  , dPollen <- 100
  , dSeed <- 100
  , avSeed <- c(50,70)
  , maxAge <- c(20,30,40)
  , adultAge <- 10
  , maxFathers <- 20
  , germ <- c(0.5, 0.6)
  , selection <- c(0.4,0.6)
  , output_years <- c(10, 30, 50, 100, 150, 200, 250, 300, 400, 500)
  )

begin(parms)



