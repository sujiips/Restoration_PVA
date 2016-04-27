begin  <-  function( parms ){  

#------------------------------------------------------#    
# fsRate
#   function for age dependent survival probability
#------------------------------------------------------#
    fsRate <- function(age){
    ifelse(age >= 1 & age <= 2, 0.29, 
           ifelse(age < 6, 0.6,
                ifelse(age < 10, 0.7, 
                    ifelse(age < 20, 0.75,
                        ifelse(age < 40, 0.85,
                  ifelse(age >= maxAge[m], 0, 0.95))))))
  }
  
  
  #--------------------------------------------------------------------------------------#
  # f_reproduction
  #    function for reproduction
  # notes: This function sets the coordinates (x,y) for new individuals; defines their
  #        alleles from the parents' genotypes; and the germination probabilty from
  #        the inividual heterozygosity
  #---------------------------------------------------------------------------------------#
  
  f_reproduction <- function( ind1 = NULL, ind2 = NULL, seedInd){
    # defining new individual's coordinates (x, y)
    theta <- runif(seedInd*fathers, 0, 2*pi)
    disp <- 1:dSeed
    probDisp <- (0.6942213/81.4806772)*((disp/81.4806772)^(0.6942213-1))*exp(-(disp/81.4806772)^0.6942213)
    x <- as.numeric(dat[ind1,1]) + round(cos(theta)*sample(disp, seedInd*fathers, prob=probDisp, replace = T))
    y <- as.numeric(dat[ind1,2]) + round(sin(theta)*sample(disp, seedInd*fathers, prob=probDisp, replace = T))
    
    # Sampling alleles from parents genotypes for new individuals 
    genotypes <- matrix(rep(0,seedInd*fathers*nloci*2),ncol=nloci*2)
    for(locus in 1:nloci){
      genotypes[,c(2*locus-1, 2*locus)] <- matrix(unlist(
        c(sample(dat[ind1,(2*locus+1):(2*locus+2)], seedInd*fathers, replace=T), 
          as.vector(t(sample(dat[ind2,(2*locus+1):(2*locus+2)], seedInd, replace=T)))
        )), byrow=F)
    }
    genotypes <- data.frame(genotypes)
    colnames(genotypes) <- paste0("l",rep(1:nloci,each=2),"_",1:2)
    age <- rep(0, seedInd*fathers)
    
    # Defining probability of germination as function of heterozygosty  
    a1 <- genotypes[,seq( 1,(nloci*2),2 )]
    a2 <- genotypes[,seq(2,(nloci*2),2 )]
    homo <- a1==a2
    ho <- apply(homo, 1, sum)
    Ho <- 1-(ho/nloci)
    Fis <- 1-(Ho/He)
    sRate <- ifelse(Fis <= 0, germ[g], germ[g]*(1-selection[s]*Fis)) 
    
    # Kill individuals whose coordinates are the same as adult individuals
    parent1 <- ind1
    parent2 <- rep(ind2, each=seedInd)
    old <- c(ind1, ind2)
    x1 <- adult$x
    y1 <- adult$y
    f_sameSpot <- function(o){
      which( x1[o] == x & y1[o] == y )}
    sameSpot <- unique(unlist(lapply(old, f_sameSpot)))
    age[sameSpot] <- NA
    
    # Define new individual data  
    new_ind <<- cbind(x, y, genotypes, age, parent1, parent2, sRate)
  }
  
  for(g in 1:length(germ)){
    for(s in 1: length(selection)){ 
      for(m in 1:length(maxAge)){
        for(av in 1:length(avSeed)){  
          
          ## Start simulations
          for (r in Start:rep){
            setwd(Dir)
            dat <- read.table( 'datain.txt', header = TRUE, sep = '\t')
            
            sRate <- fsRate(dat$age)  
            dat <- cbind(dat, sRate)
            
            new_dir <- paste0(Dir,"/out")
            dir.create(new_dir)
            setwd(new_dir)
            
            for(t in 1:time){
              if( time%%10 ==0){print(t)}
              
              ## Define new sRate 
              dat$sRate <- fsRate(dat$age)  
              
              ### Reproduction
              ## Calculate expected heterozygosity under Hardy-Weinberg Equilibrium
              He_loci <- rep(0, nloci)
              for(locus in 1:nloci){
                genotypes <- c(dat[,(2*locus+1)], dat[,(2*locus+2)])
                He_loci[locus] <- 1-sum((table(genotypes)/length(genotypes))^2)
              }
              He <- mean(He_loci)
              
              
              ## Find reproductive pairs and produce seeds
              adult <- subset(dat, dat$age >= adultAge)
              
              if(length(adult[,1]) > 0){
                for(parent in 1:length(adult[,1])){
                  pairs <- expand.grid(parent, 1:length(adult[,1]))
                  neigh <- which(apply(pairs, 1, function(p){
                    (adult$x[p[1]] - adult$x[p[2]])^2 + (adult$y[p[1]]-adult$y[p[2]])^2 <= dPollen^2}))
                  
                  if(length(neigh) < maxFathers){
                    seedInd <- round(avSeed[av]/length(neigh))
                    sexPair <- neigh
                    fathers <- length(sexPair)
                  } else {
                    seedInd <- round(avSeed[av]/maxFathers)
                    sexPair <- sample(neigh, maxFathers, replace=T)
                    fathers <- maxFathers
                  }
                  f_reproduction(parent, sexPair, seedInd)
                  dat <- rbind(dat, new_ind) 
                }
              }
              
              ### Checking if the individuals will survive to the next time
              rand <- runif(length(dat[,1]), 0, 1)
              sRate_death <- which( rand > dat$sRate )
              dat$age[sRate_death] <- NA
              
              out <- which(dat$x < 1 | dat$x > xDim | dat$y < 1 | dat$y > yDim)
              dat$age[out] <- NA
              
              dat <- dat[!is.na(dat$age),]
              
              ## Individuals that survived get 1 year older
              dat$age <- dat$age + 1
              
              
              ### Save output
              if(t %in% output_years){
                mypath <- file.path(getwd(), paste("dataout_", t, ".txt", sep=""))
                write.table( dat, file=mypath, row.names= FALSE, sep = '\t' )
              }
            }
            
## Estimate genetic parameters
  setwd(new_dir)

# List output files  
  filename1 <- list.files(new_dir, full.names=F)
  if(length(filename1)==0){break}
      else if(length(filename1)!=0){
      file.n <- (length(filename1))
      filename <- filename1[1:file.n]
      n1 <- gsub("dataout_", "", filename)
      n <- as.integer(gsub(".txt", "", n1)) 
      npop <- (1:length(filename))

# Convert outputfiles to fstat format, as defined in the hierfstat documentation                    
  fstat_dir <- paste0(new_dir,"/fstat")
  dir.create(fstat_dir)
  
  for(f in 1:length(filename)){
    setwd(new_dir)
    DataTrees <- read.table(filename[f], header=TRUE, sep="\t")
    DataTrees <- subset(DataTrees, age > 1)
    DataTrees <- DataTrees[,c(3:(nloci*2+2))]
                
  for (i in seq(1,nloci*2,2)){
    DataTrees[,i]<-as.character(paste(DataTrees[,i],DataTrees[,i+1],sep= ""))}
    DataTrees2<-DataTrees[,-seq(2,20,2)]   
    colnames(DataTrees2)<- gsub("_1","",colnames(DataTrees2));  
    pop <- rep(n[f], length(DataTrees2[,1]))
    Data_fstat <- cbind(pop, DataTrees2)
                
    setwd(fstat_dir)
    write.table(Data_fstat, file = filename[f], quote=F, row.names=F)
  }
              
# Merge all output files in one
  all <- function(id) {
    files <- list.files(fstat_dir, full.names=F)
    all_data <-data.frame()
    for (i in 1:length(files)){
      all_data <- rbind(all_data, read.table(files[i], head=T))
    }
    alldata <- all_data[order(all_data[,1]),]
  }
              
  data <- all(npop)
  pop <- length(npop)
  pop.names <- unique(data$pop)
  alleles <- allele.count(data, diploid=T)
              
  # Allelic richness per locus
   allelic.r <- allelic.richness(data, min.n=NULL, diploid=T)
   mean.rich <- round(apply(allelic.r$Ar, 2, mean), digits=1)
              
  # Mean number of alleles per locus
    a <- nb.alleles(data)
    mean.a <- round(apply(a, 2, mean), digits=1)

  ## Descriptive genetic statistics                
    basic <- basic.stats1(data, diploid=T, digits=3)
              
    # number of individuals per locus
    n.ind <- basic$n.ind.samp
    mean.ind <- round(apply(n.ind, 2, mean), digits=1) 
              
    # Fis per pop
    fis <- data.frame(basic$Fis)
    mean.fis <- round(apply(fis, 2, mean), digits=3)
              
    ## Ho per pop
    Ho <- basic$Ho
    mean.Ho <- round(apply(Ho, 2, mean), digits=3)
    sd.Ho <- round(apply(Ho, 2, sd), digits=3)
              
   ## Hs per pop
   Hs <- basic$Hs
   mean.Hs <- round(apply(Hs, 2, mean), digits=3)
   sd.Hs <- round(apply(Ho, 2, sd), digits=3)
              
              
 ## Create summary table
    basic.table <- cbind(pop.names, mean.ind, mean.a, mean.rich, mean.Ho, sd.Ho, mean.Hs, sd.Hs, mean.fis)
    colnames(basic.table) <- c("year", "ind", "A", "R", "Ho", "Ho(SD)", "Hs", "Hs(SD)", "fis")  
    write.table(basic.table, file=paste0(fstat_dir,"/basic.txt"), quote=F, row.names=F, sep="\t") 
}

# Change the name of the folder with output files 
  setwd(Dir)
  N <- paste(g,s,m,av, sep="")
  name <- paste("s", N, r, sep="_")
  file.rename("out", name) 
}
}
}
}
}
}


#---------------------------------------------------------------------------#
# basic.stats1
#   function for basic genetic descriptive statistics
#
# notes: this function was modified from package hierfstat (Goudet 2014).  
#        The modification made sets Fis to 1, when there is only one allele 
#        in a locus.
#---------------------------------------------------------------------------#

basic.stats1 <- function (data, diploid = TRUE, digits = 4) {
  loc.names <- names(data)[-1]
  if (length(table(data[, 1])) < 2) 
    data[dim(data)[1] + 1, 1] <- 2
  p <- pop.freq(data, diploid)
  n <- t(ind.count(data))
  if (diploid) {
    dum <- getal.b(data[, -1])
    Ho <- dum[, , 1] == dum[, , 2]
    sHo <- (1 - t(apply(Ho, 2, fun <- function(x) tapply(x, 
                                                         data[, 1], mean, na.rm = TRUE))))
    mHo <- apply(sHo, 1, mean, na.rm = TRUE)
  }
  else {
    sHo <- NA
    mHo <- NA
  }
  sp2 <- lapply(p, fun <- function(x) apply(x, 2, fun2 <- function(x) sum(x^2)))
  sp2 <- matrix(unlist(sp2), nrow = dim(data[, -1])[2], byrow = TRUE)
  if (diploid) {
    Hs <- (1 - sp2 - sHo/2/n)
    Hs <- n/(n - 1) * Hs
    Fis  <- ifelse(Hs==0, Fis <- 1, Fis <- 1 - sHo/Hs)
  }
  else {
    Hs <- n/(n - 1) * (1 - sp2)
    Fis <- NA
  }
  np <- apply(n, 1, fun <- function(x) sum(!is.na(x)))
  mn <- apply(n, 1, fun <- function(x) {
    np <- sum(!is.na(x))
    np/sum(1/x[!is.na(x)])
  })
  msp2 <- apply(sp2, 1, mean, na.rm = TRUE)
  mp <- lapply(p, fun <- function(x) apply(x, 1, mean, na.rm = TRUE))
  mp2 <- unlist(lapply(mp, fun1 <- function(x) sum(x^2)))
  if (diploid) {
    mHs <- mn/(mn - 1) * (1 - msp2 - mHo/2/mn)
    Ht <- 1 - mp2 + mHs/mn/np - mHo/2/mn/np
    mFis <- ifelse(Ht==0, mFis <- 1, mFis  <- 1 - mHo/mHs)
  }
  else {
    mHs <- mn/(mn - 1) * (1 - msp2)
    Ht <- 1 - mp2 + mHs/mn/np
    mFis <- NA
  }
  Dst <- Ht - mHs
  Dstp <- np/(np - 1) * Dst
  Htp = mHs + Dstp
  Fst = Dst/Ht
  Fstp = Dstp/Htp
  Dest <- Dstp/(1 - mHs)
  res <- data.frame(cbind(mHo, mHs, Ht, Dst, Htp, Dstp, Fst, 
                          Fstp, mFis, Dest))
  names(res) <- c("Ho", "Hs", "Ht", "Dst", "Htp", "Dstp", "Fst", 
                  "Fstp", "Fis", "Dest")
  if (diploid) {
    rownames(sHo) <- loc.names
    rownames(Fis) <- loc.names
  }
  overall <- apply(res, 2, mean, na.rm = TRUE)
  overall[7] <- overall[4]/overall[3]
  overall[8] <- overall[6]/overall[5]
  overall[9] <- 1 - overall[1]/overall[2]
  overall[10] <- overall[6]/(1 - overall[2])
  names(overall) <- names(res)
  all.res <- list(n.ind.samp = n, pop.freq = lapply(p, round, 
                                                    digits), Ho = round(sHo, digits), Hs = round(Hs, digits), 
                  Fis = round(Fis, digits), perloc = round(res, digits), 
                  overall = round(overall, digits))
  class(all.res) <- "bas.stats"
  all.res
}
