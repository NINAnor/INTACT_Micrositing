# 
# kxv.R
#
# k-fold cross-validation for present/available probability models
#
# Copyright (c) 2005 by John Brzustowski (jbrzusto@fastmail.fm)
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# A copy of the GNU General Public License is available via WWW at
# http://www.gnu.org/copyleft/gpl.html.  You can also obtain it by
# writing to the Free Software Foundation, Inc., 59 Temple Place,
# Suite 330, Boston, MA  02111-1307  USA.
#
# k-fold cross-validation for present/available probability
# models. Partitions the "present" records of a dataset D into k
# subsets, then for each such subset S, fits a model to D-S (i.e. all
# data not in S) and tests the model against S.  The models estimate
# the relative probability of presence given input variable values.

# kxv: partition a dataset D.  Fit a model to each reduced dataset
# D-S, and then evaluate it on S, where S is a subset in the
# partition.  Returns a list of pairs, one pair per partition.  Each
# pair has an element "model", which is the result of fitting the
# model to D-S, and an element "pred" which is the set of 
# values predicted by that model for the data in S.
# Output from this function will
# typically be fed into something like binrank() below.
#
# There are helper functions, such as kxvglm(), for
# fixed, mixed, and discrete-choice models.
#
# To have details of each model printed as it is fit, 
# set the following flag to TRUE

kxvPrintFlag <- FALSE

# To turn off plotting, set the following flag to FALSE

kxvPlotFlag <- FALSE

# To allow for fuzzy binning, set the following
# to a small positive number (e.g. 0.001)

kxvFuzzFactor <- 0

# To turn off the load-time help message, comment out
# the following call to cat()

cat("Module kxv.R loaded.\
\
Type kxvtest() to run the following example of a fixed-effects\
model with sample Elk data:\
\
elk <- read.csv(file=\"HR_Summer.csv\", header=T)\
kxvglm( PTTYPE~DEM30M+DEM30M^2+VARCES1+FORGSUM+FORGSUM^2, elk, k=5,nbin=10)\
\
The test will do kxvPrintFlag <- TRUE in order to\
print a summary of each fitted model, but this\
flag is FALSE by default.\
\
To turn off plotting, do kxvPlotFlag <- FALSE.\
")

kxv <-  function(
                 
                 ## The presence/available indicator.
                 ## This is a vector with the same number of
                 ## elements as there are observations in
                 ## the dataframe.
                 ## TRUE indicates "present" observations,
                 ## and FALSE indicates "available" observations.
                 
                 present,

                 ## the number of subsets into which to divide the
                 ## data (not used if partition != NULL)
                 
                 k = 5,

                 ## The name of a factor by which to partion the data.
                 ## For each level of the factor, data with that level
                 ## are used to test a model built with the rest of
                 ## the dataset.  If partition is not given, the data
                 ## are randomly partitioned into k equal-sized (or as
                 ## close as possible) subsets.  If given, this must
                 ## be a member of names(data).
                 
                 partition = NULL,
                 
                 ## The dataframe

                 data,

                 ## The modelling function to call.
                 ## It must accept a single parameter, namely
                 ## the datset.

                 modfun,

                 ## The prediction function to call.
                 ## It must accept two parameters, namely
                 ## the output from modfun, and the
                 ## excluded dataset.

                 predfun
                 
                 )
{

  # In case the partitioning variable isn't already a factor
  # make it one, if it is supplied.

  if(! is.null(partition)) {
    data[,partition] <- as.factor(data[,partition])
  }

  # split the dataset into "present" and "available" parts
  
  pres <- data[present, ]
  avail <- data[!present, ]
  
  if(is.null(partition)) {
    
  # No partitioning factor has been specified, so partition the
  # "present" data by adding a randomly-generated factor.  If k
  # divides evenly into the number of samples, then each subset will
  # have the same size.  Otherwise, the sizes will differ by at most
  # one.  Add this factor, but empty, to the available data so
  # that the two sets can be glued back together below.

    partition <- "__kxv_bin"
    pres[,partition] <- as.factor(1+order(runif(dim(pres)[1])) %% k)
    avail[,partition] <- NA
    
  } 
  # a list for holding results
  
  res <- list()

  # We now have a factor used to partition the dataset,
  # and build models excluding one subset at a time of
  # "present" data.
  # Because all levels of the partition factor
  # might not occur in just the "present" data,
  # determine exactly which levels do occur there
  # in order to prevent empty models.

  used.levels <- levels(as.factor(as.character(pres[,partition])))
  
  for (subs in used.levels) {
    
    # build a model using the "available" data and the "present" data
    # from all but this subset
    
    model <- modfun(rbind(avail, pres[pres[,partition] != subs,]))

    # use the model to predict probabilities for the excluded subset
    
    pred <- predfun(model, pres[pres[,partition] == subs,])

    # record this output

    res[[as.character(subs)]] <- list(model=model, pred=pred)

    if(kxvPrintFlag) {
      print(summary(model))
    }
  }
  return(res)
}

respvar <- function(formula) {
  # return the name of the response variable from a formula
  all.vars(formula)[[attr(terms(formula), "response")]]
}

predvars <- function(formula){
  # return the names of the predictor variables from a formula
  all.vars(formula)[-attr(terms(formula),"response")]
}

clumps <- function(formula, data, n=10) {
  # look for the largest n clumps in a model
  # i.e. find the most common
  # ordered lists given by c(predictors, response)
  # for a given formula and dataset

  # In principle, this could be done using xtabs(formula, data),
  # but that generates fully-factorial output which is likely to
  # be larger than the dataset itself.
  
  vars <- c(respvar(formula), predvars(formula))

  # paste together string representations of response and predictors
  # for each observation, then turn this into a factor (which uses
  # fast internal hashing) and summarize the n largest (by count)
  # levels
  
  counts <- sort(summary(
                         as.factor(
                                   apply(
                                         data[, vars],
                                         1,
                                         function(x) paste(x, collapse=",")
                                         )
                                   ),
                         maxsum=n+1
                         )[1:n],
                 decreasing=TRUE)
  
  # split the string representations back into
  # their components
  
  clumps <- data.frame(cbind(matrix(nrow=min(n, length(counts)),
                                    c(lapply(
                                             names(counts),
                                             function(x)
                                             eval(
                                                  parse(text=paste("c(\"", gsub(",","\",\"",x), "\")",  sep=""))
                                                  )
                                             ),
                                      recursive=TRUE
                                      ),
                                    byrow=TRUE
                                    ),
                             counts
                             ),
                       row.names=NULL
                       )

  # fix column and row names
  
  colnames(clumps) <- c(vars, "(count)")
#  rownames(clumps) <- NULL

  # make numeric vectors back into such
  for (i in vars)
    if(is.numeric(data[, i]))
      clumps[, i] <- as.numeric(as.vector(clumps[, i]))

  return(clumps)
}
  
  

binrank <- function (
                     ## The output from the kxv function
                     out,
                     
                     ## The number of bins into which predicted
                     ## probabilities should be grouped.
                     
                     nbin = 10,
                     
                     ## A function returning fitted model
                     ## values from the model output.
                     
                     fitted,
                     
                     ## names of the excluded subsets
                     ## (NULL corresponds to the usual
                     ## case of k random subsets)
                     
                     subnames = NULL,
                     
                     ## A fuzz factor for breaking up
                     ## clumps in the data.  Predicted and
                     ## fitted model values will have normally-
                     ## distributed error with this relative
                     ## standard deviation added to them.
                     
                     fuzzfactor = 0
                     
                     )
{

  # initialize empty lists for summarizing test results
  
  y <- c()
  tab <- c()
  
  for (subs in names(out)) {

    # unpack the results for this subset
    
    model <- out[[subs]]$model
    pred <- out[[subs]]$pred
    avail <- fitted(model)

    # add fuzz to allow clean binning; this is
    # normally-distributed random error with SD
    # equal to fuzzfactor times the mean
    # magnitude among the fitted and predicted
    # model values

    fuzz <- fuzzfactor * mean(abs(c(pred, avail)))
    pred <- pred + rnorm(length(pred), 0, fuzz)
    avail <- avail + rnorm(length(avail), 0, fuzz)
    
    # find the (right) boundaries for bins holding equal
    # (or as near as possible) numbers of RSF values
    # from only the "available" observations

    binbound <- sort(avail)[ceiling((1:nbin)/nbin*length(avail))]

    # if unable to find distinct bin boundaries, data are
    # too clumped for the specified number of bins.  Fail
    # with an error message.  This should almost never
    # happen with a non-zero fuzz-factor.

    if (! all(diff(binbound)))
      stop(paste("\nUnable to create", nbin, "distinct bins.\n"),
           "Your model and data yield so many identical fitted values\n",
           "that some adjacent bins would have the same endpoint.\n",
           "Try a smaller number for nbin or\n",
           "specify a fuzz factor (e.g. kxvFuzzFactor <- 0.001).\n",
           "You can use clumps(formula, data) to see a possible cause of this problem.\n"
           )
    
    # For the subset of "present" values, assign each RSF
    # value to its appropriate bin.  Bins include
    # their right but not their left endpoints.
    # The leftmost bin is extended to minus infinity.
    # The approx function here generates a piecewise-constant
    # function (i.e. bin number as a function of RSF value)

    binnum <- approx(x=binbound, y=1:nbin, xout=pred,
                     method="constant", f=1, rule=2)$y

    # compute the counts per bin
    # (a bit ugly since the "summary" function
    # leaves out levels with counts of zero)
    counts <- rep(0, nbin)
    prescounts <- summary(as.factor(binnum), maxsum=nbin)
    for (i in names(prescounts))
      counts[[as.numeric(i)]] <- prescounts[[i]]

    
    # "area-adjust" counts:  we want 1 to mean what would
    # be expected if "present" points were selected at random.
    # Since each bin is (nearly) equal in size, we just
    # scale the counts uniformly so that 1 has the desired
    # meaning; i.e. for each bin, we want fraction of "present"
    # in that bin points divided by the 
    # fraction of all points in that bin
    
    counts <- (counts / sum(counts)) / (1 / nbin)
    
    # record the y-coordinates for graphing

    y <- cbind(y, counts)
    
    # test the association between area-adjusted counts and
    # bin rank
    
    test <- cor.test(counts, 1:nbin, method="spearman")

    # record test results for printing

    tab <- rbind(tab, c(r=test$estimate, p=test$p.value))
  }

  countsMean <- apply(y, 1, mean)
  countsSD <- apply(y, 1, sd)
  
  kxvPlotFlag <- FALSE
  
  if(kxvPlotFlag) {
    # do both plots in one page
    par(mfcol=c(2,1))
    
    # plot the curves for each partition
    matplot(1:nbin, y, xlab="Binned RSF Score",
            ylab="Area-adjusted Frequency",
            type="b", pch=1:nbin,
            ylim=c(0,max(c(y))))

    # if subset names are given, show a legend on the graph
    if(!is.null(subnames))
      legend (x="topleft", legend=subnames, col=1:length(subnames),
              pch=1:length(subnames), bty="n", cex=0.5)

    # plot the curve for the mean, with "arrows"
    # showing +/- SD
    
    plot(1:nbin, countsMean, xlab="Binned RSF Score",
         ylab=c("Area-adjusted Frequency"," (mean+/-SD)"),
         type="b",
         ylim=c(0, max(countsMean+countsSD)))
    arrows(1:nbin, countsMean - countsSD, 1:nbin, countsMean + countsSD,
            code = 3, angle = 90, length=0.05)

  }
  
  # do the test for the mean area-adjusted counts

  test <- cor.test(countsMean, 1:nbin, method="spearman")
  tab <- rbind(tab, c(r=test$estimate, p=test$p.value))
  rownames(tab) <- c(names(out), "(meanfreq)")
  return(tab)
}  
  
    

# kxvglm: a helper function to call kxv with sensible defaults for a
# fixed-effects model using glm, and then feed the output into binrank

kxvglm <- function(
                   ## formula for the model
                   
                   formula,

                   ## the dataset to use
                   
                   data,

                   ## the number of subsets into which to divide the
                   ## data (not used if partition != NULL)
                   
                   k = 5,
                   
                   ## the number of bins to use in binrank

                   nbin = 10,
                   
                   ## the value of the response variable indicating
                   ## a "present" (vs. "available") observation

                   presval = 1,
                   
                   ## The name of a factor by which to partion the
                   ## data.  For each level of the factor, data with
                   ## that level are used to test a model built with
                   ## the rest of the dataset.  If not given, the data
                   ## are randomly partitioned into k equal-sized (or
                   ## as close as possible) subsets
                   
                   partition = NULL,

                   ## the family of link/error function to use

                   family = binomial(link="logit"),
                   
                   ## additional arguments to glm
                   
                   ...
                   
                   )
{
  presvar <- respvar(formula)
  # recode the presence variable so that
  # T is present, F is available

  present <- data[, presvar] <- data[, presvar] == presval
  binrank(   kxv( present=present,
                 k=k,
                 partition=partition,
                 data=data,
                 modfun=function(dat) {
                   glm(formula=formula, family=family,
                       data=dat, ...)
                 },
                 predfun=function(mod, dat) {
                   predict.glm(object=mod, newdata=dat, type="response",allow.new.levels=T)
                 }
                 ),
          nbin=nbin,
          fitted=function(mod){
            mod$fitted.values[1:sum(! present)]
          },
          subnames=if(is.null(partition)) NULL else
  paste(partition,levels(as.factor(data[,partition])),sep="="),
          fuzzfactor=kxvFuzzFactor
          
          )
}


kxvtest <- function() {
  # run with sample dataset
  # the print command is only needed so that when "source"ing
  # this file, the results of running the model are printed.

  elk <- read.csv(file="HR_Summer.csv", header=T)
  print("Read file HR_Summer.csv ")
  kxvPrintFlag <<- TRUE
  print(kxvglm( PTTYPE~DEM30M+DEM30M^2+VARCES1+FORGSUM+FORGSUM^2, elk, k=5, nbin=10))
  kxvPrintFlag <<- FALSE

  # The other two are implemented, but I don't know enough about
  # these types of models to give them a proper test.
  # In particular, there are likely problems with exactly what is
  # being used to get fitted and predicted values from the model structure.
  #
  # kxvglmmPQL( fixed=PTTYPE~DEM30M+DEM30M^2+FORGSUM+FORGSUM^2, random=~1|ELKNUM, elk)
  #
  # kxvclog(PTTYPE~DEM30M+DEM30M^2+VARCES1+FORGSUM+FORGSUM^2, elk)

}

## The following helper functions will almost certainly not work as
## written.
## They are provided as templates which might help you implement kxv
## for other types of models.

## kxvglmmPQL: a helper function to call kxv with sensible defaults for a
## mixed-effects model using glmmPQL, and then feed the output into binrank

 kxvlmer <- function(
                    ## formula for the model
                   
                    formula,
                   
                    ## the dataset to use
                   
                    data,

                    ## the number of subsets into which to divide the
                    ## data (not used if partition != NULL)
                   
                    k = 5,
                   
                    ## the number of bins to use in binrank

                    nbin = 10,
                   
                    ## the value of the response variable indicating
                    ## a "present" (vs. "available") observation

                    presval = 1,
                       
                    ## The name of a factor by which to partion the
                    ## data.  For each level of the factor, data with
                    ## that level are used to test a model built with
                    ## the rest of the dataset.  If not given, the data
                    ## are randomly partitioned into k equal-sized (or
                    ## as close as possible) subsets
                   
                    partition = NULL,

                    ## the family of link/error function to use

                    family = poisson,
                   
                    ## additional arguments to lmer
                   
                    ...
                   
                    )
 {
   library("lme4")
   presvar <- respvar(formula)
   present <- data[, presvar] <- data[, presvar] == presval
   binrank( kxv(present=present, 
                k=k,
                partition=partition,
                data=data,
                modfun=function(dat) {
                  lmer(formula=formula, family=family,
                      data=dat, ...)
                },
                predfun=function(mod, dat) {
                  predict(object=mod, newdata=dat,allow.new.levels=T)
                }
                ),
           nbin=nbin,
           fitted=function(mod){
             mod$fitted[1:sum(! present)]

           }
           )
 }

 kxvglmer <- function(
   ## formula for the model
   
   formula,
   
   ## the dataset to use
   
   data,
   
   ## the number of subsets into which to divide the
   ## data (not used if partition != NULL)
   
   k = 5,
   
   ## the number of bins to use in binrank
   
   nbin = 10,
   
   ## the value of the response variable indicating
   ## a "present" (vs. "available") observation
   
   presval = 1,
   
   ## The name of a factor by which to partion the
   ## data.  For each level of the factor, data with
   ## that level are used to test a model built with
   ## the rest of the dataset.  If not given, the data
   ## are randomly partitioned into k equal-sized (or
   ## as close as possible) subsets
   
   partition = NULL,
   
   ## the family of link/error function to use
   
   family = binomial,
   
   ## additional arguments to lmer
   
   ...
   
 )
 {
   library("lme4")
   presvar <- respvar(formula)
   present <- data[, presvar] <- data[, presvar] == presval
   binrank( kxv(present=present, 
                k=k,
                partition=partition,
                data=data,
                modfun=function(dat) {
                  glmer(formula=formula, family=family,
                       data=dat, ...)
                },
                predfun=function(mod, dat) {
                  predict(object=mod, newdata=dat,type="response",allow.new.levels=T)
                }
   ),
   nbin=nbin,
   fitted=function(mod){
     fitted(mod)[1:sum(! present)]
     
   }
   )
 }
 
 
## # kxvclog: a helper function to call kxv with sensible defaults for a
## # discrete choice model using clogit, and then feed the output into binrank

 kxvclog <- function(
                    ## formula for the model
                   
                    formula,

                    ## the dataset to use
                   
                    data,

                    ## the number of subsets into which to divide the
                    ## data (not used if partition != NULL)
                  
                    k = 5,
                   
                    ## the number of bins to use in binrank

                    nbin = 10,
                   
                    ## the value of the response variable indicating
                    ## a "present" (vs. "available") observation

                    presval = 1,
                    
                    ## The name of a factor by which to partion the
                    ## data.  For each level of the factor, data with
                    ## that level are used to test a model built with
                    ## the rest of the dataset.  If not given, the data
                    ## are randomly partitioned into k equal-sized (or
                    ## as close as possible) subsets
                   
                    partition = NULL,

                    ## additional arguments to clogit
                   
                    ...
                   
                    )
 {
   library("survival")
#   library("Design")
   presvar <- respvar(formula)
   present <- data[, presvar] <- data[, presvar] == presval
   binrank( kxv( present=present,
                k=k,
                partition=partition,
                data=data,
                modfun=function(dat) {
                  clogit(formula=formula, data=dat, ...)
                },
                predfun=function(mod, dat) {
                  predict(object=mod, newdata=dat, type="risk",allow.new.levels=T)
                }
                ),
           nbin=nbin,
           fitted=function(mod){
             predict(object=mod, type="risk",allow.new.levels=T)[1:sum(! present)]
           }
           )
 }



