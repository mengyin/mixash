
################################## GENERIC FUNCTIONS ############################
# find matrix of densities at y, for each component of the mixture
# INPUT y is an n-vector
# OUTPUT k by n matrix of densities
compdens = function(x,y,log=FALSE){
  UseMethod("compdens")
}
compdens.default = function(x,y,log=FALSE){
  stop(paste("Invalid class", class(m), "for first argument in",  match.call()))  
}

#standard deviations
comp_sd = function(m){
  UseMethod("comp_sd")
}
comp_sd.default = function(m){
  stop("method comp_sd not written for this class")
}

#second moments
comp_mean2 = function(m){
  UseMethod("comp_mean2")
}
comp_mean2.default = function(m){
  comp_sd(m)^2 + comp_mean(m)^2
}


#return the overall mean of the mixture
calc_mixmean = function(m){
  UseMethod("calc_mixmean")
}
calc_mixmean.default = function(m){
  sum(m$pi * comp_mean(m))
}

#return the overall second moment of the mixture
mixmean2 = function(m){
  UseMethod("mixmean2")
}
mixmean2.default = function(m){
  sum(m$pi * comp_mean2(m))
}

#return the overall sd of the mixture
calc_mixsd = function(m){
  UseMethod("calc_mixsd")
}
calc_mixsd.default = function(m){
  sqrt(mixmean2(m)-calc_mixmean(m)^2)
}

#means
comp_mean = function(m){
  UseMethod("comp_mean")
}
comp_mean.default = function(m){
  stop("method comp_mean not written for this class")
}

#number of components
ncomp = function(m){
  UseMethod("ncomp")
}
ncomp.default = function(m){
  return(length(m$pi))
}

#return mixture proportions, a generic function
mixprop = function(m){
  UseMethod("mixprop")
}
mixprop.default = function(m){
  m$pi
}

#' @title mixcdf
#'
#' @description Returns cdf for a mixture (generic function)
#' 
#' @details None
#' 
#' @param x a mixture (eg of type normalmix or unimix)
#' @param y locations at which cdf to be computed
#' @param lower.tail: boolean indicating whether to report lower tail
#' 
#' @return an object of class normalmix
#' 
#' @examples mixcdf(normalmix(c(0.5,0.5),c(0,0),c(1,2)),seq(-4,4,length=100))
#' 
mixcdf = function(x,y,lower.tail=TRUE){
  UseMethod("mixcdf")
}

mixcdf.default = function(x,y,lower.tail=TRUE){
  x$pi %*% comp_cdf(x,y,lower.tail)
}

#find cdf for each component, a generic function
comp_cdf = function(x,y,lower.tail=TRUE){
  UseMethod("comp_cdf")
}
comp_cdf.default = function(x,y,lower.tail=TRUE){
  stop("comp_cdf not implemented for this class")
}


#find density at y, a generic function
dens = function(x,y){
  UseMethod("dens")
}
dens.default = function(x,y){
  return (x$pi %*% compdens(x, y))
}

#find log likelihood of data in x (a vector) for mixture in m
loglik = function(m,x){
  UseMethod("loglik")
}
loglik.default = function(m,x){
  sum(log(dens(m,x)))
}

#find log likelihood of data in betahat, when 
#the mixture m is convolved with a normal with sd betahatsd
#betahatsd is an n vector
#betahat is an n vector
#v is the degree of freedom
loglik_conv = function(m,betahat,betahatsd,v,FUN="+"){
  UseMethod("loglik_conv")
}

loglik_conv.default = function(m,betahat,betahatsd,v,FUN="+"){
  sum(log(dens_conv(m,betahat,betahatsd,v,FUN)))
}

#compute the density of the components of the mixture m
#when convoluted with a normal with standard deviation s
#or a scaled (se) student.t with df v
#the density is evaluated at x
#x and s are n-vectors
#m is a mixture with k components
#output is a k by n matrix of densities
compdens_conv = function(m,x,s,v,FUN="+"){
  UseMethod("compdens_conv")
}
compdens_conv.default = function(m,x,s,v,FUN="+"){
  stop(paste("Invalid class", class(m), "for first argument in",  match.call()))  
}

#compute density of mixture m convoluted with normal of sd (s) or student t with df v
#at locations x
#m is a mixture
#x is an n vector
#s is an n vector or integer
dens_conv = function(m,x,s,v,FUN="+"){
  UseMethod("dens_conv")
}
dens_conv.default = function(m,x,s,v,FUN="+"){
  colSums(m$pi * compdens_conv(m,x,s,v,FUN))
}

#compute the posterior prob that each observation
#came from each component of the mixture m
#output a k by n vector of probabilities
#computed by weighting the component densities by pi
#and then normalizing
comppostprob=function(m,x,s,v){
  UseMethod("comppostprob") 
}

comppostprob.default = function(m,x,s,v){
  tmp= (t(m$pi * compdens_conv(m,x,s,v))/dens_conv(m,x,s,v))
  ismissing = (is.na(x) | is.na(s))
  tmp[ismissing,]=m$pi
  t(tmp)
}
# evaluate cdf of posterior distribution of beta at c
# m is the prior on beta, a mixture
# c is location of evaluation
# assumption is betahat | beta \sim N(beta,sebetahat)
# m is a mixture with k components
# c a scalar
# betahat, sebetahat are n vectors 
# output is a k by n matrix
compcdf_post=function(m,c,betahat,sebetahat,v){
  UseMethod("compcdf_post")
}
compcdf_post.default=function(m,c,betahat,sebetahat,v){
  stop("method compcdf_post not written for this class")
}


cdf_post = function(m,c,betahat,sebetahat,v){
  UseMethod("cdf_post")
}
cdf_post.default=function(m,c,betahat,sebetahat,v){
  colSums(comppostprob(m,betahat,sebetahat,v)*compcdf_post(m,c,betahat,sebetahat,v))
}

#output posterior mean for beta for prior mixture m,
#given observations betahat, sebetahat, df v
postmean = function(m, betahat,sebetahat,v){
  UseMethod("postmean")
}
postmean.default = function(m,betahat,sebetahat,v){
  colSums(comppostprob(m,betahat,sebetahat,v) * comp_postmean(m,betahat,sebetahat,v))
}



#output posterior mean-squared value for beta for prior mixture m,
#given observations betahat, sebetahat, df v
postmean2 = function(m, betahat,sebetahat,v){
  UseMethod("postmean2")
}
postmean2.default = function(m,betahat,sebetahat,v){
  colSums(comppostprob(m,betahat,sebetahat,v) * comp_postmean2(m,betahat,sebetahat,v))
}

#output posterior sd for beta for prior mixture m,
#given observations betahat, sebetahat, df v
postsd = function(m,betahat,sebetahat,v){
  UseMethod("postsd")
}
postsd.default = function(m,betahat,sebetahat,v){
  sqrt(postmean2(m,betahat,sebetahat,v)-postmean(m,betahat,sebetahat,v)^2)
}

#output posterior mean-squared value for beta for prior mixture m,
#given observations betahat, sebetahat, df v
comp_postmean2 = function(m,betahat,sebetahat,v){
  UseMethod("comp_postmean2")
}
comp_postmean2.default = function(m,betahat,sebetahat,v){
  comp_postsd(m,betahat,sebetahat,v)^2 + comp_postmean(m,betahat,sebetahat,v)^2
}


#output posterior mean for beta for each component of prior mixture m,
#given observations betahat, sebetahat, df v
comp_postmean = function(m, betahat,sebetahat,v){
  UseMethod("comp_postmean")
}
comp_postmean.default = function(m,betahat,sebetahat,v){
  stop("method comp_postmean not written for this class")
}




#output posterior sd for beta for each component of prior mixture m,
#given observations betahat, sebetahat, df v
comp_postsd = function(m, betahat,sebetahat,v){
  UseMethod("comp_postsd")
}
comp_postsd.default = function(m,betahat,sebetahat,v){
  stop("method comp_postsd not written for this class")
}

#find nice limits of mixture m for plotting
min_lim = function(m){
  UseMethod("min_lim")
}
min_lim.default=function(m){
  -5
}

max_lim = function(m){
  UseMethod("max_lim")
}
max_lim.default=function(m){
  5
}


#plot density of mixture
plot_dens = function(m,npoints=100,...){
  UseMethod("plot_dens")
}
plot_dens.default = function(m,npoints=100,...){
  x = seq(min_lim(m),max_lim(m),length=npoints)
  plot(x,dens(m,x),type="l",xlab="density",ylab="x",...)
}

plot_post_cdf = function(m,betahat,sebetahat,v,npoints=100,...){
  UseMethod("plot_post_cdf")
}
plot_post_cdf.default = function(m,betahat,sebetahat,v,npoints=100,...){
  x = seq(min_lim(m),max_lim(m),length=npoints)
  x_cdf = vapply(x,cdf_post,FUN.VALUE=betahat,m=m,betahat=betahat,sebetahat=sebetahat,v=v)
  plot(x,x_cdf,type="l",xlab="x",ylab="cdf",...)
  # for(i in 2:nrow(x_cdf)){
  #   lines(x,x_cdf[i,],col=i)
  # }
}

############################### METHODS FOR normalmix class ###########################

#' @title Constructor for normalmix class
#'
#' @description Creates an object of class normalmix (finite mixture of univariate normals)
#' 
#' @details None
#' 
#' @param pi vector of mixture proportions
#' @param mean vector of means
#' @param sd: vector of standard deviations
#' 
#' @return an object of class normalmix
#' 
#' @examples normalmix(c(0.5,0.5),c(0,0),c(1,2))
#' 
normalmix = function(pi,mean,sd){
  structure(data.frame(pi,mean,sd),class="normalmix")
}

comp_sd.normalmix = function(m){
  m$sd
}

comp_mean.normalmix = function(m){
  m$mean
}

compdens.normalmix = function(x,y,log=FALSE){
  k=ncomp(x)
  n=length(y)
  d = matrix(rep(y,rep(k,n)),nrow=k)
  return(matrix(dnorm(d, x$mean, x$sd, log),nrow=k))  
}

#density of convolution of each component of a normal mixture with N(0,s^2) or s*t(v) at x
# x an n-vector at which density is to be evaluated
#return a k by n matrix
#Note that convolution of two normals is normal, so it works that way
compdens_conv.normalmix = function(m,x,s,v,FUN="+"){
  if(!is.null(v)){
    stop("method comp_postsd of normal mixture not written for df!=NULL")
  }
  if(length(s)==1){s=rep(s,length(x))}
  sdmat = sqrt(outer(s^2,m$sd^2,FUN)) #n by k matrix of standard deviations of convolutions
  return(t(dnorm(outer(x,m$mean,FUN="-")/sdmat)/sdmat))
}


comp_cdf.normalmix = function(x,y,lower.tail=TRUE){
  vapply(y,pnorm,x$mean,x$mean,x$sd,lower.tail)
}

#c is a scalar
#m a mixture with k components
#betahat a vector of n observations
#sebetahat an n vector of standard errors
#return a k by n matrix of the posterior cdf
compcdf_post.normalmix=function(m,c,betahat,sebetahat,v){
  if(!is.null(v)){
    stop("Error: normal mixture for student-t likelihood is not yet implemented")
  }  
  k = length(m$pi)
  n=length(betahat)
  #compute posterior standard deviation (s1) and posterior mean (m1)
  s1 = sqrt(outer(sebetahat^2,m$sd^2,FUN="*")/outer(sebetahat^2,m$sd^2,FUN="+"))
  ismissing = (is.na(betahat) | is.na(sebetahat))
  s1[ismissing,]=m$sd
  
  m1 = t(comp_postmean(m,betahat,sebetahat,v))
  t(pnorm(c,mean=m1,sd=s1))
}

#return posterior mean for each component of prior m, given observations betahat and sebetahat
#input, m is a mixture with k components
#betahat, sebetahat are n vectors
#output is a k by n matrix
comp_postmean.normalmix = function(m,betahat,sebetahat,v){
  if(!is.null(v)){
    stop("method comp_postmean of normal mixture not written for df!=NULL")
  }
  tmp=(outer(sebetahat^2,m$mean, FUN="*") + outer(betahat,m$sd^2, FUN="*"))/outer(sebetahat^2,m$sd^2,FUN="+")
  ismissing = (is.na(betahat) | is.na(sebetahat))
  tmp[ismissing,]=m$mean #return prior mean when missing data
  t(tmp)
}


#return posterior standard deviation for each component of prior m, given observations betahat and sebetahat
#input, m is a mixture with k components
#betahat, sebetahat are n vectors
#output is a k by n matrix
comp_postsd.normalmix = function(m,betahat,sebetahat,v){
  if(!is.null(v)){
    stop("method comp_postsd of normal mixture not written for df!=NULL")
  }
  t(sqrt(outer(sebetahat^2,m$sd^2,FUN="*")/outer(sebetahat^2,m$sd^2,FUN="+")))
}

#find log likelihood of data, when 
#the mixture m is convolved with l-comp normal mixture
#in betahat with mixture sd betahatsd, mixture proportion pilik
#betahatsd is an n by l matrix
#betahat is an n vector
#v is the degree of freedom
#pilik is an l vector
loglik_conv_mixlik = function(m,betahat,betahatsd,v,pilik,FUN="+"){
  UseMethod("loglik_conv_mixlik")
}
loglik_conv_mixlik.default = function(m,betahat,betahatsd,v,pilik,FUN="+"){
  sum(log(dens_conv_mixlik(m,betahat,betahatsd,v,pilik,FUN)))
}

#compute the density of the components of the mixture m
#when convoluted with l-components normal mixture with standard deviation s
#or C (scale vector) multiplies scaled (se) student.t l-components mixture with df v
#with mixture proportion pilik
#the density is evaluated at x
#x is an n-vector
#s and pilik are n by l matrices
#v and c are l-vectors
#m is a mixture with k components
#output is a (k*l) by n matrix of densities
compdens_conv_mixlik = function(m,x,s,v,pilik,FUN="+"){
  UseMethod("compdens_conv_mixlik")
}
compdens_conv_mixlik.default = function(m,x,s,v,pilik,FUN="+"){
  dens=NULL
  for (i in 1:dim(pilik)[2]){
    dens=rbind(dens,pilik[,i]*compdens_conv(m,x,s[,i],v[i],FUN))
  }
  return(dens)
}

#compute density of mixture m convoluted with l-components
#normal mixture of sd (s) or student t mixture with df v
#with mixture proportion pilik
#at locations x
#m is a mixture
#x is an n vector
#s and pilik are n by l matrices
dens_conv_mixlik = function(m,x,s,v,pilik,FUN="+"){
  UseMethod("dens_conv_mixlik")
}
dens_conv_mixlik.default = function(m,x,s,v,pilik,FUN="+"){
  l=dim(pilik)[2]
  colSums(rep(m$pi,l) * compdens_conv_mixlik(m,x,s,v,pilik,FUN))
}

#compute the posterior prob that each observation
#came from each component of the mixture m
#output a k by n vector of probabilities
#computed by weighting the component densities by pi
#and then normalizing
#when likelihood is an l-components normal mixture or student t mixture
#with mixture proportion pilik
comppostprob_mixlik=function(m,x,s,v,pilik){
  UseMethod("comppostprob_mixlik") 
}
comppostprob_mixlik.default = function(m,x,s,v,pilik){
  l=dim(pilik)[2]
  k=length(m$pi)
  tmp= (t(rep(m$pi,l) * compdens_conv_mixlik(m,x,s,v,pilik))/dens_conv_mixlik(m,x,s,v,pilik))
  ismissing = (is.na(x) | apply(is.na(s),1,sum))
  tmp[ismissing,]=rep(m$pi,l)/l
  group=rep(1:k,l)
  return(rowsum(t(tmp),group))
}

#compute the posterior prob that each observation
#came from each component of the mixture m and the likelihood mixture
#output a (k*l) by n vector of probabilities
#computed by weighting the component densities by pi
#and then normalizing
#when likelihood is an l-components normal mixture or student t mixture
#with mixture proportion pilik
comppostprob_mixlik2=function(m,x,s,v,pilik){
  UseMethod("comppostprob_mixlik2") 
}
comppostprob_mixlik2.default = function(m,x,s,v,pilik){
  l=dim(pilik)[2]
  k=length(m$pi)
  tmp= (t(rep(m$pi,l) * compdens_conv_mixlik(m,x,s,v,pilik))/dens_conv_mixlik(m,x,s,v,pilik))
  ismissing = (is.na(x) | apply(is.na(s),1,sum))
  tmp[ismissing,]=rep(m$pi,l)/l
  return(t(tmp))
}

# evaluate cdf of posterior distribution of beta at c
# m is the prior on beta, a mixture
# c is location of evaluation
# assumption is betahat | beta \sim 
# l-components normal mixture with sd sebetahat
# and mixture proportion pilik
# m is a mixture with k components
# c a scalar
# betahat is an n vector 
# sebetahat and pilik are n by l matrices 
# output is a (k*l) by n matrix
compcdf_post_mixlik=function(m,c,betahat,sebetahat,v,pilik){
  UseMethod("compcdf_post_mixlik")
}
compcdf_post_mixlik.default=function(m,c,betahat,sebetahat,v,pilik){
  cdf=NULL
  for (i in 1:dim(pilik)[2]){
    cdf=rbind(cdf,compcdf_post(m,c,betahat,sebetahat[,i],v[i]))
  }
  cdf
}

cdf_post_mixlik = function(m,c,betahat,sebetahat,v,pilik){
  UseMethod("cdf_post_mixlik")
}
cdf_post_mixlik.default=function(m,c,betahat,sebetahat,v,pilik){
  colSums(comppostprob_mixlik2(m,betahat,sebetahat,v,pilik)*
            compcdf_post_mixlik(m,c,betahat,sebetahat,v,pilik))
}

#output posterior mean for beta for prior mixture m,
#given observations betahat, sebetahat, df v
#from l-components mixture likelihood
#with mixture proportion pilik
postmean_mixlik = function(m, betahat,sebetahat,v,pilik){
  UseMethod("postmean_mixlik")
}
postmean_mixlik.default = function(m,betahat,sebetahat,v,pilik){
  colSums(comppostprob_mixlik2(m,betahat,sebetahat,v,pilik) * comp_postmean_mixlik(m,betahat,sebetahat,v,pilik))
}

#output posterior mean-squared value for beta for prior mixture m,
#given observations betahat, sebetahat, df v
#from l-components mixture likelihood
#with mixture proportion pilik
postmean2_mixlik = function(m, betahat,sebetahat,v,pilik){
  UseMethod("postmean2_mixlik")
}
postmean2_mixlik.default = function(m,betahat,sebetahat,v,pilik){
  colSums(comppostprob_mixlik2(m,betahat,sebetahat,v,pilik) * comp_postmean2_mixlik(m,betahat,sebetahat,v,pilik))
}

#output posterior sd for beta for prior mixture m,
#given observations betahat, sebetahat, df v
#from l-components mixture likelihood
#with mixture proportion pilik
postsd_mixlik = function(m,betahat,sebetahat,v,pilik){
  UseMethod("postsd_mixlik")
}
postsd_mixlik.default = function(m,betahat,sebetahat,v,pilik){
  sqrt(postmean2_mixlik(m,betahat,sebetahat,v,pilik)-postmean_mixlik(m,betahat,sebetahat,v,pilik)^2)
}

#output posterior mean-squared value for beta for prior mixture m,
#given observations betahat, sebetahat, df v
#from l-components mixture likelihood
#with mixture proportion pilik
comp_postmean2_mixlik = function(m,betahat,sebetahat,v,pilik){
  UseMethod("comp_postmean2_mixlik")
}
comp_postmean2_mixlik.default = function(m,betahat,sebetahat,v,pilik){
  comp_postsd_mixlik(m,betahat,sebetahat,v,pilik)^2 + 
    comp_postmean_mixlik(m,betahat,sebetahat,v,pilik)^2
}

#output posterior mean for beta for each component of prior mixture m,
#given observations betahat, sebetahat, df v
#from l-components mixture likelihood
#with mixture proportion pilik
comp_postmean_mixlik=function(m,betahat,sebetahat,v,pilik){
  UseMethod("comp_postmean_mixlik")
}
comp_postmean_mixlik.default=function(m,betahat,sebetahat,v,pilik){
  mean=NULL
  for (i in 1:dim(pilik)[2]){
    mean=rbind(mean,comp_postmean(m,betahat,sebetahat[,i],v[i]))
  }
  return(mean)
}

#output posterior sd for beta for each component of prior mixture m,
#given observations betahat, sebetahat, df v
#from l-components mixture likelihood
#with mixture proportion pilik
comp_postsd_mixlik=function(m,betahat,sebetahat,v,pilik){
  UseMethod("comp_postsd_mixlik")
}
comp_postsd_mixlik.default=function(m,betahat,sebetahat,v,pilik){
  sd=NULL
  for (i in 1:dim(pilik)[2]){
    sd=rbind(sd,comp_postsd(m,betahat,sebetahat[,i],v[i]))
  }
  return(sd)
}


############################### METHODS FOR unimix class ###########################

#constructor; pi, a and b are vectors; kth component is Uniform(a[k],b[k])
unimix = function(pi,a,b){
  structure(data.frame(pi,a,b),class="unimix")
}

comp_cdf.unimix = function(m,y,lower.tail=TRUE){
  vapply(y,punif,m$a,min=m$a,max=m$b,lower.tail)
}

comp_sd.unimix = function(m){
  (m$b-m$a)/sqrt(12)
}

comp_mean.unimix = function(m){
  (m$a+m$b)/2
}



compdens.unimix = function(x,y,log=FALSE){
  k=ncomp(x)
  n=length(y)
  d = matrix(rep(y,rep(k,n)),nrow=k)
  return(matrix(dunif(d, x$a, x$b, log),nrow=k))  
}

#density of convolution of each component of a unif mixture with s*t_nu() at x
# x an n-vector
#return a k by n matrix
compdens_conv.unimix = function(m,x,s,v, FUN="+"){
  if(FUN!="+") stop("Error; compdens_conv not implemented for uniform with FUN!=+")
  if(is.null(v)){
    compdens= t(pnorm(outer(x,m$a,FUN="-")/s)-pnorm(outer(x,m$b,FUN="-")/s))/(m$b-m$a)
    compdens[m$a==m$b,]=t(dnorm(outer(x,m$a,FUN="-")/s)/s)[m$a==m$b,]
  }
  else{
    compdens= t(pt(outer(x,m$a,FUN="-")/s,df=v)-pt(outer(x,m$b,FUN="-")/s,df=v))/(m$b-m$a)
    compdens[m$a==m$b,]=t(dt(outer(x,m$a,FUN="-")/s,df=v)/s)[m$a==m$b,]
  }
  return(compdens)
}


#c is a scalar
#m a mixture with k components
#betahat a vector of n observations
#sebetahat an n vector of standard errors
#return a k by n matrix of the posterior cdf
compcdf_post.unimix=function(m,c,betahat,sebetahat,v){
  k = length(m$pi)
  n=length(betahat)
  tmp = matrix(1,nrow=k,ncol=n)
  tmp[m$a > c,] = 0
  subset = m$a<=c & m$b>c # subset of components (1..k) with nontrivial cdf
  if(sum(subset)>0){
    if(is.null(v)){
      pna = pnorm(outer(betahat,m$a[subset],FUN="-")/sebetahat)
      pnc = pnorm(outer(betahat,rep(c,sum(subset)),FUN="-")/sebetahat)
      pnb = pnorm(outer(betahat,m$b[subset],FUN="-")/sebetahat)
    }else{
      pna = pt(outer(betahat,m$a[subset],FUN="-")/sebetahat, df=v)
      pnc = pt(outer(betahat,rep(c,sum(subset)),FUN="-")/sebetahat, df=v)
      pnb = pt(outer(betahat,m$b[subset],FUN="-")/sebetahat, df=v)
    }
    tmp[subset,] = t((pnc-pna)/(pnb-pna))
  }
  subset = (m$a == m$b) #subset of components with trivial cdf
  tmp[subset,]= rep(m$a[subset] <= c,n)
  tmpnaive=matrix(rep((c-m$a)/(m$b-m$a),length(betahat)),nrow=k,ncol=n)
  tmp[is.nan(tmp)]= tmpnaive[is.nan(tmp)]
  tmp
}

my_etruncnorm= function(a,b,mean=0,sd=1){
  alpha = (a-mean)/sd
  beta =  (b-mean)/sd
  #Flip the onese where both are positive, as the computations are more stable
  #when both negative
  flip = (alpha>0 & beta>0)
  flip[is.na(flip)]=FALSE #deal with NAs
  alpha[flip]= -alpha[flip]
  beta[flip]=-beta[flip]
  
  #Fix a bug of quoting the truncnorm package
  #E(X|a<X<b)=a when a==b as a natural result
  #while etruncnorm would simply return NaN,causing PosteriorMean also NaN
  tmp1=etruncnorm(alpha,beta,0,1)
  isequal=(alpha==beta)
  tmp1[isequal]=alpha[isequal]
  
  tmp= (-1)^flip * (mean+sd*tmp1)
  
  max_alphabeta = ifelse(alpha<beta, beta,alpha)
  max_ab = ifelse(alpha<beta,b,a)
  toobig = max_alphabeta<(-30)
  toobig[is.na(toobig)]=FALSE 
  tmp[toobig] = max_ab[toobig]
  tmp
}


#return posterior mean for each component of prior m, given observations betahat and sebetahat
#input, m is a mixture with k components
#betahat, sebetahat are n vectors
#output is a k by n matrix
#note that with uniform prior, posterior is truncated normal, so
#this is computed using formula for mean of truncated normal 
comp_postmean.unimix = function(m,betahat,sebetahat,v){
  #   k= ncomp(m)
  #   n=length(betahat)
  #   a = matrix(m$a,nrow=n,ncol=k,byrow=TRUE)
  #   b = matrix(m$b,nrow=n,ncol=k,byrow=TRUE)
  #   matrix(etruncnorm(a,b,betahat,sebetahat),nrow=k,byrow=TRUE)
  #note: etruncnorm is more stable for a and b negative than positive
  #so maybe use this, and standardize to make the whole more stable.
  
  alpha = outer(-betahat, m$a,FUN="+")/sebetahat
  beta = outer(-betahat, m$b, FUN="+")/sebetahat
  if(is.null(v)){
    tmp = betahat + sebetahat*my_etruncnorm(alpha,beta,0,1)
  }else{
    tmp = betahat + sebetahat*my_etrunct(alpha,beta,v)
  }
  ismissing = is.na(betahat) | is.na(sebetahat)
  tmp[ismissing,]= (m$a+m$b)/2
  t(tmp)
  #   t(
  #     betahat + sebetahat* 
  #       exp(dnorm(alpha,log=TRUE)- pnorm(alpha,log=TRUE))
  #    * 
  #       (-expm1(dnorm(beta,log=TRUE)-dnorm(alpha,log=TRUE)))
  #     /
  #       (expm1(pnorm(beta,log=TRUE)-pnorm(alpha,log=TRUE)))
  #   )
}

#not yet implemented!
#just returns 0s for now
comp_postsd.unimix = function(m,betahat,sebetahat,v){
  print("Function ashci() is provided for computing the credible interval(HPD),see documentation for usage and example.")
  k= ncomp(m)
  n=length(betahat)
  return(matrix(NA,nrow=k,ncol=n)) 
}

# the mean of a truncated student.t
# the result is from the paper 'Moments of truncated Student-t distribution' by H.-J Kim 

my_etrunct= function(a,b,v){
  A = v+a^2
  B = v+b^2
  F_a = pt(a,df=v)
  F_b = pt(b,df=v)
  G = gamma((v-1)/2)*v^(v/2)/(2*(F_b-F_a)*gamma(v/2)*gamma(1/2))
  tmp = ifelse(a==b,a,G*(A^(-(v-1)/2)-B^(-(v-1)/2)))
  tmp
}


############################### METHODS FOR igmix class ###########################

#' @title Constructor for igmix class
#'
#' @description Creates an object of class igmix (finite mixture of univariate inverse-gammas)
#' 
#' @details None
#' 
#' @param pi: vector of mixture proportions
#' @param alpha: vector of shape parameters
#' @param beta: vector of rate parameters
#' 
#' @return an object of class igmix
#' 
#' @examples igmix(c(0.5,0.5),c(1,1),c(1,2))
#' 
igmix = function(pi,alpha,beta){
  structure(data.frame(pi,alpha,beta),class="igmix")
}

comp_sd.igmix = function(m){
  m$beta/(m$alpha-1)/sqrt(m$alpha-2)
}

comp_mean.igmix = function(m){
  m$beta/(m$alpha-1)
}

compdens.igmix = function(x,y,log=FALSE){
  k=ncomp(x)
  n=length(y)
  d = matrix(rep(y,rep(k,n)),nrow=k)
  return(matrix(dgamma(1/d, shape=x$alpha, rate=outer(m$beta,1/y^2),log),nrow=k))  
}

#density of product of each component of a inverse-gamma mixture with Gamma(v/2,v/2) at s
# s an n-vector at which density is to be evaluated
#return a k by n matrix
compdens_conv.igmix = function(m,x,s,v,FUN="+"){  
  k=ncomp(m)
  n=length(s)
  dens = t(exp(v/2*log(v/2)-lgamma(v/2)
               +(v/2-1)*outer(log(s^2),rep(1,k))
               +outer(rep(1,n),m$alpha*log(m$beta)-lgamma(m$alpha)+lgamma(m$alpha+v/2))
               -outer(rep(1,n),m$alpha+v/2)*log(outer(v/2*s^2,m$beta,FUN="+"))))
  return(dens)
  
}

comp_cdf.igmix = function(x,y,lower.tail=TRUE){
  vapply(y,pigamma,x$alpha,x$alpha,x$beta,lower.tail)
}

#c is a scalar
#m a mixture with k components
#betahat a vector of n observations
#sebetahat an n vector of standard errors
#return a k by n matrix of the posterior cdf
compcdf_post.igmix=function(m,c,betahat,sebetahat,v){
  #compute posterior shape (alpha1) and rate (beta1)
  alpha1 = m$alpha+v/2
  beta1 = outer(m$beta,v/2*sebetahat^2,FUN="+")
  ismissing = is.na(sebetahat)
  beta1[,ismissing]=m$beta
  
  return(t(pigamma(c,alpha1,beta1)))
}

#return posterior mean for each component of prior m, given observations sebetahat
#input, m is a mixture with k components
#betahat, sebetahat are n vectors
#output is a k by n matrix
comp_postmean.igmix = function(m,betahat,sebetahat,v){
  k = length(m$pi)
  n=length(sebetahat)
  tmp=outer(v/2*sebetahat^2,m$beta,FUN="+")/outer(rep(1,n),m$alpha+v/2-1)
  ismissing = is.na(sebetahat)
  tmp[ismissing,]=m$beta/(m$alpha-1) #return prior mean when missing data
  t(tmp)
}


#return posterior standard deviation for each component of prior m, given observations betahat and sebetahat
#input, m is a mixture with k components
#betahat, sebetahat are n vectors
#output is a k by n matrix
comp_postsd.igmix = function(m,betahat,sebetahat,v){
  k = length(m$pi)
  n=length(sebetahat)
  alpha1 = outer(rep(1,n),m$alpha+v/2-1)
  beta1 = outer(v/2*sebetahat^2,m$beta,FUN="+")
  return(t(beta1/(alpha1-1)/sqrt(alpha1-2)))
}
