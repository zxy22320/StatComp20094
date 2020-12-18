#' @title Upper bound for the correlation of correlated binary distribution
#' @description The upper bound for the pairwise correlation in exchange correlated binary distribution. 
#' @param p marginal frequencies (numeric vector)
#' @return The upper bound for the pairwise correlation \code{rhoMax}
#' @examples
#' \dontrun{
#' p=c(0.1,0.2,0.3)
#' rhoMaxEx(p)
#' }
#' @export


rhoMaxEx<-function(p){
	#The upper bound for the pairwise correlation in exchange correlated binary distribution with marginal frequencies p
	if(!is.atomic(p) || typeof(p)!='double') return(NaN)
	if(sum((p<0) | (p>1))!=0) return(NaN)
	minP<-min(p)
	maxP<-max(p)
	rhoMax<-sqrt((minP/(1-minP))/(maxP/(1-maxP)))
	return(rhoMax)
}


#' @title Generating One Dimension Correlated Binary Variables
#' @description A function to generate one dimension correlated bernoulli distribution variables with certain correlation and success probability. 
#' @importFrom stats rbinom
#' @param n length of variables (numeric)
#' @param p success probility (numeric)
#' @param rho correlationship (numeric)
#' @return Correlated bernoulli distribution variables \code{X}
#' @examples
#' \dontrun{
#' rcbern0(10,0.5,0.2)
#' }
#' @export


rcbern0<-function(n, p, rho){
	#Generate Correlated Bernoulli Distribution
	#See 'A Note on Generating Correlated Binary Variables'
	#A. D. Lunn and S. J. Davies
	#Biometrika
	#Vol. 85, No. 2 (Jun., 1998), pp. 487-490 
	if(!is.atomic(p) || typeof(p)!='double') return(NaN)
	if(length(p)!=1) return(NaN)
	if((p<0) | (p>1)) return(NaN)
	
	if((rho<0) || (rho>1)) return(NaN)
	Z<-rbinom(1, 1, p)
	Y<-rbinom(n, 1, p)
	U<-rbinom(n, 1, sqrt(rho))
	X<-(1-U)*Y+U*Z

	return(X)
}


#' @title Generating High Dimension Correlated Binary Variables
#' @description A function to generate high dimension correlated bernoulli distribution variables with certain correlation and success probability. 
#' @importFrom stats rbinom
#' @param n length of variables (numeric)
#' @param p success probility (numeric vector)
#' @param rho correlationship (numeric)
#' @return Matrix of correlated bernoulli distribution variables \code{X}
#' @examples
#' \dontrun{
#' p<-c(0.1,0.2,0.3)
#' rcbern(10,p,0.2)
#' }
#' @export


rcbern<-function(n, p, rho){
  #Generate Correlated Bernoulli Distribution
  #See "A Note on Generating Correlated Binary Variables"
  #A. D. Lunn and S. J. Davies
  #Biometrika
  #Vol. 85, No. 2 (Jun., 1998), pp. 487-490 
  if(!is.atomic(p) || typeof(p)!='double') return(NaN)
  if(sum((p<0) | (p>1))!=0) return(NaN)
  
  m<-length(p)
  minP<-min(p)
  maxP<-max(p)
  rhoLimit<-sqrt((minP/(1-minP))/(maxP/(1-maxP)))
  if((rho<0) || (rho>rhoLimit)){
    cat(paste('The range of rho is [',0,',',round(rhoLimit,3),']\n', sep=''))
    cat('rho is out-of-range\n')
    return(NaN)
  }
  Pc<-sqrt(minP*maxP)/(sqrt(minP*maxP)+sqrt((1-minP)*(1-maxP)))
  Pa<-sqrt(rho*p*(1-p)/(Pc*(1-Pc)))
  if(max(Pa)==1){
    X<-replicate(n, rep(rbinom(1, 1, p[1]), m))
  }else{
    Pb<-(p-Pa*Pc)/(1-Pa)
    X<-replicate(n, {
      U<-rbinom(m, 1, Pa)
      Y<-rbinom(m, 1, Pb)
      Z<-rbinom(1, 1, Pc)
      (1-U)*Y+U*Z
    })
  }
  X<-t(X)
  return(X)
}



#' @title Generating correlated genotypes
#' @description A function to generate genotypes for correlated individuals at SNPs with certain MAFs.
#' @importFrom stats rbinom
#' @useDynLib StatComp20094
#' @param n sample size (numeric)
#' @param f MAFs (numeric vector)
#' @param rho correlationship (numeric)
#' @return Genotype Matrix \code{X}
#' @examples
#' \dontrun{
#' f<-rep(0.2,10)
#' simGenoData(2,f,0.1)
#' }
#' @export



simGenoData<-function(n, f, rho){
	m<-length(f)
	X<-matrix(0, n, m)
	for(i in 1:m){
		if(rho!=0){
			#X[,i]<-colSums(rmvbin(2, rep(f[i], n), bincorr=(1-rho)*diag(n)+rho))		
			X[,i]<-rcbern0(n, f[i], rho)+rcbern0(n, f[i], rho)
		}else{
			X[,i]<-rbinom(n,2,f[i])
		}
	}
	return(X)
}


#' @title Gwas simulation data
#' @description A function to generate genotypes when both cases and controls are correlated.
#' @importFrom stats rbinom
#' @useDynLib StatComp20094
#' @param n0 case size (numeric)
#' @param n1 control size (numeric)
#' @param f0 MAFs of cases (numeric vector)
#' @param f1 MAFs of controls (numeric vector) 
#' @param rho correlationship between all individuals (numeric)
#' @return Genotype Matrix \code{X}
#' @return label of case or control \code{y}
#' @examples
#' \dontrun{
#' f0<-f1<-rep(0.2,100)
#' simGwasData(10,10,f0,f1,0.1)
#' }
#' @export


#Both cases and controls are correlated
simGwasData<-function(n0, n1, f0, f1, rho){
  n<-n0+n1
  m<-length(f0)
  X<-matrix(0, n, m)
  for(i in 1:m){
    if(rho!=0){
      #X[,i]<-colSums(rmvbin(2, rep(f[i], n), bincorr=(1-rho)*diag(n)+rho))		
      f<-c(rep(f0[i],n0), rep(f1[i], n1))
      X[,i]<-rcbern(1, f, rho)+rcbern(1, f, rho)  
    }else{
      X[,i]<-c(rbinom(n0, 2, f0[i]), rbinom(n1, 2, f1[i]))
    }
  }
  y<-c(rep(0, n0), rep(1, n1))
  return(list(X=X, y=y))
}


#' @title Gwas simulation data I
#' @description A function to generate genotypes when the correlation is the same among cases and controls, but no correlations between a case and a control.
#' @importFrom stats rbinom
#' @useDynLib StatComp20094
#' @param n0 case size (numeric)
#' @param n1 control size (numeric)
#' @param f0 MAFs of cases (numeric vector)
#' @param f1 MAFs of controls (numeric vector) 
#' @param rho correlationship within cases and controls (numeric)
#' @return Genotype Matrix \code{X}
#' @return label of case or control \code{y}
#' @examples
#' \dontrun{
#' f0<-f1<-rep(0.2,100)
#' simGwasDataI(10,10,f0,f1,0.1)
#' }
#' @export



#No correlations between cases and controls
simGwasDataI<-function(n0, n1, f0, f1, rho){
  n<-n0+n1
  m<-length(f0)
  X<-matrix(0, n, m)
  for(i in 1:m){
    if(rho!=0){
      #X[,i]<-colSums(rmvbin(2, rep(f[i], n), bincorr=(1-rho)*diag(n)+rho))
      f0i<-rep(f0[i], n0)
      f1i<-rep(f1[i], n1)
      X[,i]<-c(rcbern(1, f0i, rho)+rcbern(1, f0i, rho), rcbern(1, f1i, rho)+rcbern(1, f1i, rho))
    }else{
      X[,i]<-c(rbinom(n0, 2, f0[i]), rbinom(n1, 2, f1[i]))
    }
  }
  y<-c(rep(0, n0), rep(1, n1))
  return(list(X=X, y=y))
}



#' @title Gwas simulation data II
#' @description A function to generate genotypes when different correlations are used among cases and controls, but no correlations between a case and a control.
#' @importFrom stats rbinom
#' @useDynLib StatComp20094
#' @param n0 case size (numeric)
#' @param n1 control size (numeric)
#' @param f0 MAFs of cases (numeric vector)
#' @param f1 MAFs of controls (numeric vector) 
#' @param rho0 correlationship between cases (numeric)
#' @param rho1 correlationship between controls (numeric)
#' @return Genotype Matrix \code{X}
#' @return label of case or control \code{y}
#' @examples
#' \dontrun{
#' f0<-f1<-rep(0.3,100)
#' simGwasDataII(10,10,f0,f1,0.1,0.2)
#' }
#' @export


#Different correlations used among cases and controls, but no correlations between a case and a control
simGwasDataII<-function(n0, n1, f0, f1, rho0, rho1){
  n<-n0+n1
  m<-length(f0)
  
  X<-matrix(0, n, m)
  for(i in 1:m){
    if(rho0!=0 || rho1!=0){
      #X[,i]<-colSums(rmvbin(2, rep(f[i], n), bincorr=(1-rho)*diag(n)+rho))
      f0i<-rep(f0[i], n0)
      f1i<-rep(f1[i], n1)
      X[,i]<-c(rcbern(1, f0i, rho0)+rcbern(1, f0i, rho0), rcbern(1, f1i, rho1)+rcbern(1, f1i, rho1))
    }
    else{
      X[,i]<-c(rbinom(n0, 2, f0[i]), rbinom(n1, 2, f1[i]))
    }
  }
  y<-c(rep(0, n0), rep(1, n1))
  return(list(X=X, y=y))
  # return(X)
}







