library(splines)
library(sandwich)
library(lmtest)
library(Matrix)

# lagmatrix <- function( X , lag ){
#   X <- as.matrix(X)
#   N <- ncol(X)
#   T <- nrow(X)
#   X <- rbind( matrix(NA,lag,N) , X )
#   X <- embed(X,lag)
#   #X[ is.na(X) ] <- 0
#   X[1:T,]
# }

lagmatrix <- function( X , lag ){
  if (lag>0) {
    X <- as.matrix(X)
    N <- ncol(X)
    T <- nrow(X)
    X <- rbind( matrix(NA,lag,N) , X )
    X <- embed(X,lag)
    #X[ is.na(X) ] <- 0
    X[1:T,]
  }
}


#
lproj <- function( y , x , w=NULL , const=TRUE , type='reg' , H , h1 , r=0 , zero=FALSE , lambda=0, trace=1 ){
  
  # dimensions
  T   <- length(y)
  
  # construct basis
  if( type=='smooth' ){
    h.range <- h1:H
    bdeg    <- 3
    knots   <- seq(-bdeg+h1,H+bdeg,1)
    basis   <- spline.des(knots, h.range, bdeg + 1, 0 * h.range ,outer.ok=TRUE)$design
  }
  
  # construct w
  if( const==TRUE ){
    w <- cbind( rep(1,T) , w )
  }
  
  # 1 shock std dev
  if( !is.null(w) ){
    delta <- summary(lm(x ~ 0+w))$sigma
  }
  else{
    delta <- sd(x)
  }  
  
  # dimensions
  HR  <- H+1-h1
  TS  <- T*HR
  if( type=='reg' ){
    XS  <- HR
  }
  else{
    XS  <- ncol(basis) 
  }
  WS  <- HR  
  if( !is.null(w) | const==TRUE ){
    NW  <- ncol(w)  
  }
  else{
    NW  <- 0
  }
  
  # y
  IDX <- matrix(0,TS,2)
  Y   <- rep(NA,TS)
  Xb  <- matrix(0,TS,XS)  
  Xc  <- array(0,dim=c(TS,HR,NW))
  
  # construct Y and X
  for( t in 1:(T-h1) ){      
    idx.beg <- (t-1)*HR + 1
    idx.end <- t*HR
    
    IDX[ idx.beg:idx.end , 1 ] <- t
    IDX[ idx.beg:idx.end , 2 ] <- h1:H
    
    # y
    y.range <- (t+h1) : min((t+H),T)
    Y [ idx.beg:idx.end ] <- c( y[ y.range ] , rep(NA,HR-length(y.range)) )
    
    # x
    if( type=='reg' ){ 
      Xb[ idx.beg:idx.end , ] <- diag(HR)*x[t]      
    }
    else {
      Xb[ idx.beg:idx.end , ] <- basis*x[t]
    }
    
    # w
    for( i in seq_len(NW) ){
      Xc[ idx.beg:idx.end ,  , i ] <- diag(HR)*w[t,i]
    }
  }
  
  X   <- cbind(Xb)
  for( i in seq_len(NW)){
    X <- cbind(X,Xc[,,i])
  }
  X  <- Matrix::Matrix(X,sparse=TRUE)
  
  sel <- !is.na(Y)
  IDX<- IDX[sel,]
  Y  <- Y[sel]
  X  <- X[sel,]
  TS <- length(Y)
  
  XX <- t(X)%*%X
  XY <- t(X)%*%Y
  
  # penalty
  P <- matrix(0,ncol(X),ncol(X))
  P <- Matrix(P,sparse=TRUE)
  
  if( type=='smooth' ){
    D   <- diag(XS) 
    for (k in seq_len(r)) D <- diff(D)
    
    if( zero ){
      DP     <- rep(0,XS)
      DP[XS] <- 1
      D      <- rbind(D,DP)
    }
    
    P[1:XS,1:XS] <- t(D) %*% D 
  }
  
  ir    <- matrix(0,H+1,length(lambda))
  theta <- matrix(0,ncol(X),length(lambda))
  mul   <- matrix(0,HR,length(lambda))
  
  for( i in 1:length(lambda) ){
    if (trace==1) {cat('.')}
    A         <- XX + lambda[i]*TS*P 
    b         <- XY
    theta[,i] <- as.vector( Matrix::solve( A , b ) )
    
    if( type=='reg' ){
      mul[,i]   <- theta[1:XS,i]
    }
    else{
      beta  <- theta[1:XS,i]
      mul[,i]   <- as.matrix(basis) %*% as.vector(beta)     
    }
    
    ir[(h1+1):(H+1),i]   <- mul[,i]*delta
  }
  if (trace==1) {cat('\n')}
  
  obj <- list()
  obj$type<- type
  if( type=='smooth'){
    obj$basis <- basis
  }
  obj$h1     <- h1
  obj$H      <- H
  obj$XS     <- XS
  obj$HR     <- HR
  obj$T      <- T
  obj$H      <- H
  obj$TS     <- TS
  obj$IDX    <- IDX
  obj$y      <- y 
  obj$x      <- x
  obj$w      <- w
  obj$Y      <- Y
  obj$X      <- X
  obj$theta  <- theta
  obj$mul    <- mul
  obj$lambda <- lambda
  obj$P      <- P
  obj$ir     <- ir
  obj$delta  <- delta
  
  obj
}

lproj.conf <- function( obj , l=1 ){
  
  u <- obj$Y - obj$X %*% obj$theta[,l];    
  S <- obj$X * ( u %*% t(rep(1,ncol(obj$X))) )
  
  # BREAD
  bread <- solve( t(obj$X)%*%obj$X + obj$lambda[l]*obj$TS*obj$P )
  
  # MEAT
  nlag    = min(floor(1.2*(obj$T)**(1/3)),obj$T)
  nlag    = obj$H
  weights = (nlag+1-(0:nlag))/(nlag+1)
  V <- t(S) %*% S
  for( i in 1:nlag ){
    Gammai      <- t( S[(i+1):obj$T,] ) %*% S[1:(obj$T-i),]
    GplusGprime <- Gammai+t(Gammai)
    V           <- V+weights[i+1]*GplusGprime
  }
  
  meat <- V
  
  V  <- bread %*% meat %*% bread
  
  if( obj$type == 'reg' ){
    se   <- sqrt( diag( V[ 1:(obj$H+1-obj$h1) , 1:(obj$H+1-obj$h1) ] ) )
    conf <- matrix( 0 , length(se) , 2 )
    conf[,1] <- obj$mul[,l] + se*qnorm(0.10)      
    conf[,2] <- obj$mul[,l] + se*qnorm(0.90)      
    
  }
  else{
    V    <- as.matrix(obj$basis) %*% V[ 1:obj$XS , 1:obj$XS ] %*% t(as.matrix(obj$basis))
    se   <- sqrt( diag( V ) )
    conf <- matrix( 0 , length(se) , 2 )
    conf[,1] <- obj$mul[,l] + se*qnorm(0.10)      
    conf[,2] <- obj$mul[,l] + se*qnorm(0.90)      
    #conf[nrow(conf),] <- NA
  }
  
  irc <- matrix(NA,obj$H+1,2)
  irc[(1+obj$h1):(obj$H+1),] <- conf*obj$delta
  
  obj$se  <- se
  obj$irc <- irc
  
  obj
}

lproj.cv <- function( obj , K, trace=1 ){
  
  T   <- obj$T
  L   <- length(obj$lambda)
  
  ind <- ceiling( (obj$IDX[,1]/T)*K )
  rss <- rep(0,L)
  
  for( l in 1:L ){
    if (trace==1) {cat('.')}
    
    rss.l <- rep(0, K)
    
    for( i in 1:K ){
      
      #Y.in   <- obj$Y[ obj$IDX[,1] != i ] 
      #X.in   <- obj$X[ obj$IDX[,1] != i , ]
      #Y.out  <- obj$Y[ obj$IDX[,1] == i ] 
      #X.out  <- obj$X[ obj$IDX[,1] == i , ]
      
      Y.in   <- obj$Y[ ind != i ] 
      X.in   <- obj$X[ ind != i , ]
      Y.out  <- obj$Y[ ind == i ] 
      X.out  <- obj$X[ ind == i , ]
      
      A <- t(X.in)%*%X.in + obj$lambda[l] * obj$TS * ((K-1)/K) * obj$P
      b <- t(X.in)%*%Y.in  
      beta     <- Matrix::solve(A,b)  
      rss.l[i] <- mean( ( Y.out - X.out %*% beta )**2 )
    }
    
    rss[l] <- mean(rss.l)
  }
  if (trace==1) {cat('\n')}
  
  obj$rss     <- rss
  obj$idx.opt <- tail(which(min(rss)==rss),1)
  obj$ir.opt  <- obj$ir[ , tail(which(min(rss)==rss),1) ]
  
  obj
}
