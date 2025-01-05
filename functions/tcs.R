
 # Description:   A trend cycle season filter, translated from Matlab code by Mohr 
 #                translated using mat2r() from library(matconv)

 # Arguments:
 #       y        - time series to be filtered
 #       nd       - level of the trend (usually in the range of 0 - 2)
 #       nc       - level of the 1st cyclical process (usually in the range of 0 - 2)
 #       cyc      - the critical cycle for the first stochastic AR-cycle model, for nc = 0
 #       s        - the season: 4 for quarterly data, 12 for monthly data, etc., s = 1 (default), no season
 #       du       - dummy matrix for breaks 
 #       rho      - dampening parameter for the stochastic AR-cycle model 
 
 tcs <- function(y, nd = 2, nc = 0, cyc = 1600, s = 1, du = NULL, rho = 0.975){
 
 N   <- length(y)
 trs <- s # trs: seasonal difference; trs=1 -> trend does not contain seasonal root
 U <- matrix(0, N, 0) # U contains dummies and/or unity vector if nd = 1
 if (nd > 0L){# order of trend > 0 -> stochastic trend is specified
 D <- dmat(N,trs,nd)
 if (nd == 1L){# nd = 1 -> first order stochastic trend with deterministic drift
 U <- matrix(1, N, 1)# U contains unity vector to capture the deterministic drift constant
 U[1] <- 0
 }
 if (!is.null(du)){ # structural breaks in the trend specified by dummy matrix du
 U <- cbind(U, D%*%du) # write dummy matrix du into U
 }
 if (ncol(U) > 0L){# U exists and is not empty -> either structural breaks or constant drift term
 W <- U%*%solve(t(U)%*%U)%*%t(U)
 T <- t(D)%*%(diag(N) - W)%*%D
 } else {# neither structural breaks nor drift
 T <- t(D)%*%D
 }
 } else {
 D <- matrix(0, N, N)# order of trend = 0 -> no trend specified
 T <- matrix(0, N, N)
 }
 if (nc == 0L & cyc > 0){# this captures the case of the HP filter or the EES
 T <- cyc * T# if nc== 0 and cyc > 0, cyc contains the value of lambda of the HP filter or the EES
 }
 if (nc > 0){# nc > 0 -> stochastic cycle is specified
 mue  <- 2*pi/cyc  # mean length of cycle in the same frequency as data (quarters, months, etc.)
 arma <- matrix(c(1, -2*rho*cos(mue), rho^2, 
 1, -rho*cos(mue), 0), 2, 3, byrow = T) # create A and B, the ar ant ma matrices of the stochastic cycle
 arma.mat <- armatomat(arma,N)
 A <- arma.mat$armat
 B <- arma.mat$mamat
 A <- pomat(A, nc)
 B <- pomat(B, nc)
 A <- A[-c(1:2*nc),]
 B <- B[-c(1:2*nc),]
 C <- t(A)%*%ginv(B%*%t(B))%*%A # create C
 } else {
 C <- matrix(0, N, N)# order of cycle =0 -> no stochastic cycle specified
 }
 if (s > 1){# if a sesonsal process was specified, compute the matrix S for the stochastic seasonal process
 seasm      <- seasonalmat(s)
 seasm.mat  <- armatomat(seasm, N)
 P          <- seasm.mat$armat[-c(1:(s-1)),]
 Q          <- seasm.mat$mamat[-c(1:(s-1)),]
 S          <- t(P)%*%ginv(Q%*%t(Q))%*%P
 } else {
 S <- matrix(0, N, N)# s<1 -> 2-component seasonal filters collapse to NxN zero matrix
 }
 # now compute the 3-component filter matrices
 MTCS <- f3mat(T,C,S)
 MSTC <- f3mat(S,T,C)
 if (nc > 0L){
 MCTS <- f3mat(C,T,S)
 } else {# order of stochastic cycle zero,
 if (cyc > 0){# but cycle > 0: this is the case of the HP filter or the EES
 MCTS <- diag(N) - MTCS - MSTC# which means that the cycle is just the residual
 } else {# and cycle <=0: there is no cycle
 MCTS <- matrix(0, N, N)
 }
 }
 trend  <- MTCS%*%y  # compute the trend component
 cycle  <- MCTS%*%y  # compute the cyclical component
 season <- MSTC%*%y  # compute the seasonal component
 data   <- ts.union(y, trend, cycle, season)
 if (ncol(U) > 0L){# no-empty U means that
 # there is a constant drift (i.e. nd=1) or there are structural breaks
 # -> compute the drift (if there is one) and the size of  the structural breaks
 v <-  solve(t(U)%*%U)%*%t(U)%*%D%*%trend # -> drift and structural breaks written into v
 } else {
 v <- NULL
 }
 return(list(data = data, v = v))
 }
 #============================================================================
 
 pomat <- function(A, n){
 if(n == 0L) return(diag(ncol(A)))
 if(n == 1L) return(A)
 PA <- A # raise matrix A to power n
 for(f in 1:(n-1)){
 PA <- PA%*%A
 }
 return(PA)
 }
 
 
 #============================================================================
 f1mat <- function(A){
 #----------------------------------------------------------------------------
 # provides 1-component filter matrix MA
 #----------------------------------------------------------------------------
 N <- ncol(A)
 if (all(A == 0)){
 f1m <- matrix(0, N, N)
 } else {
 f1m <- ginv(diag(N) + A)
 }
 return(f1m)
 }
 #============================================================================
 
 #============================================================================
 f2mat <- function(A,B){
 #----------------------------------------------------------------------------
 # provides 2-component filter matrix MAB
 #----------------------------------------------------------------------------
 N <- ncol(A)
 if (all(A == 0)){
 f2m <- matrix(0, N, N)
 } else {
 MA <- f1mat(A)
 MB <- f1mat(B)
 f2m <- ginv(diag(N)-MA%*%MB)%*%MA%*%(diag(N)-MB)
 }
 return(f2m)
 }
 #============================================================================
 
 #============================================================================
 f3mat <- function(A,B,C){
 #----------------------------------------------------------------------------
 # provides 3-component filter matrix MABC
 #----------------------------------------------------------------------------
 N <- ncol(A)
 if (all(A == 0)){
 f3m <- matrix(0, N, N)
 } else {
 MAB <- f2mat(A,B)
 MCB <- f2mat(C,B)
 f3m <- ginv(diag(N)-MAB%*%MCB)%*%MAB%*%(diag(N)-MCB)
 }
 return(f3m)
 }
 #============================================================================
 
 #============================================================================
 seasonalmat <- function(s){
 #----------------------------------------------------------------------------
 # Builds a 2xs matrix with the parameters of the Pauly/Schlicht (Empirica 1983) seasonal ARMA process
 # The first row of this matrix contains the vector of ar-parameters and
 # the second row contains the vector of ma-parameters.
 # The implied ARMA process is p(L)x_s=q(L)e, where e is white noise and x_s is the seaonal process
 # p(L):=sum[L^i,{i,0,s-1}], q(L):=(1/s)Sum[(s - i) L^(i - 1), {i, 1, s - 1}]
 # The parameters for lag zero appear in the first column of the arma matrix.
 # s[scalar]: the seasonal frequency (i.e. s for quarterly data, 12 for monthly data)
 # seasmmatrix]: the ARMA process output matrix The first row contains the AR parameters,
 # the second row contains the MA params. Note that the parameters
 # for lag zero appear in the first column of ARMA
 # Needs:
 # nothing
 #--------------------------------------------------------------------------
 if (s==0){
 seasm <- 0
 } else {
 seasm <- rbind(rep(1,s), (s-(1:s))/s)
 }
 return(seasm)
 }
 #============================================================================
 
 #============================================================================
 armatomat <- function(arma, T){
 #----------------------------------------------------------------------------
 # Builds an ar matrix and an ma matrix out of a 2xk matrix
 # which first row contains the vector of ar-parameters and
 # which second row contains the vector of ma-parameters.
 # The implied ARMA process is A(L)y=B(L)e, where e is white noise
 # Note that the arma matrix includes the parameters for lag zero in the first column
 # which - in most cases - will be equal to one
 # arma[matrix]: 2xk matrix. The first row contains the AR params,
 # the second row contains the MA params. Note that the parameters
 # include the parameters for lag zero in the first column of arma
 # and will ususally be equal to 1
 # T[scalar]: Integer, size of the ar and ma matrices
 # armat[matrix]: the ar output matrix
 # mamat[matrix]: the ma output matrix
 #--------------------------------------------------------------------------
 n  <- nrow(arma)
 k  <- ncol(arma)
 ar <- armat <- c(arma[1, k:1], rep(0, T - k)) #revert order
 ma <- mamat <- c(arma[2, k:1], rep(0, T - k))
 for(i in 1:(T - k)){
 ar <- c(0, ar[-T])
 ma <- c(0, ma[-T])
 armat <- rbind(armat, ar)
 mamat <- rbind(mamat, ma)
 }
 zeros <- matrix(0, nrow = ncol(armat) - nrow(armat), ncol = T)
 armat <- rbind(zeros, armat)
 mamat <- rbind(zeros, mamat)
 rownames(armat) <- rownames(mamat) <- NULL
 return(list(armat = armat, mamat = mamat))
 }
 #============================================================================
 
 #============================================================================
 lmat <- function(T, n){
 #--------------------------------------------------------------------------
 # Description
 # Gives a TxT lag matrix Lmat with lag n such that
 # lmat(T)*X gives n-lagged values of X (X being an arbitrary vector)
 # Usage:
 # lmat(T, n)
 # T[scalar]: Integer, Dimension of the TxT output matrix
 # n[scalar]: Integer, lag length
 # Needs:
 # Nothing
 #--------------------------------------------------------------------------
 if(T <= n) {
 warning("Time series are too short!")
 return(0*diag(T))}
 
 if(n == 0){return(diag(T))}
 
 lm <- cbind(rbind(matrix(0, n, T - n), diag(T - n)), matrix(0, T, n))
 return(lm)
 
 }
 #============================================================================
 
 #============================================================================
 dmat <- function(T, s, n){
 #----------------------------------------------------------------------------
 # Gives a TxT nth-difference matrix dmat such that
 # dmat*X produces the n-th Differences of season sss of X (X being an arbitrary vector): (1-L^[sss])^n X[t]
 # T[scalar]: Integer, dimension of the TxT output matrix
 # s[scalar]: Integer, order of season (e.g.: 1 for annual data, 4 for quarterly data)
 # n[scalar]: Integer, order of difference
 # Needs:
 # lmat[from here]
 #--------------------------------------------------------------------------
 if(T <= s*n) {
 warning("Time series are too short!")
 return(0*diag(T))}
 
 dd    <- diff(diag(T), lag = s, d = n)
 zeros <- matrix(0, nrow = n*s, ncol = T)
 dm    <- rbind(zeros, dd)
 return(dm)
 }
 #============================================================================
