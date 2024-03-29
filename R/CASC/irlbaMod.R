
# a matrix
# ---------------------------------------------------------------------

if (!require(irlba)) {
  install.packages('irlba', dependencies = T)
  require(irlba)
}

irlbaMod <-
  function (A, nu=5, nv=5, adjust=3, aug=c("ritz","harm"), sigma=c("ls","ss"), 
            maxit=1000, m_b=20, reorth=2, tol=1e-6, V=NULL,
            matmul=NULL)
  {
    # ---------------------------------------------------------------------
    # Check input parameters
    # ---------------------------------------------------------------------
    eps <- .Machine$double.eps
    # Profiling option
    options(digits.secs=3)
    m <- nrow(A[[1]])
    n <- ncol(A[[1]])
    k <- max(nu,nv)
    interchange <- FALSE
    sigma = match.arg(sigma)
    aug   = match.arg(aug)
    # Interchange dimensions m,n so that dim(A'A) = min(m,n) when seeking
    # the smallest singular values. This avoids finding zero smallest 
    # singular values.
    if (n>m && sigma=="ss") {
      t <- m
      m <- n
      n <- t
      interchange <- TRUE
    }
    
    # Increase the number of desired signular values by 'adjust' to
    # help convergence. k is re-adjusted as vectors converge--this is
    # only an initial value;
    k_org <- k;
    k <- k + adjust;
    if (k<=0)  stop ("k must be positive")
    if (k>min(m,n)) stop ("k must be less than min(m,n)+adjust")
    if (m_b<=1) stop ("m_b must be greater than 1")
    if (tol<0) stop ("tol must be non-negative")
    if (maxit<=0) stop ("maxit must be positive")
    if (m_b>= min(n,m)) {
      m_b <- floor(min(n,m)-0.1)
      if (m_b-k-1<0) {
        adjust <- 0
        k <- m_b-1
      }
    }
    if (m_b-k-1<0) m_b <- ceiling(k+1+0.1)
    if (m_b>=min(m,n)) {
      m_b <- floor(min(m,n)-0.1)
      adjust <- 0
      k <- m_b - 1
    }
    if (tol<eps) tol <- eps
    
    # Allocate memory for W and F:
    W <- matrix(0.0,m,m_b) 
    F <- matrix(0.0,n,1)
    # If starting matrix V is not given then set V to be an
    # (n x 1) matrix of normally distributed random numbers.
    # In any case, allocate V appropriate to problem size:
    if (is.null(V)) {
      V <- matrix(0.0,n,m_b)
      V[,1] <- rnorm(n)
    }
    else {
      V <- cbind(V, matrix(0.0,n,m_b-ncol(V)))
    }
    
    
    # ---------------------------------------------------------------------
    # Initialize local variables
    # ---------------------------------------------------------------------
    
    B <- NULL                  # Bidiagonal matrix
    Bsz <- NULL                # Size of B
    eps23 <- eps^(2/3)         # Used for Smax/avoids using zero
    I <- NULL                  # Indexing
    J <- NULL                  # Indexing
    iter <- 1                  # Man loop iteration count
    mprod <- 0                 # Number of matrix-vector products
    R_F <- NULL                # 2-norm of residual vector F
    sqrteps <- sqrt(eps)       #
    Smax <- 1                  # Max value of all computed singular values of
    # B est. ||A||_2
    Smin <- NULL               # Min value of all computed singular values of
    # B est. cond(A)
    SVTol <- max(sqrteps,tol)  # Tolerance for singular vector convergence
    S_B <- NULL                # Singular values of B
    U_B <- NULL                # Left singular vectors of B
    V_B <- NULL                # Right singular vectors of B
    V_B_last <- NULL           # last row of modified V_B 
    S_B2 <- NULL               # S.V. of [B ||F||]
    U_B2 <- NULL               # 
    V_B2 <- NULL               #  
    
    # ---------------------------------------------------------------------
    # Basic functions
    # ---------------------------------------------------------------------
    
    # Euclidean norm
    norm2 <- function (x) return(as.numeric(sqrt(crossprod(x))))
    
    # Orthogonalize vectors Y against vectors X. Y and X must be R matrix
    # objects (they must have a dim attribute).
    # Note: this function unnecessarily copies the contents of Y
    orthog <- function (Y,X)
    {
      if (dim(X)[2] < dim(Y)[2]) dotY <- crossprod (X,Y)
      else dotY <- t (crossprod(Y,X))
      return (Y - X %*% dotY)
    }
    
    # Convergence tests
    # Input parameters
    # Bsz            Number of rows of the bidiagonal matrix B
    # tol
    # k_org
    # U_B            Left singular vectors of small matrix B
    # S_B            Singular values of B
    # V_B            Right singular vectors of B
    # residuals      
    # k
    # SVTol
    # Smax
    #
    # Output parameter list
    # converged      TRUE/FALSE
    # U_B            Left singular vectors of small matrix B
    # S_B            Singular values of B
    # V_B            Right singular vectors of B
    # k              Number of singular vectors returned 
    convtests <- function (Bsz, tol, k_org, U_B, S_B, V_B, 
                           residuals, k, SVTol, Smax)
    {
      Len_res <- sum(residuals[1:k_org] < tol*Smax)
      if (Len_res == k_org) {
        return (list(converged=TRUE, U_B=U_B[,1:k_org, drop=FALSE], 
                     S_B=S_B[1:k_org, drop=FALSE], V_B=V_B[,1:k_org, drop=FALSE], k=k) )
      } 
      #   Not converged yet...
      #   Adjust k to include more vectors as the number of vectors converge.
      Len_res <- sum(residuals[1:k_org] < SVTol*Smax)
      k <- max(k, k_org + Len_res)
      if (k > Bsz -3) k <- Bsz -3
      return (list(converged=FALSE, U_B=U_B, S_B=S_B, V_B=V_B, k=k) )
    }
    
    # ---------------------------------------------------------------------
    # Main iteration
    # ---------------------------------------------------------------------
    
    while (iter <= maxit) {
      
      # ---------------------------------------------------------------------
      # Lanczos bidiagonalization iteration
      # Compute the Lanczos bidiagonal decomposition:
      # AV  = WB
      # t(A)W = VB + Ft(E)
      # with full reorthogonalization.
      # This routine updates W,V,F,B,mprod
      # ---------------------------------------------------------------------
      j <- 1
      #   Normalize starting vector:
      if (iter==1) V[,1] <- V[,1, drop=FALSE]/norm2(V[,1, drop=FALSE]) 
      else j <- k + 1
      
      #   Compute W=AV (the use of as.matrix here converts Matrix class objects)
      if(!is.null(matmul)) {
        #     User-specified matrix multiply function
        #print(matmul(A, V[,j,drop=FALSE], transpose=TRUE))
        if(interchange)
          W[,j] <- matmul(A, V[,j,drop=FALSE], transpose=TRUE)
        else
          W[,j] <- matmul(A, V[,j,drop=FALSE])
      }
      else {
        if (interchange)  W[,j] <- t (as.matrix(crossprod (V[,j,drop=FALSE], A)))
        else              W[,j] <- as.matrix(A %*% V[,j, drop=FALSE])
      }
      mprod <- mprod + 1
      
      #   Orthogonalize
      if (iter != 1) {
        W[,j] <- orthog (W[,j, drop=FALSE], W[,1:(j-1), drop=FALSE])
      }
      
      S <- norm2(W[,j, drop=FALSE])
      #   Check for linearly dependent vectors
      if ((S < SVTol) && (j==1)) stop ("Starting vector near the null space")
      if (S < SVTol) {
        W[,j] <- rnorm(nrow(W))
        W[,j] <- orthog(W[,j, drop=FALSE],W[,1:(j-1), drop=FALSE])
        W[,j] <- W[,j, drop=FALSE]/norm2(W[,j, drop=FALSE])
        S <- 0 
      }
      else W[,j] <- W[,j, drop=FALSE]/S
      
      #   Lanczos process
      while (j <= m_b) {
        if(!is.null(matmul)) {
          #       User-specified matrix multiply function
          if(interchange)
            F <- matmul(A, W[,j,drop=FALSE])
          else
            F <- matmul(A, W[,j,drop=FALSE], transpose=TRUE)
        }
        else{
          if (interchange) F <- as.matrix(A %*% W[,j, drop=FALSE])
          else F <- t(as.matrix(crossprod(W[,j,drop=FALSE],A)))
        }
        
        mprod <- mprod + 1
        F <- F - S*V[,j, drop=FALSE]
        #     Orthogonalize
        F <- orthog(F,V[,1:j, drop=FALSE])
        
        if (j+1 <= m_b) {
          R <- norm2(F)
          #       Check for linear dependence
          if (R<=SVTol) {
            F <- matrix(rnorm(dim(V)[1]),dim(V)[1],1)
            F <- orthog(F, V[,1:j, drop=FALSE])
            V[,j+1] <- F/norm2(F)
            R <- 0 
          }
          else V[,j+1] <- F/R
          
          #       Compute block diagonal matrix 
          if (is.null(B)) B <- cbind(S, R)
          else            B <- rbind(cbind(B,0),c(rep(0,j-1),S,R))
          
          if(!is.null(matmul)) {
            #         User-specified matrix multiply function
            if(interchange)
              W[,j+1] <- matmul(A, V[,j+1,drop=FALSE], transpose=TRUE)
            else
              W[,j+1] <- matmul(A, V[,j+1,drop=FALSE])
          }
          else{
            if (interchange) 
              W[,j+1] <- t (as.matrix(crossprod (V[,j+1, drop=FALSE],A)))
            else             W[,j+1] <- as.matrix(A %*% V[,j+1, drop=FALSE])
          }
          mprod <- mprod + 1
          
          #       One step of the classical Gram-Schmidt process
          W[,j+1] <- W[,j+1, drop=FALSE] - W[,j, drop=FALSE]*R
          
          #       Full reorthogonalization of W
          if (iter==1 || reorth==2)
            W[,j+1] <- orthog(W[,j+1, drop=FALSE],W[,1:j, drop=FALSE])
          S <- norm2(W[,j+1, drop=FALSE])
          if (S<=SVTol) {
            W[,j+1] <- rnorm(nrow(W))
            W[,j+1] <- orthog(W[,j+1, drop=FALSE],W[,1:j, drop=FALSE])
            W[,j+1] <- W[,j+1, drop=FALSE]/norm2(W[,j+1, drop=FALSE])
            S <- 0
          }
          else W[,j+1] <- W[,j+1, drop=FALSE]/S
        }
        else {
          #       Add a last block to matrix B
          B <- rbind(B,c(rep(0,j-1),S))
        }
        j <- j + 1
      }
      #cat ("iter = ",iter," j = ",j-1, "mprod = ",mprod,"\n",file=stderr())
      # ---------------------------------------------------------------------
      # (End of the Lanczos bidiagonalization part)
      # ---------------------------------------------------------------------
      
      Bsz <- nrow(B)
      R_F <- norm2(F)
      F <- F/R_F
      #   Compute singular triplets of B. Expect svd to return s.v.s in order
      #   from largest to smallest.
      Bsvd <- svd(B)
      
      #   Estimate ||A|| using the largest singular value over all iterations
      #   and estimate the cond(A) using approximations to the largest and 
      #   smallest singular values. If a small singular value is less than sqrteps
      #   use only Ritz vectors to augment and require two-sided reorthogonalization.
      if (iter ==1) {
        Smax <- Bsvd$d[1]
        Smin <- Bsvd$d[Bsz]
      }
      else {
        Smax <- max(Smax, Bsvd$d[1])
        Smin <- min(Smin, Bsvd$d[Bsz])
      }
      Smax <- max(eps23,Smax)
      if ((Smin/Smax < sqrteps) && reorth <2) {
        warning ("The matrix is ill-conditioned. Each basis will be reorthogonalized.")
        reorth <- 2
        aug <- "ritz"
      }
      
      #   Re-order the singular values accordingly.
      if (sigma=="ss") {
        jj <- seq (ncol (Bsvd$u), 1, by=-1)
        Bsvd$u <- Bsvd$u[,jj]
        Bsvd$d <- Bsvd$d[jj]
        Bsvd$v <- Bsvd$v[,jj]
      }
      
      #   Compute the residuals
      R <- R_F * Bsvd$u[Bsz,, drop=FALSE]
      
      #   Check for convergence
      ct <- convtests(Bsz, tol, k_org, Bsvd$u, Bsvd$d, Bsvd$v, abs(R), k, SVTol, Smax)
      k <- ct$k
      
      #   If all desired singular values converged, then exit main loop
      if (ct$converged) break
      
      if (iter>=maxit) break
      
      #   Compute the starting vectors and first block of B[1:k,1:(k+1), drop=FALSE]
      if (aug=="harm") {
        #     Update the SVD of B to be the svd of [B ||F||E_m]
        Bsvd2.d <- Bsvd$d
        Bsvd2.d <- diag(Bsvd2.d, nrow=length(Bsvd2.d))
        Bsvd2 <- svd (cbind (Bsvd2.d, t (R)))
        if (sigma=="ss") {
          jj <- seq (ncol (Bsvd2$u), 1, by=-1)
          Bsvd2$u <- Bsvd2$u[,jj]
          Bsvd2$d <- Bsvd2$d[jj]
          Bsvd2$v <- Bsvd2$v[,jj]
        }
        Bsvd$d <- Bsvd2$d
        Bsvd$u <- Bsvd$u %*% Bsvd2$u
        Bsvd$v <- cbind (rbind (Bsvd$v, rep (0,Bsz)), 
                         c (rep (0,Bsz), 1)) %*% Bsvd2$v
        V_B_last <- Bsvd$v [Bsz + 1, 1:k, drop=FALSE]
        s <- R_F * solve (B, cbind (c (rep(0,Bsz-1), 1)))
        Bsvd$v <- Bsvd$v[1:Bsz, , drop=FALSE] + s %*% Bsvd$v[Bsz+1, ,drop=FALSE]
        
        qrv <- qr (cbind ( rbind (Bsvd$v[,1:k], 0), rbind (-s, 1)))
        Bsvd$v <- qr.Q(qrv)
        R <- qr.R(qrv)
        V[,1:(k+1)] <- cbind(V, F) %*% Bsvd$v
        
        #     Update and compute the k x k+1 part of B
        UT <- t(R[1:(k+1), 1:k, drop=FALSE] + R[,k+1,drop=FALSE] %*% V_B_last)
        B <- diag(Bsvd$d[1:k],nrow=k) %*% (UT*upper.tri(UT,diag=TRUE))
      }
      else {
        #     Use the Ritz vectors
        V[,1:(k + dim(F)[2])] <- cbind(V[,1:(dim(Bsvd$v)[1]), drop=FALSE] %*% Bsvd$v[,1:k, drop=FALSE], F)
        B <- cbind( diag(Bsvd$d[1:k],nrow=k), R[1:k, drop=FALSE])
      }
      
      #   Update the left approximate singular vectors
      W[,1:k] <- W[, 1:(dim(Bsvd$u)[1]), drop=FALSE] %*% Bsvd$u[,1:k, drop=FALSE]
      
      iter <- iter + 1
    }
    # ---------------------------------------------------------------------
    # End of the main iteration loop
    # Output results
    # ---------------------------------------------------------------------
    d <- Bsvd$d[1:k_org]
    u <- W[, 1:(dim(Bsvd$u)[1]), drop=FALSE] %*% Bsvd$u[,1:k_org, drop=FALSE]
    v <- V[,1:(dim(Bsvd$v)[1]), drop=FALSE] %*% Bsvd$v[,1:k_org, drop=FALSE]
    # Adjust ordering if smallest singular values selected so that singular values
    # are reported in non-increasing order.
    if(sigma=="ss") {
      reverse <- seq(length(d),1)
      d <- d[reverse]
      u <- u[,reverse,drop=FALSE]
      v <- v[,reverse,drop=FALSE]
    }
    return (list(d=d, u=u[,1:nu,drop=FALSE], v=v[,1:nv,drop=FALSE], iter=iter,mprod=mprod))
  }