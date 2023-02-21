# FGL new unconnected = 1

  flsa.general_copula <-
    function(A,L,lam1,lam2,penalize.diagonal)
    {
      trueA = A
      if(is.matrix(A[[1]])) {p=dim(A[[1]])[1]}
      if(is.vector(A[[1]])) {p=length(A[[1]])}
      K = length(A)
      # results matrices:
      X = list()
      #for(k in 1:K) {X[[k]] = matrix(NA,p,p)}
      for(k in 1:K) {X[[k]] = A[[1]]*NA}
      if(is.matrix(A[[1]])) {fusions = array(FALSE,dim=c(K,K,p,p))}
      if(is.vector(A[[1]])) {fusions = array(FALSE,dim=c(K,K,p,1))}

      # get starting newc: list of matrices.  newc[[k]][i,j] gives the (2*ordermats[[k]]-K-1) adjustment for lam2 at that element.  Basically, how many rank higher vs how many rank lower?
      newc = list()
      for(k in 1:K)
      {
        others = setdiff(1:K,k)
        others.smaller.k = 1:(k-1)
        newc[[k]] = A[[1]]*0
        for(o in others)
        {
          # fix the errors that occur when working with complex numbers
          if(is.complex(A[[o]]))
          {
            A[[o]] <- matrix(as.numeric(A[[o]]),nrow(A[[o]]),ncol(A[[o]]))
          }
          if(is.complex(A[[k]]))
          {
            A[[k]] <- matrix(as.numeric(A[[k]]),nrow(A[[k]]),ncol(A[[k]]))
          }

          #print(A[[o]]) # op zoek naar complex values!
          #print(A[[k]]) # op zoek naar complex values!

          newc[[k]] = newc[[k]] + (A[[o]]-A[[k]]< -1e-4) - (A[[o]]-A[[k]]>1e-4)
        }

        #for(o in others) {newc[[k]] = newc[[k]] + (A[[o]]-A[[k]]< -1e-4) - (A[[o]]-A[[k]]>1e-4)}
      }

      ######### start the loop here:
      for(iter in 1:(K-1))
      {

        # create order matrices:
        ordermats = list()
        for(k in 1:K)
        {
          others = setdiff(1:K,k)
          others.smaller.k = 1:(k-1)
          ordermats[[k]] = A[[1]]*0
          for(o in others) {ordermats[[k]] = ordermats[[k]] + (A[[k]]-A[[o]]>1e-4)}
          # to deal with ties, also add a unit to ordermat[[k]] if a class with a lower k has a tie at element i,j:
          if(k>1)
          {
            for(o in others.smaller.k) {ordermats[[k]] = ordermats[[k]] + (abs(A[[o]]-A[[k]])<1e-4)}
          }
          ordermats[[k]] = ordermats[[k]] + 1
        }

        # create beta.g matrices, holding the solution to Holger's "unconstrained problem"
        #  (prending we're not constraining the order of the solution to match the order of the A matrices)
        betas.g = list()
        for(k in 1:K)
        {
          betas.g[[k]] = A[[k]] - lam2/L*newc[[k]]
        }

        # identify and fuse all elements for which the betas.g are out of order:
        new.ordermats = list()
        for(k in 1:K)
        {
          others = setdiff(1:K,k)
          others.smaller.k = 1:(k-1)
          new.ordermats[[k]] = A[[1]]*0

          ##### hier gaat het weer mis met complex values, dus we doen dezelfde fix
          for(o in others)
          {
            if(is.complex(betas.g[[o]]))
            {
              betas.g[[o]] <- matrix(as.numeric(betas.g[[o]]),nrow(betas.g[[o]]),ncol(betas.g[[o]]))
            }
            if(is.complex(betas.g[[k]]))
            {
              betas.g[[k]] <- matrix(as.numeric(betas.g[[k]]),nrow(betas.g[[k]]),ncol(betas.g[[k]]))
            }
            new.ordermats[[k]] = new.ordermats[[k]] + (betas.g[[k]]-betas.g[[o]]>1e-4)
          }
          # for(o in others) {new.ordermats[[k]] = new.ordermats[[k]] + (betas.g[[k]]-betas.g[[o]]>1e-4)}

          # to deal with ties, also add a unit to ordermat[[k]] if a class with a lower k has a tie at element i,j:
          if(k>1)
          {
            for(o in others.smaller.k) {new.ordermats[[k]] = new.ordermats[[k]] + (abs(betas.g[[o]]-betas.g[[k]])<1e-4)}
          }
          new.ordermats[[k]] = new.ordermats[[k]] + 1
        }

        # identify neighboring fusions:  "fusions": K x K x p x p array: K x K matrices, T/F for fusions
        for(k in 1:K){
          for(kp in 1:K){
            #given k,kp, declare a fusion when their ordermats entries are adjacent, and their new.ordermats entries have the opposite direction:
            fusions[k,kp,,] = fusions[k,kp,,]+
              ((ordermats[[k]]-1==ordermats[[kp]])&(new.ordermats[[k]]<new.ordermats[[kp]]))+
              ((ordermats[[k]]+1==ordermats[[kp]])&(new.ordermats[[k]]>new.ordermats[[kp]]))+
              (abs(A[[k]]-A[[kp]])<1e-4)
            #(existing fusions, neighboring fusions, and ties)
            fusions = (fusions>0)*1
          }}


        # now we've noted fusions between all entries with adjacent ordermats entries and reversed new.ordermats entries
        # next: extend fusions to non-adjecent entries: if a-b and b-c, then connect a-c:
        for(k in 1:K){
          for(kp in 1:K){
            others = setdiff(1:K,c(k,kp))
            for(o in others)
            {
              #identify elements in o which are fused with the same element in both k and kp, then add them to the list of k-kp fusions:
              bothfused = fusions[k,o,,] & fusions[kp,o,,]
              fusions[k,kp,,] = fusions[k,kp,,] | bothfused
            }
          }}

        # now recalculate A with the new fused entries:
        # to recalculate A, for each non-zero entry, identify the classes k with which it must be fused, and get their average:
        for(k in 1:K)
        {
          others = setdiff(1:K,k)
          #fusemean and denom: the mean value of all the trueA to be fused, and the number of values to be fused:
          fusemean = trueA[[k]]
          denom = A[[1]]*0+1
          for(o in others)
          {
            fusemean = fusemean+fusions[k,o,,]*trueA[[o]]  #add the values of the elements which must be fused to fusemean
            denom = denom+fusions[k,o,,]
          }
          # now redefine A[[k]]: unchanged from trueA if there's no fusion, and the mean of the fused elements when there is fusion:
          A[[k]] = fusemean/denom
        }

        #newc: list of matrices.  newc[[k]][i,j] gives the (2*ordermats[[k]]-K-1) adjustment for lam2 at that element.  Basically, how many rank higher vs how many rank lower?
        newc = list()
        for(k in 1:K)
        {
          others = setdiff(1:K,k)
          others.smaller.k = 1:(k-1)
          newc[[k]] = A[[1]]*0
          for(o in others)
          {
            # fix the errors that occur when working with complex numbers
            if(is.complex(A[[o]]))
            {
              A[[o]] <- matrix(as.numeric(A[[o]]),nrow(A[[o]]),ncol(A[[o]]))
            }
            if(is.complex(A[[k]]))
            {
              A[[k]] <- matrix(as.numeric(A[[k]]),nrow(A[[k]]),ncol(A[[k]]))
            }

            #print(A[[o]]) # op zoek naar complex values!
            #print(A[[k]]) # op zoek naar complex values!

            newc[[k]] = newc[[k]] + (A[[o]]-A[[k]]< -1e-4) - (A[[o]]-A[[k]]>1e-4)
          }

          #for(o in others) {newc[[k]] = newc[[k]] + (A[[o]]-A[[k]]< -1e-4) - (A[[o]]-A[[k]]>1e-4)}
        }

      } #end loop here

      # final version of betas.g:
      for(k in 1:K)
      {
        betas.g[[k]] = A[[k]] - lam2/L*newc[[k]]
        if(is.complex(betas.g[[k]]))
        {
          betas.g[[k]] <- matrix(as.numeric(betas.g[[k]]),nrow(betas.g[[k]]),ncol(betas.g[[k]]))
        }
      }
      # now soft-threshold the solution matices:
      for(k in 1:K)
      {
        X[[k]] = soft_copula(betas.g[[k]],lam=lam1/L,penalize.diagonal=penalize.diagonal)
      }
      return(X)
    }

  flsa2_copula <-
    function(A,L,lam1,lam2,penalize.diagonal)  #A is a list of 2 matrices from which we apply an L2 penalty to departures
    {
      # 1st apply fused penalty:
      S1 = abs(A[[1]]-A[[2]])<=2*lam2/L
      X1 = (A[[1]]+A[[2]])/2
      Y1 = X1

      S2 = (A[[1]] > A[[2]]+2*lam2/L)
      X2 = A[[1]] - lam2/L
      Y2 = A[[2]] + lam2/L

      S3 = (A[[2]] > A[[1]]+2*lam2/L)
      X3 = A[[1]] + lam2/L
      Y3 = A[[2]] - lam2/L

      X = soft_copula(a = S1*X1 + S2*X2 + S3*X3, lam = lam1/L, penalize.diagonal=penalize.diagonal)
      Y = soft_copula(a = S1*Y1 + S2*Y2 + S3*Y3, lam = lam1/L, penalize.diagonal=penalize.diagonal)

      return(list(X,Y))
    }


  soft_copula <-
    function(a,lam,penalize.diagonal){ # if last argument is FALSE, soft-threshold a matrix but don't penalize the diagonal
      out <- sign(a)*pmax(0, abs(a)-lam)
      if(!penalize.diagonal) diag(out) <- diag(a)
      return(out)
    }


  penalty.as.matrix_copula <-
    function(lambda,p,penalize.diagonal)
    {
      # for matrix penalties:  check dim and symmetry:
      if(is.matrix(lambda))
      {
        if(sum(lambda!= t(lambda))>0) {stop("error: penalty matrix is not symmetric")}
        if(sum(abs(dim(lambda)-p))!=0 ) {stop("error: penalty matrix has wrong dimension")}
      }
      # for scalar penalties: convert to matrix form:
      if(length(lambda)==1) {lambda=matrix(lambda,p,p)}
      # apply the penalize.diagonal argument:
      if(!penalize.diagonal) {diag(lambda)=0}
      return(lambda)
    }


  #
  # admm.iters.unconnected_copula = function(S, unconnected, lambda1,lambda2,rho=1,rho.increment=1,weights,maxiter = 1000,tol=1e-5)
  # {
  #   K = length(S)
  #   for(k in 1:K){S[[k]]=as.matrix(S[[k]])}
  #   p = dim(S[[1]])[2]
  #   n=weights
  #
  #
  #   ### DIT IS EEN VERANDERING!!! NORMAAL ZAT DIT STUKJE CODE ER NIET IN
  #   #for(k in 1:K){
  #   #  S[[k]] <- S[[k]][unconnected]
  #   #}
  #
  #   # They use the variances for the unconnected variable for S, whilst we obtain the correlation.
  #
  #   # initialize theta:
  #   theta = list()
  #   for(k in 1:K){theta[[k]] = 1/S[[k]]}
  #   # initialize Z:
  #   Z = list(); for(k in 1:K){Z[[k]] = rep(0,p)}
  #   # initialize W:
  #   W = list();	for(k in 1:K){W[[k]] = rep(0,p)}
  #
  #   # initialize lambdas:
  #   lam1 = lambda1
  #   lam2 = lambda2
  #
  #   # iterations:
  #   iter=0
  #   diff_value = 10
  #   while((iter==0) || (iter<maxiter && diff_value > tol))
  #   {
  #     # update theta to minimize -logdet(theta) + <S,theta> + rho/2||theta - Z + W ||^2_2:
  #     theta.prev = theta
  #     for(k in 1:K)
  #     {
  #       B = n[k]*S[[k]] - rho*(Z[[k]] - W[[k]])
  #       theta[[k]] = 1/(2*rho) * ( -B + sqrt(B^2 + 4*rho*n[k]) )
  #     }
  #
  #     # update Z:
  #     # define A matrices:
  #     A = list()
  #     for(k in 1:K){ A[[k]] = theta[[k]] + W[[k]] }
  #     # use flsa to minimize rho/2 ||Z-A||_F^2 + P(Z):
  #     if(K==2){Z = flsa2_copula(A,rho,lam1,lam2,penalize.diagonal=TRUE)}
  #     if(K>2){Z = flsa.general_copula(A,rho,lam1,lam2,penalize.diagonal=TRUE)}
  #
  #
  #     # update the dual variable W:
  #     for(k in 1:K){W[[k]] = W[[k]] + (theta[[k]]-Z[[k]])}
  #
  #     # bookkeeping:
  #     iter = iter+1
  #     diff_value = 0
  #     for(k in 1:K) {diff_value = diff_value + sum(abs(theta[[k]] - theta.prev[[k]])) / sum(abs(theta.prev[[k]]))}
  #     # increment rho by a constant factor:
  #     rho = rho * rho.increment
  #   }
  #   diff = 0; for(k in 1:K){diff = diff + sum(abs(theta[[k]]-Z[[k]]))}
  #   out = list(theta=theta,Z=Z,diff=diff,iters=iter)
  #   return(out)
  # }


  admm.iters_copula = function(S, lambda1,lambda2,rho=1,rho.increment=1,weights,penalize.diagonal,maxiter = 1000,tol=1e-5)
  {
    K = length(S)
    p = dim(S[[1]])[2]
    # We moeten de dimensie van S aanpassen, omdat je anders in het geval van unconnected nodes een situatie krijgt
    # Waarbij je aantal rows gelijk is aan het aantal rows van je originele covariance matrix
    # Dit komt omdat het algoritme is ingesteld op een data matrix X, waarbij het aantal rijen niet uitmaakt!
    # We willen alleen de connected rows
    #for(k in 1:K){
    #  S[[k]] = S[[k]][connect_TF,]
    #}
    n=weights

    ns = c(); for(k in 1:K){ns[k] = dim(S[[k]])[1]}

    # initialize theta:
    theta = list()
    for(k in 1:K){theta[[k]] = diag(1/diag(S[[k]]))}
    # initialize Z:
    Z = list(); for(k in 1:K){Z[[k]]=matrix(0,p,p)}
    # initialize W:
    W = list();	for(k in 1:K) {W[[k]] = matrix(0,p,p) }

    # initialize lambdas:  (shouldn't need to do this if the function is called from the main wrapper function, JGL)
    lam1 = penalty.as.matrix_copula(lambda1,p,penalize.diagonal=penalize.diagonal)
    lam2 = penalty.as.matrix_copula(lambda2,p,penalize.diagonal=TRUE)
    # iterations:
    iter=0
    diff_value = 10
    while((iter==0) || (iter<maxiter && diff_value > tol))
    {
      # reporting
      if(FALSE)
      {
        print(paste("iter=",iter))
        print(paste("crit=",crit(theta,S,n=rep(1,K),lam1,lam2,penalize.diagonal=penalize.diagonal)))
        print(paste("crit=",crit(Z,S,n=rep(1,K),lam1,lam2,penalize.diagonal=penalize.diagonal)))
      }

      # update theta:
      theta.prev = theta
      for(k in 1:K){
        edecomp = eigen(S[[k]] - rho*Z[[k]]/n[k] + rho*W[[k]]/n[k])
        D = edecomp$values
        V = edecomp$vectors
        D2 = n[k]/(2*rho) * ( -D + sqrt(D^2 + 4*rho/n[k]) )
        theta[[k]] = V %*% diag(D2) %*% t(V)
      }

      # update Z:
      # define A matrices:
      A = list()
      for(k in 1:K){ A[[k]] = theta[[k]] + W[[k]] }
      # use flsa to minimize rho/2 ||Z-A||_F^2 + P(Z):
      if(K==2){Z = flsa2_copula(A,rho,lam1,lam2,penalize.diagonal=TRUE)}
      if(K>2){Z = flsa.general_copula(A,rho,lam1,lam2,penalize.diagonal=TRUE)}  # the option to not penalize the diagonal is exercised when we initialize the lambda matrices


      # update the dual variable W:
      for(k in 1:K){W[[k]] = W[[k]] + (theta[[k]]-Z[[k]])}

      # bookkeeping:
      iter = iter+1
      diff_value = 0
      for(k in 1:K) {diff_value = diff_value + sum(abs(theta[[k]] - theta.prev[[k]])) / sum(abs(theta.prev[[k]]))}
      # increment rho by a constant factor:
      rho = rho*rho.increment
    }
    diff = 0; for(k in 1:K){diff = diff + sum(abs(theta[[k]]-Z[[k]]))}
    out = list(theta=theta,Z=Z,diff=diff,iters=iter)
    return(out)
  }

  fused_gl_copula <- function(S, Y, lambda1,lambda2,rho=1,weights="equal",penalize.diagonal=FALSE,maxiter=1000,tol=1e-5,warm=NULL,return.whole.theta=TRUE,truncate=1e-5)
  {
    ## initialize:
    p = dim(S[[1]])[2]
    K = length(S)
    n = rep(0,K)
    for(k in 1:K) {
      n[k] = dim(S[[k]])[1]
    }

    # assign feature names if none exist:
    if(length(dimnames(S[[1]])[[2]])==0)
    {
      for(k in 1:K)
      {
        dimnames(S[[k]])[[2]]=paste("V",1:p,sep="")
      }
    }

    # set weights:
    weights = rep(1,K)

    ## Identify connected nodes
    connected = rep(TRUE,p)

    lam1 = lambda1
    lam2 = lambda2

    ## examine criteria:  (value 0 where S allows theta=0 to satisfy KKT; value 1 where theta must be connected)
    if(K==2)  #use bi-conditional screening rule to identify block structure exactly
    {
      crit1 = list()
      for(k in 1:K) { crit1[[k]] =  abs(S[[k]])*weights[k] > lam1 + lam2 }
      S.sum = matrix(0,sum(connected),sum(connected))
      for(k in 1:K) {S.sum = S.sum + weights[k]*S[[k]]}
      S.sum = abs(S.sum)
      crit2 = S.sum > 2*lam1
    }

    if(K>2)  #use sufficient screening rule to identify larger-grained block structure
    {
      crit1 = list()
      for(k in 1:K) { crit1[[k]] =  abs(S[[k]])*weights[k] > lam1 }
      crit2 = matrix(0,sum(connected),sum(connected))
    }

    # are both criteria met?
    critboth = crit2
    for(k in 1:K) {critboth = critboth + crit1[[k]]}
    critboth = (critboth!=0)
    diag(critboth) = 1



    ## now identify block structure using igraph:
    g1 <- graph.adjacency(critboth)
    cout = clusters(g1)
    blocklist = list()
    # identify unconnected elements, and get blocks:
    unconnected = c()

    # adapt cout$membership to start with index 1:
    if(min(cout$membership)==0){cout$membership=cout$membership+1}
    for(i in 1:(cout$no))
    {
      if(sum(cout$membership==i)==1) { unconnected <- c(unconnected,which(cout$membership==i)) }
      if(sum(cout$membership==i)>1) { blocklist[[length(blocklist)+1]] <- which(cout$membership==i) }
    }

    # final set of connected nodes
    connected[unconnected] = FALSE
    # connected indices of connected nodes:  0 for unconnected nodes, and 1:length(connected) for the rest.
    # maps features 1:p to their order in the connected features
    connected.index = rep(0,p)
    connected.index[connected] = 1:sum(connected)
    # regular indices of connected nodes: map connected nodes onto 1:p indexing:

    # redefine unconnected as !connected (up until now it's been extra nodes caught as unconnected)
    unconnected=!connected

    ## define theta on all connected:   (so theta is really theta.connected).
    theta = list()
    for(k in 1:K)
    {
      theta[[k]] = matrix(0,sum(connected),sum(connected))
      if(sum(connected)>0)
      {
        dimnames(theta[[k]])[[1]]=dimnames(theta[[k]])[[2]]=dimnames(Y[[k]])[[2]][connected]
      }
    }

    ## get solution on unconnected nodes
    # data:
    Su = list()
    for(k in 1:K){Su[[k]] = S[[k]][,unconnected]}
    # penalty vectors:
    # note: for admm.iters.unconnected, we use the penalize.diagonal argument before calling the function.  for admm.iters, we use it IN the function.
    lam1.unconnected = lambda1
    lam2.unconnected = lambda2

    # if penalize.diagonal==FALSE, then set the appropriate penalty vectors to zero:
    if(!penalize.diagonal)
    {
      lam1.unconnected = lam1.unconnected * 0
    }
    # get the unconnected portion of theta:

    theta.unconnected = vector("list", K)
    #for(k in 1:K){
    #  theta.unconnected[[k]] <- diag(p)
    #  theta.unconnected[[k]] <- theta.unconnected[[k]][,unconnected]
    #}
    for(k in 1:K){
      theta.unconnected[[k]] <- rep(1 , sum(unconnected))
    }
    for(k in 1:K) { names(theta.unconnected[[k]])=dimnames(S[[k]])[[2]][!connected] }
    if(sum(unconnected)==0) {theta.unconnected = NULL}

    ## now run JGL on each block of the connected nodes to fill in theta:
    if(length(blocklist)>0){
      for(i in 1:length(blocklist)){
        # the variables in the block
        bl <- blocklist[[i]]
        Sbl = list()
        # get the data on only those variables
        for(k in 1:K)
        {
          Sbl[[k]] = S[[k]][bl,bl]
        }
        # penalty matrices:
        lam1.bl = lambda1
        lam2.bl = lambda2
        # initialize lambdas:
        lam1.bl = penalty.as.matrix_copula(lam1.bl,dim(Sbl[[1]])[2],penalize.diagonal=penalize.diagonal)
        lam2.bl = penalty.as.matrix_copula(lam2.bl,dim(Sbl[[1]])[2],penalize.diagonal=TRUE)

        # implement warm start if desired
        # run JGL on the block:
        Thetabl = admm.iters_copula(Sbl, lam1.bl,lam2.bl,rho=rho,weights=weights,penalize.diagonal=TRUE,maxiter=maxiter,tol=tol)

        # Hier gaat het mis: Error in S[[k]] - rho * Z[[k]]/n[k] : non-conformable arrays
        # Dit moeten we dus opzoeken!!!

        # update Theta with Thetabl's results:
        for(k in 1:K) {theta[[k]][connected.index[bl],connected.index[bl]] = Thetabl$Z[[k]]}
      }}

    # round very small theta entries down to zero:
    if(dim(theta[[1]])[1]>0)
    {
      for(k in 1:K)
      {
        rounddown = abs(theta[[k]])<truncate; diag(rounddown)=FALSE
        theta[[k]]=theta[[k]]*(1-rounddown)
      }}


    # return output: theta on connected nodes, diagonal theta on unconnected nodes, and the identities of the connected nodes
    if(!return.whole.theta)
    {
      out = list(theta=theta,theta.unconnected=theta.unconnected,connected=connected)
    }
    if(return.whole.theta)
    {
      whole.theta = list()
      for(k in 1:K)
      {
        whole.theta[[k]] = matrix(0,p,p)
        diag(whole.theta[[k]])[unconnected] = theta.unconnected[[k]]
        whole.theta[[k]][connected,connected] = theta[[k]]
        dimnames(whole.theta[[k]])[[1]] = dimnames(whole.theta[[k]])[[2]] = dimnames(Y[[k]])[[2]]
      }
      out = list(theta=whole.theta,connected=connected)
    }
    class(out)="jgl"
    return(out)
  }
