#initialize lattice values
LT <- 16
LX <- 8
LY <- 8
LZ <- 8
vol <- LT*LX*LY*LZ
nloop <- vol*12*12

#initialize gamma matrices
gam <- array( NA, dim = c(4,4,4))
gam[1,1,] <-  c(0,0,-1,0)
gam[1,2,] <- c(0,0,0,-1)
gam[1,3,] <- c(-1,0,0,0)
gam[1,4,] <- c(0,-1,0,0)
gam[2,1,] <- c(0,0,0,-1.i)
gam[2,2,] <- c(0,0,-1.i,0)
gam[2,3,] <- c(0,1.i,0,0)
gam[2,4,] <- c(1.i,0,0,0)
gam[3,1,] <- c(0,0,0,-1)
gam[3,2,] <- c(0,0,1,0)
gam[3,3,] <- c(0,1,0,0)
gam[3,4,] <- c(-1,0,0,0)
gam[4,1,] <- c(0,0,-1.i,0)
gam[4,2,] <- c(0,0,0,1.i)
gam[4,3,] <- c(1.i,0,0,0)
gam[4,4,] <- c(0,-1.i,0,0)

#gamma_5
gam5 <- array(0,dim=c(4,4))

for (a in 1:4) {
  for (b in 1:4) {
    for (c in 1:4) {
      for (d in 1:4) {
        for ( e in 1:4) {
          gam5[a,b] <- gam5[a,b] + gam[1,a,c]*gam[2,c,d]*gam[3,d,e]*gam[4,e,b]
        }
      }
    }
  }
}

#gamma_5 * gamma_mu
gam5mu <- array(0,dim=c(4,4,4))
for (mu in 1:4) {
  for (a in 1:4) {
    for (b in 1:4) {
      for (c in 1:4) {
        gam5mu[mu,a,b] <- gam5mu[mu,a,b] + gam[mu,a,b]*gam5[b,c]
      }
    }
  }
}

gamma <- array(NA,dim=c(10,4,4))
gamma[1,,] <- gam[1,,]
gamma[2,,] <- gam[2,,]
gamma[3,,] <- gam[3,,]
gamma[4,,] <- gam[4,,]
gamma[5,1,] <- c(1,0,0,0)
gamma[5,2,] <- c(0,1,0,0)
gamma[5,3,] <- c(0,0,1,0)
gamma[5,4,] <- c(0,0,0,1)
gamma[6,,] <- gam5
gamma[7,,] <- gam5mu[1,,]
gamma[8,,] <- gam5mu[2,,]
gamma[9,,] <- gam5mu[3,,]
gamma[10,,] <- gam5mu[4,,]

#read propagator
getprop <- function(file_prefix) {
  prop <- array( NA, dim = c( 12, vol, 12 ) )
  for ( i in 1:12 ) {
    filename <- paste(file_prefix,i-1,".inverted.ascii",sep="")
    read_prop <- readLines(filename)
    s <- array( NA, dim = c(vol,12) )
    for ( j in 1:vol ) {
      for ( k in 1:12 ) {
        l <- (j-1)*12+k+1*j
        eval(parse(text=read_prop[l]))
      }
    }
    prop[i,,] <- s
  }
  return(prop)
}

#read sequential source
getseq <- function(file_prefix) {
  seq_done <- array( NA, dim = c( 12, vol, 12 ) )
  for ( i in 1:12 ) {
    filename <- paste(file_prefix,i-1,".ascii",sep="")
    read_seq <- readLines(filename)
    s <- array( NA, dim = c(vol,12) )
    for ( j in 1:vol ) {
      for ( k in 1:12 ) {
        l <- (j-1)*12+k+1*j
        eval(parse(text=read_seq[l]))
      }
    }
    seq_done[i,,] <- s
  }
  return(seq_done)
}

#read stochastic noise
getnoise <- function(filename) {
  scfld <- as.numeric( readLines( filename ) )
  return(scfld)
}

#prop * gamma * noise // The id here is same as gamma_id in Marcus' code
seq_W <- function(prop,id) {
  scfld <- getnoise(paste("scalar_field_",id,".txt",sep=""))
  seq_source <- array(0,dim=c(12,vol,12))
  for (x in 1:vol) {
    for (isc in 1:12) {
      for (alpha in 1:4) {
        for (a in 1:3) {
          ka <- (alpha-1)*3+a
          for (beta in 1:4) {
            kb <- (beta-1)*3+a
            seq_source[kb,x,isc] <- seq_source[kb,x,isc] + prop[isc,x,ka]*gamma[id+1,alpha,beta]*scfld[x]
          }
        }
      }
    }
  }
  return(seq_source)
}

#compare read and calculated values of sequential source
compare <- function(seq_source,seq_done) {
  flag <- c()
  ind <- 0
  for ( i in 1:12 ) {
    for ( x in 1:vol ) {
      for ( j in 1:12 ) {
        t1 <- signif(seq_source[j,x,i], digits = 6 )
        t2 <- signif(seq_done[i,x,j], digits = 6 )
        if ( t1 != t2) {
          ind <- ind+1
          print(c(ind,t1,t2))
        }
      }
    }
  }
  return(flag)
}

compareabs <- function(seq_source,seq_done) {
  flag <- c()
  ind <- 0
  for ( i in 1:12 ) {
    for ( x in 1:vol ) {
      for ( j in 1:12 ) {
        t11 <- abs(signif(Re(seq_source[j,x,i]), digits = 6 ))
        t12 <- abs(signif(Im(seq_source[j,x,i]), digits = 6 ))
        t21 <- abs(signif(Re(seq_done[i,x,j]), digits = 6 ))
        t22 <- abs(signif(Im(seq_done[i,x,j]), digits = 6 ))
        if ( t11 != t21 | t12 != t22) {
          ind <- ind+1
          print(c(ind,t11+t12*i,t21+t22*i))
        }
      }
    }
  }
}