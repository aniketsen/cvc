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

#gamma_mu * gamma_5
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

#read loop
getloop <- function(filename) {
  loop <- array( NA, dim = c(12,vol,12) )
  read_loop <- as.numeric( readLines( filename ) )
  for ( x in 1:vol ) {
    for (i in 1:12 ) {
      for (j in 1:12) {
        loop[i,x,j] <- complex( real = read_loop[2*(((x-1)*12+i-1)*12+j-1)+1], imaginary = read_loop[2*(((x-1)*12+i-1)*12+j-1)+2] )
      }
    }
  }
  return(loop)
}

#gamma_mu * loop * gamma_mu
gloopg <- function(loop) {
  gloopg <- array(0,dim=c(12,vol,12))
  for ( mu in 1:4 ) {
    for ( x in 1:vol ) {
      for (alpha in 1:4 ) {
        for ( a in 1:3 ) {
          ka <- (alpha - 1)*3 + a
          for ( beta in 1:4 ) {
            for ( b in 1:3 ) {
              kb <- (beta - 1)*3 + b
              for ( gamma in 1:4 ) {
                kc <- (gamma - 1)*3 + a
                for ( delta in 1:4 ) {
                  kd <- (delta - 1)*3 + b
                  gloopg[ka,x,kb] <- gloopg[ka,x,kb] + gam[mu,alpha,gamma]*loop[kc,x,kd]*gam[mu,delta,beta]
                }
              }
            }
          }
        }
      }
    }
  }
  return(gloopg)
}

#gamma_5 * loop * gamma_5
g5loopg5 <- function(loop) {
  gloopg <- array(0,dim=c(12,vol,12))
  for ( x in 1:vol ) {
    for (alpha in 1:4 ) {
      for ( a in 1:3 ) {
        ka <- (alpha - 1)*3 + a
        for ( beta in 1:4 ) {
          for ( b in 1:3 ) {
            kb <- (beta - 1)*3 + b
            for ( gamma in 1:4 ) {
              kc <- (gamma - 1)*3 + a
              for ( delta in 1:4 ) {
                kd <- (delta - 1)*3 + b
                gloopg[ka,x,kb] <- gloopg[ka,x,kb] + gam5[alpha,gamma]*loop[kc,x,kd]*gam5[delta,beta]
              }
            }
          }
        }
      }
    }
  }
  return(gloopg)
}

#gamma_mu * gamma_5 * loop * gamma_mu * gamma_5
g5gloopg5g <- function(loop) {
  gloopg <- array(0,dim=c(12,vol,12))
  for ( mu in 1:4 ) {
    for ( x in 1:vol ) {
      for (alpha in 1:4 ) {
        for ( a in 1:3 ) {
          ka <- (alpha - 1)*3 + a
          for ( beta in 1:4 ) {
            for ( b in 1:3 ) {
              kb <- (beta - 1)*3 + b
              for ( gamma in 1:4 ) {
                kc <- (gamma - 1)*3 + a
                for ( delta in 1:4 ) {
                  kd <- (delta - 1)*3 + b
                  gloopg[ka,x,kb] <- gloopg[ka,x,kb] + gam5mu[mu,alpha,gamma]*loop[kc,x,kd]*gam5mu[mu,delta,beta]
                }
              }
            }
          }
        }
      }
    }
  }
  return(gloopg)
}

#change flavor of loop using gamma_5
loop_chng_flv <- function(loop) {
  loopd <- array(0,dim=c(12,vol,12))
  for ( x in 1:vol ) {
    for (alpha in 1:4 ) {
      for ( a in 1:3 ) {
        ka <- (alpha - 1)*3 + a
        for ( beta in 1:4 ) {
          for ( b in 1:3 ) {
            kb <- (beta - 1)*3 + b
            for ( gamma in 1:4 ) {
              kc <- (gamma - 1)*3 + a
              for ( delta in 1:4 ) {
                kd <- (delta - 1)*3 + b
                loopd[ka,x,kb] <- loopd[ka,x,kb] + gam5[alpha,gamma]*Conj(loop[kd,x,kc])*gam5[delta,beta]
              }
            }
          }
        }
      }
    }
  }
  return(loopd)
}

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

#loop * propagator
loop_ti_prop <- function(loop,prop) {
  seq_source <- array( 0, dim = c( 12, vol, 12 ) )
  for ( x in 1:vol ) {
    for ( i in 1:12 ) {
      for (j in 1:12 ) {
        for ( k in 1:12 ) {
          seq_source[i,x,j] <- seq_source[i,x,j] + loop[i,x,k]*prop[j,x,k]
        }
      }
    }
  }
  return(seq_source)
}

#compares read and calculated values of sequential source // outputs values at points of inequality
#for side by side comparison
compare <- function(seq_source,seq_done) {
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
}