tt <- seq(0, 4500)
rho1_t <- numeric(length(tt))
rho2_t <- numeric(length(tt))


for (t in seq_along(tt)) {
  
  vax_starts = c(3630-1144, 3630-1127, 3630-588, 3630-438, 3630-417, 3630-256, 3630-221, 3630-100, 3630-70, 3630+165, 3630+287, 3630+374, 3630+411, 3630+467, 3630+485)
  vax2D = c(TRUE, TRUE, TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE)
  vax_dur1 = c(3,4,8,6,7,4,9,7,8,7,9,10,7,10,6)
  pvax1 = c(72500, 412767, 1544464, 1167939, 840774, 101245, 1036869, 160340, 265188, 100713, 1909405, 2229941, 1858472, 1522407, 862326)
  vax_dur2 = c(3, 10, 6, 7, NA, 4, 6, NA, NA, NA, NA, NA, NA, NA, NA)
  vax_int = c(14, 168, 63, 84, NA, 14, 91, NA, NA, NA, NA, NA, NA, NA, NA)
  pvax2 = c(72500, 286203, 1515214, 1164639, NA, 101245, 1036859, NA, NA, NA, NA, NA, NA, NA, NA)
  
#vaccination <- function(t, vax_starts, vax2D, vax_dur1, pvax1, vax_dur2, vax_int, pvax2) {
  
  if (length(vax_starts) == 0 ) {
    return(c(rho1_t, rho2_t))
  } else {
    for (j in seq_along(vax_starts)) {
      rho1_jt <- numeric(length(tt))
      rho2_jt <- numeric(length(tt))
      
      rho1_jt[vax_starts[j] : (vax_starts[j]+vax_dur1[j])] <- 1
      
      rho1_j <- (-(log(1-(pvax1[j]/initP))) / vax_dur1[j])
      rho1_jt <- rho1_jt * rho1_j
      
      rho1_t <- rho1_t + rho1_jt
      
      if (vax2D[j] == TRUE) {

        rho2_jt[(vax_starts[j]+vax_dur1[j]+vax_int[j]): (vax_starts[j]+vax_dur1[j]+vax_int[j]+vax_dur2[j])] <- 1

        cov2 <- pvax2[j]/initP
        cov2_capped <- ifelse(cov2 >= 1, 0.99, cov2)
        
        rho2_j <- (-(log(1-cov2_capped)) / vax_dur2[j])
        rho2_jt <- rho2_jt * rho2_j
        
        rho2_t <- rho2_t + rho2_jt
        
      }
    }
  }
  rho1 <- rho1_t[t]
  rho2 <- rho2_t[t]
  
  #return(c(rho1_t,rho2_t))
#}

#vals <- vaccination(t, vax_starts, vax2D, vax_dur1, pvax1, vax_dur2, vax_int, pvax2)
#rho1 <- vals[1]
#rho2 <- vals[2]

}


plot(rho1_t)
plot(rho2_t)





