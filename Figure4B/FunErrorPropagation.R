# http://lectureonline.cl.msu.edu/~mmp/labs/error/e2.htm
# this is for x/y (not for x-y)
funErrorProp <- function(x, y, z, x.error, y.error, z.error){
  
  R1 <- x/x
  R2 <- y/x
  R3 <- z/x
  
  # absolute values of the ratios
  R1.abs <- abs(R1)
  R2.abs <- abs(R2)
  R3.abs <- abs(R3)
  
  r.x <- (x.error/x)^2
  r.y <- (y.error/y)^2
  r.z <- (z.error/z)^2
  
  new.error.1 <- R1.abs*(r.x + r.x)^(1/2)
  new.error.2 <- R2.abs*(r.x + r.y)^(1/2)
  new.error.3 <- R3.abs*(r.x + r.z)^(1/2)
  
  new.error <- data.frame(new.error.1, new.error.2, new.error.3)
    
  return(new.error)
  
}
