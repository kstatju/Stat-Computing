f3 <- function(x1, x2){
  if (((x1-1)^2+(x2-1)^2 <= 1) & ((x1-1)^2+(x2+1)^2 <= 1)){
    return(x1^2+x2^2)
    
  }else {return(NA)}
}

f1 <- function(x1, x2){
  if ((x1-1)^2+(x2-1)^2 <= 1){
    return(x1^2+x2^2)
    
  }else {return(NA)}
}

f2 <- function(x1, x2){
  if ( (x1-1)^2+(x2+1)^2 <= 1){
    return(x1^2+x2^2)
    
  }else {return(NA)}
}


f <- function(x1, x2){
   return(x1^2+x2^2)

}

x1r = seq(-2, 2, by = 0.1)
x2r = seq(-2, 2, by = 0.1)


for (i in 1:length(x1r)){
  y = NULL
    for (j in 1:length(x2r)){
    x = f(x1r[i], x2r[j])
    y = rbind(y, x)
  }
  if (i == 1) {z = y; dim(z)
  } else {dim (z); z = cbind(z, y)}
}

for (i in 1:length(x1r)){
  y = NULL
  for (j in 1:length(x2r)){
    x = f1(x1r[i], x2r[j])
    y = rbind(y, x)
  }
  if (i == 1) {z1 = y
  } else {z1 = cbind(z1, y)}
}
for (i in 1:length(x1r)){
  y = NULL
  for (j in 1:length(x2r)){
    x = f2(x1r[i], x2r[j])
    y = rbind(y, x)
  }
  if (i == 1) {z2 = y;
  } else { z2 = cbind(z2, y)}
}

for (i in 1:length(x1r)){
  y = NULL
  for (j in 1:length(x2r)){
    x = f3(x1r[i], x2r[j])
    y = rbind(y, x)
  }
  if (i == 1) {z3 = y;
  } else { z3 = cbind(z3, y)}
}

contour(x2r, x1r, as.matrix(z))
contour(x2r, x1r, as.matrix(z1), add = TRUE, col = 2)
contour(x2r, x1r, as.matrix(z2), add = TRUE, col = 3)
contour(x2r, x1r, as.matrix(z3), add = TRUE, col = 4)

