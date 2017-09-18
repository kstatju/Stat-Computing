# Maximize x1,x2,x3
# 20x1 + 16x2 ???2x2 1 ???x2 2 ???x2 3,

# Maximize x1,x2,x3
# 2x2 1 +x2 2 +x2 3 - 20x1 - 16x2,
# subject to x1 + x2 ??? 5 x1 + x2 ???x3 = 0, x1 ??? 0, x2 ??? 0, x3 ??? 0.

require(quadprog)


k.matrix = 2*diag (c (2, 1, 1));
c.vector = c (-20, -16, 0);
a.matrix = matrix (c(1, 1, -1, -1, -1, 0, diag(1, 3, 3)), nrow=5, ncol=3, byrow = TRUE);
d.vector = c(0, -5, rep(0,3));
xHat = solve.QP(k.matrix, -c.vector, t(a.matrix), d.vector, meq = 1)
xHat






k.mat = 2*diag (c (1, 2));
k.mat
c.vecr = c (2, 8);
a.mat = matrix (c(-1, -2, -1, 0, 0, -1), nrow=3, ncol=2, byrow = TRUE);
a.mat
d.vec = c(10, rep(0,2));
k.lag = (a.mat %*% solve(k.mat) %*% t(a.mat))
k.lag
diag(k.lag) <- diag(k.lag)+10^-10
k.lag
c.lag = a.mat%*%solve(k.mat)%*%c.vec-d.vec
c.lag
a.mat = diag(1, 3, 3);
d.vec = c(rep(0,3))
xHat = solve.QP(k.lag, -c.lag, t(a.mat), d.vec, meq = 0)
xHat


