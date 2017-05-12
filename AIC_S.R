## Robust AIC (based on Claeskens et al.)

require(robustbase)

AIC_S <- function(y, X, beta, scale, cc) {
    re  <- drop((y - X %*% beta) / scale)
    U2  <- Mchi(re, cc = cc, psi = 'bisquare', deriv = 2)
    U1  <- Mchi(re, cc = cc, psi = 'bisquare', deriv = 1)
    J   <- crossprod(X, U2 * X)
    K   <- crossprod(U1 * X)
    tr  <- sum(diag(solve(J, K)))
    gof <- 2 * length(y) * log(scale)
    return(list(AIC = gof + 2 * tr, traza = tr))
}
