################################################################################

library(msm)      # dpexp, ppexp
library(pracma)   # hessian

################################################################################
# Distribucion power piecewise exponential (PPE)
################################################################################
dPPE <- function(y, rate = 1, alpha = 1, t = 0, log = FALSE) {
  if (any(rate <= 0)) stop("rate must be greater than 0")
  if (any(alpha <= 0)) stop("alpha must be greater than 0")
  if (!isTRUE(all.equal(0, t[1]))) stop("first element of t should be 0")
  if (is.unsorted(t)) stop("t should be in increasing order")
  if (any(y <= 0)) stop("y must be greater than 0")
  logfy <- log(alpha) + dpexp(y, rate = rate, t = t, log = TRUE) +
    (alpha - 1) * ppexp(y, rate = rate, t = t, log = TRUE)
  if (log) logfy else exp(logfy)
}

pPPE <- function(q, rate = 1, alpha = 1, t = 0, log = FALSE) {
  if (any(rate <= 0)) stop("rate must be greater than 0")
  if (any(alpha <= 0)) stop("alpha must be greater than 0")
  if (!isTRUE(all.equal(0, t[1]))) stop("first element of t should be 0")
  if (is.unsorted(t)) stop("t should be in increasing order")
  if (any(q <= 0)) stop("q must be greater than 0")
  logFy <- alpha * ppexp(q, rate = rate, t = t, log = TRUE)
  if (log) logFy else exp(logFy)
}

################################################################################
# Log-verosimilitudes (escala interna: beta, log lambda, log alpha, log(q+1))
# nu = probabilidad de cura = exp(z'beta) / (1 + exp(z'beta))
################################################################################
llike.PoPE <- function(psi, data, t) {
  y <- data$y; delta <- data$delta
  s <- nrow(data$z); L <- length(t); z <- matrix(data$z, nrow = s)
  beta <- matrix(psi[1:s], ncol = 1)
  lambda <- exp(psi[(s + 1):(s + L)])
  nu <- exp(t(z) %*% beta - log1p(exp(t(z) %*% beta)))
  logvf <- dpexp(y, rate = lambda, t = t, log = TRUE)
  vF <- ppexp(y, rate = lambda, t = t, log = FALSE)
  logfpop <- log(-log(nu)) + vF * log(nu) + logvf
  logSpop <- vF * log(nu)
  -sum(delta * logfpop + (1 - delta) * logSpop)
}

llike.BerPE <- function(psi, data, t) {
  y <- data$y; delta <- data$delta
  s <- nrow(data$z); L <- length(t); z <- matrix(data$z, nrow = s)
  beta <- matrix(psi[1:s], ncol = 1)
  lambda <- exp(psi[(s + 1):(s + L)])
  nu <- exp(t(z) %*% beta - log1p(exp(t(z) %*% beta)))
  logvf <- dpexp(y, rate = lambda, t = t, log = TRUE)
  vF <- ppexp(y, rate = lambda, t = t, log = FALSE)
  logfpop <- log(1 - nu) + logvf
  logSpop <- log(nu + (1 - nu) * (1 - vF))
  -sum(delta * logfpop + (1 - delta) * logSpop)
}

llike.NBPE <- function(psi, data, t) {
  y <- data$y; delta <- data$delta
  s <- nrow(data$z); L <- length(t); z <- matrix(data$z, nrow = s)
  beta <- matrix(psi[1:s], ncol = 1)
  lambda <- exp(psi[(s + 1):(s + L)])
  tau <- exp(psi[s + L + 1]) - 1
  nu <- exp(t(z) %*% beta - log1p(exp(t(z) %*% beta)))
  logvf <- dpexp(y, rate = lambda, t = t, log = TRUE)
  vF <- ppexp(y, rate = lambda, t = t, log = FALSE)
  logfpop <- (-1 / tau - 1) * log1p((nu^(-tau) - 1) * vF) + log((nu^(-tau) - 1) / tau) + logvf
  logSpop <- (-1 / tau) * log1p((nu^(-tau) - 1) * vF)
  -sum(delta * logfpop + (1 - delta) * logSpop)
}

llike.PoPPE <- function(psi, data, t) {
  y <- data$y; delta <- data$delta
  s <- nrow(data$z); L <- length(t); z <- matrix(data$z, nrow = s)
  beta <- matrix(psi[1:s], ncol = 1)
  lambda <- exp(psi[(s + 1):(s + L)])
  alpha <- exp(psi[s + L + 1])
  nu <- exp(t(z) %*% beta - log1p(exp(t(z) %*% beta)))
  logvf <- dPPE(y, rate = lambda, alpha = alpha, t = t, log = TRUE)
  vF <- pPPE(y, rate = lambda, alpha = alpha, t = t, log = FALSE)
  logfpop <- log(-log(nu)) + vF * log(nu) + logvf
  logSpop <- vF * log(nu)
  -sum(delta * logfpop + (1 - delta) * logSpop)
}

llike.BerPPE <- function(psi, data, t) {
  y <- data$y; delta <- data$delta
  s <- nrow(data$z); L <- length(t); z <- matrix(data$z, nrow = s)
  beta <- matrix(psi[1:s], ncol = 1)
  lambda <- exp(psi[(s + 1):(s + L)])
  alpha <- exp(psi[s + L + 1])
  nu <- exp(t(z) %*% beta - log1p(exp(t(z) %*% beta)))
  logvf <- dPPE(y, rate = lambda, alpha = alpha, t = t, log = TRUE)
  vF <- pPPE(y, rate = lambda, alpha = alpha, t = t, log = FALSE)
  logfpop <- log(1 - nu) + logvf
  logSpop <- log(nu + (1 - nu) * (1 - vF))
  -sum(delta * logfpop + (1 - delta) * logSpop)
}

llike.NBPPE <- function(psi, data, t) {
  y <- data$y; delta <- data$delta
  s <- nrow(data$z); L <- length(t); z <- matrix(data$z, nrow = s)
  beta <- matrix(psi[1:s], ncol = 1)
  lambda <- exp(psi[(s + 1):(s + L)])
  alpha <- exp(psi[s + L + 1])
  tau <- exp(psi[s + L + 2]) - 1
  nu <- exp(t(z) %*% beta - log1p(exp(t(z) %*% beta)))
  logvf <- dPPE(y, rate = lambda, alpha = alpha, t = t, log = TRUE)
  vF <- pPPE(y, rate = lambda, alpha = alpha, t = t, log = FALSE)
  logfpop <- (-1 / tau - 1) * log1p((nu^(-tau) - 1) * vF) + log((nu^(-tau) - 1) / tau) + logvf
  logSpop <- (-1 / tau) * log1p((nu^(-tau) - 1) * vF)
  -sum(delta * logfpop + (1 - delta) * logSpop)
}

################################################################################
# Funciones auxiliares de estimacion
################################################################################

# Particiones por cuantiles de los tiempos de evento (igual que en la version original)
part.quantiles <- function(data, L) {
  if (L <= 1) return(0)
  c(0, as.numeric(quantile(data$y[data$delta == 1], (1:(L - 1)) / L)))
}

# Optimiza nll desde valores iniciales generados por ini.fun(); repite hasta
# obtener convergencia y un Hessiano invertible con varianzas positivas.
.fit.optim <- function(nll, ini.fun, data, t, max.tries = 200, better.than = Inf) {
  for (k in 1:max.tries) {
    ini <- ini.fun()
    r1 <- try(optim(ini, nll, data = data, t = t, control = list(maxit = 10000)), silent = TRUE)
    if (inherits(r1, "try-error") || r1$convergence != 0 || !is.finite(r1$value)) next
    # Refinamiento con BFGS
    r2 <- try(optim(r1$par, nll, data = data, t = t, method = "BFGS",
                    control = list(maxit = 10000, reltol = 1e-12)), silent = TRUE)
    if (!inherits(r2, "try-error") && is.finite(r2$value) && r2$value < r1$value) r1 <- r2
    if (r1$value >= better.than) next
    H <- try(hessian(nll, x0 = r1$par, data = data, t = t), silent = TRUE)
    if (inherits(H, "try-error") || any(!is.finite(H))) next
    V <- try(solve(H), silent = TRUE)
    if (inherits(V, "try-error")) next
    if (any(!is.finite(diag(V))) || min(diag(V)) <= 0) next
    r1$vcov <- V
    return(r1)
  }
  stop("No se logro convergencia con Hessiano valido")
}

# Construye la tabla de resultados en escala natural (metodo delta)
#   trans: vector con "id" (beta), "exp" (lambda, alpha) o "expm1" (q = exp(x) - 1)
.make.results <- function(fit, trans, names, data, t, n.par) {
  par <- fit$par; V <- fit$vcov; se.int <- sqrt(diag(V))
  est <- ifelse(trans == "id", par, ifelse(trans == "exp", exp(par), exp(par) - 1))
  deriv <- ifelse(trans == "id", 1, exp(par))
  se <- abs(deriv) * se.int
  z975 <- qnorm(0.975)
  lower <- ifelse(trans == "id", par - z975 * se.int,
                  ifelse(trans == "exp", exp(par - z975 * se.int), exp(par - z975 * se.int) - 1))
  upper <- ifelse(trans == "id", par + z975 * se.int,
                  ifelse(trans == "exp", exp(par + z975 * se.int), exp(par + z975 * se.int) - 1))
  pval <- ifelse(trans == "id", 2 * pnorm(-abs(par / se.int)), NA)
  res <- cbind(Estimate = est, s.e. = se, lower = lower, upper = upper, p.value = pval)
  rownames(res) <- names
  n <- length(data$y)
  list(results = res, a = t, log.likelihood = -fit$value,
       AIC = 2 * fit$value + 2 * n.par, BIC = 2 * fit$value + log(n) * n.par,
       convergence = fit$convergence, par = par, vcov = V)
}

.ini.base <- function(data, s, L) c(rep(0, s), rep(-log(mean(data$y)), L)) + rnorm(s + L, sd = 0.1)

################################################################################
# Modelos de cura con piecewise exponential (PE)
# L = numero de intervalos (particiones por cuantiles). Para particiones
# manuales usar las versiones *.est2(data, tiempo), con tiempo[1] = 0.
################################################################################

PoPE.est2 <- function(data, tiempo) {
  s <- nrow(data$z); t <- tiempo; L <- length(t)
  fit <- .fit.optim(llike.PoPE, function() .ini.base(data, s, L), data, t)
  .make.results(fit, c(rep("id", s), rep("exp", L)),
                c(paste("beta", 0:(s - 1)), paste("lambda", 1:L)), data, t, s + L)
}
PoPE.est <- function(L = 1, data) PoPE.est2(data, part.quantiles(data, L))

BerPE.est2 <- function(data, tiempo) {
  s <- nrow(data$z); t <- tiempo; L <- length(t)
  fit <- .fit.optim(llike.BerPE, function() .ini.base(data, s, L), data, t)
  .make.results(fit, c(rep("id", s), rep("exp", L)),
                c(paste("beta", 0:(s - 1)), paste("lambda", 1:L)), data, t, s + L)
}
BerPE.est <- function(L = 1, data) BerPE.est2(data, part.quantiles(data, L))

NBPE.est2 <- function(data, tiempo) {
  s <- nrow(data$z); t <- tiempo; L <- length(t)
  po <- .fit.optim(llike.PoPE, function() .ini.base(data, s, L), data, t)   # valores iniciales
  fit <- .fit.optim(llike.NBPE, function() c(po$par, 0) + rnorm(s + L + 1, sd = 0.1), data, t)
  .make.results(fit, c(rep("id", s), rep("exp", L), "expm1"),
                c(paste("beta", 0:(s - 1)), paste("lambda", 1:L), "q"), data, t, s + L + 1)
}
NBPE.est <- function(L = 1, data) NBPE.est2(data, part.quantiles(data, L))

################################################################################
# Modelos de cura con power piecewise exponential (PPE)
################################################################################

PoPPE.est2 <- function(data, tiempo) {
  s <- nrow(data$z); t <- tiempo; L <- length(t)
  po <- .fit.optim(llike.PoPE, function() .ini.base(data, s, L), data, t)
  fit <- .fit.optim(llike.PoPPE, function() c(po$par, 0) + rnorm(s + L + 1, sd = 0.1), data, t)
  .make.results(fit, c(rep("id", s), rep("exp", L), "exp"),
                c(paste("beta", 0:(s - 1)), paste("lambda", 1:L), "alpha"), data, t, s + L + 1)
}
PoPPE.est <- function(L = 1, data) PoPPE.est2(data, part.quantiles(data, L))

BerPPE.est2 <- function(data, tiempo) {
  s <- nrow(data$z); t <- tiempo; L <- length(t)
  be <- .fit.optim(llike.BerPE, function() .ini.base(data, s, L), data, t)
  fit <- .fit.optim(llike.BerPPE, function() c(be$par, 0) + rnorm(s + L + 1, sd = 0.1), data, t)
  .make.results(fit, c(rep("id", s), rep("exp", L), "exp"),
                c(paste("beta", 0:(s - 1)), paste("lambda", 1:L), "alpha"), data, t, s + L + 1)
}
BerPPE.est <- function(L = 1, data) BerPPE.est2(data, part.quantiles(data, L))

NBPPE.est2 <- function(data, tiempo) {
  s <- nrow(data$z); t <- tiempo; L <- length(t)
  po <- .fit.optim(llike.PoPE, function() .ini.base(data, s, L), data, t)
  fit <- .fit.optim(llike.NBPPE, function() c(po$par, 0, 0) + rnorm(s + L + 2, sd = 0.1), data, t)
  .make.results(fit, c(rep("id", s), rep("exp", L), "exp", "expm1"),
                c(paste("beta", 0:(s - 1)), paste("lambda", 1:L), "alpha", "q"), data, t, s + L + 2)
}
NBPPE.est <- function(L = 1, data) NBPPE.est2(data, part.quantiles(data, L))

################################################################################
# Probabilidad de cura, fraccion de cura y eventos por intervalo
################################################################################

# Probabilidad de cura de cada paciente: nu_i = logit^{-1}(z_i' beta)
cure.prob <- function(fit, data) {
  s <- nrow(data$z)
  as.numeric(plogis(t(data$z) %*% fit$par[1:s]))
}

# Fraccion de cura promedio con IC por simulacion desde la normal asintotica
cure.fraction <- function(fit, data, nsim = 5000, seed = 2021) {
  s <- nrow(data$z)
  set.seed(seed)
  B <- MASS::mvrnorm(nsim, fit$par[1:s], fit$vcov[1:s, 1:s])
  sims <- apply(B, 1, function(b) mean(plogis(t(data$z) %*% b)))
  c(estimate = mean(cure.prob(fit, data)),
    lower = unname(quantile(sims, 0.025)), upper = unname(quantile(sims, 0.975)))
}

# Numero de eventos y de pacientes en riesgo en cada intervalo de la particion
events.by.interval <- function(data, t) {
  brk <- c(t, Inf); L <- length(t)
  ev <- table(cut(data$y[data$delta == 1], brk, right = FALSE))
  risk <- sapply(1:L, function(l) sum(data$y >= brk[l]))
  data.frame(interval = names(ev), at.risk = risk, events = as.numeric(ev))
}

################################################################################
# Residuos (sin cambios; psi en escala natural = fit$results[, 1])
################################################################################

nr.NBPE <- function(psi, data, t, nrep = 5) {
  y <- data$y; delta <- data$delta; s <- nrow(data$z); L <- length(t)
  z <- matrix(data$z, nrow = s); beta <- matrix(psi[1:s], ncol = 1)
  lambda <- psi[(s + 1):(s + L)]; tau <- psi[s + L + 1]
  nu <- exp(t(z) %*% beta - log1p(exp(t(z) %*% beta)))
  vF <- ppexp(y, rate = lambda, t = t, log = FALSE)
  Spop <- exp((-1 / tau) * log1p((nu^(-tau) - 1) * vF))
  u <- delta * (1 - Spop) + (1 - delta) * runif(length(y), 1 - Spop)
  mq <- sapply(1:nrep, function(i) sort(qnorm(runif(u))))
  apply(mq, 1, median)
}

nr.PoPE <- function(psi, data, t, nrep = 5) {
  y <- data$y; delta <- data$delta; s <- nrow(data$z); L <- length(t)
  z <- matrix(data$z, nrow = s); beta <- matrix(psi[1:s], ncol = 1)
  lambda <- psi[(s + 1):(s + L)]
  nu <- exp(t(z) %*% beta - log1p(exp(t(z) %*% beta)))
  vF <- ppexp(y, rate = lambda, t = t, log = FALSE)
  Spop <- exp(vF * log(nu))
  u <- delta * (1 - Spop) + (1 - delta) * runif(length(y), 1 - Spop)
  mq <- sapply(1:nrep, function(i) sort(qnorm(runif(u))))
  apply(mq, 1, median)
}

nr.BerPE <- function(psi, data, t, nrep = 5) {
  y <- data$y; delta <- data$delta; s <- nrow(data$z); L <- length(t)
  z <- matrix(data$z, nrow = s); beta <- matrix(psi[1:s], ncol = 1)
  lambda <- psi[(s + 1):(s + L)]
  nu <- exp(t(z) %*% beta - log1p(exp(t(z) %*% beta)))
  vF <- ppexp(y, rate = lambda, t = t, log = FALSE)
  Spop <- nu + (1 - nu) * (1 - vF)
  u <- delta * (1 - Spop) + (1 - delta) * runif(length(y), 1 - Spop)
  mq <- sapply(1:nrep, function(i) sort(qnorm(runif(u))))
  apply(mq, 1, median)
}

cox.NBPE <- function(data, psi, t) {
  s <- nrow(data$z); L <- length(t); z <- matrix(data$z, nrow = s)
  beta <- matrix(psi[1:s], ncol = 1); lambda <- psi[(s + 1):(s + L)]; tau <- psi[s + L + 1]
  nu <- exp(t(z) %*% beta - log1p(exp(t(z) %*% beta)))
  vF <- ppexp(data$y, rate = lambda, t = t, log = FALSE)
  -((-1 / tau) * log1p((nu^(-tau) - 1) * vF))
}

cox.PoPE <- function(data, psi, t) {
  s <- nrow(data$z); L <- length(t); z <- matrix(data$z, nrow = s)
  beta <- matrix(psi[1:s], ncol = 1); lambda <- psi[(s + 1):(s + L)]
  nu <- exp(t(z) %*% beta - log1p(exp(t(z) %*% beta)))
  vF <- ppexp(data$y, rate = lambda, t = t, log = FALSE)
  -(vF * log(nu))
}

cox.BerPE <- function(data, psi, t) {
  s <- nrow(data$z); L <- length(t); z <- matrix(data$z, nrow = s)
  beta <- matrix(psi[1:s], ncol = 1); lambda <- psi[(s + 1):(s + L)]
  nu <- exp(t(z) %*% beta - log1p(exp(t(z) %*% beta)))
  vF <- ppexp(data$y, rate = lambda, t = t, log = FALSE)
  -log(nu + (1 - nu) * (1 - vF))
}

################################################################################
# Modelos de cura con distribucion Weibull (misma estructura que los modelos PE)
# F(t) = 1 - exp(-(t/lambda)^alpha); escala interna: beta, log alpha, log lambda
################################################################################
.weib <- function(psi, s) list(alpha = exp(psi[s + 1]), lambda = exp(psi[s + 2]))

llike.BerW <- function(psi, data, t = NULL) {
  y <- data$y; delta <- data$delta; s <- nrow(data$z); z <- matrix(data$z, nrow = s)
  w <- .weib(psi, s); nu <- plogis(drop(t(z) %*% psi[1:s]))
  logvf <- dweibull(y, w$alpha, w$lambda, log = TRUE); vF <- pweibull(y, w$alpha, w$lambda)
  -sum(delta * (log(1 - nu) + logvf) + (1 - delta) * log(nu + (1 - nu) * (1 - vF)))
}
llike.PoW <- function(psi, data, t = NULL) {
  y <- data$y; delta <- data$delta; s <- nrow(data$z); z <- matrix(data$z, nrow = s)
  w <- .weib(psi, s); nu <- plogis(drop(t(z) %*% psi[1:s]))
  logvf <- dweibull(y, w$alpha, w$lambda, log = TRUE); vF <- pweibull(y, w$alpha, w$lambda)
  -sum(delta * (log(-log(nu)) + vF * log(nu) + logvf) + (1 - delta) * vF * log(nu))
}
llike.NBW <- function(psi, data, t = NULL) {
  y <- data$y; delta <- data$delta; s <- nrow(data$z); z <- matrix(data$z, nrow = s)
  w <- .weib(psi, s); tau <- exp(psi[s + 3]) - 1; nu <- plogis(drop(t(z) %*% psi[1:s]))
  logvf <- dweibull(y, w$alpha, w$lambda, log = TRUE); vF <- pweibull(y, w$alpha, w$lambda)
  logfpop <- (-1 / tau - 1) * log1p((nu^(-tau) - 1) * vF) + log((nu^(-tau) - 1) / tau) + logvf
  logSpop <- (-1 / tau) * log1p((nu^(-tau) - 1) * vF)
  -sum(delta * logfpop + (1 - delta) * logSpop)
}
.ini.w <- function(data, s) c(rep(0, s), log(1), log(mean(data$y) * 3)) + rnorm(s + 2, sd = 0.1)

BerW.est <- function(data) {
  s <- nrow(data$z)
  fit <- .fit.optim(llike.BerW, function() .ini.w(data, s), data, 0)
  .make.results(fit, c(rep("id", s), "exp", "exp"), c(paste("beta", 0:(s - 1)), "alpha", "lambda"), data, 0, s + 2)
}
PoW.est <- function(data) {
  s <- nrow(data$z)
  fit <- .fit.optim(llike.PoW, function() .ini.w(data, s), data, 0)
  .make.results(fit, c(rep("id", s), "exp", "exp"), c(paste("beta", 0:(s - 1)), "alpha", "lambda"), data, 0, s + 2)
}
NBW.est <- function(data) {
  s <- nrow(data$z)
  po <- .fit.optim(llike.PoW, function() .ini.w(data, s), data, 0)
  fit <- .fit.optim(llike.NBW, function() c(po$par, 0) + rnorm(s + 3, sd = 0.1), data, 0)
  .make.results(fit, c(rep("id", s), "exp", "exp", "expm1"), c(paste("beta", 0:(s - 1)), "alpha", "lambda", "q"), data, 0, s + 3)
}
