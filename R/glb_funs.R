# code from psych Package:
# Revelle, W. (2018) psych: Procedures for Personality and Psychological Research,
# Northwestern University, Evanston, Illinois, USA, https://CRAN.R-project.org/package=psych Version = 1.8.4.
glb.algebraic2 <- function (Cov, LoBounds = NULL, UpBounds = NULL)
{
  if (!requireNamespace("Rcsdp")) {
    stop("Rcsdp must be installed to find the glb.algebraic")
  }
  cl <- match.call()
  p <- dim(Cov)[2]
  if (dim(Cov)[1] != p)
    Cov <- cov(Cov)
  if (any(t(Cov) != Cov))
    stop("'Cov' is not symmetric")
  if (is.null(LoBounds))
    LoBounds <- rep(0, ncol(Cov))
  if (is.null(UpBounds))
    UpBounds <- diag(Cov)
  if (any(LoBounds > UpBounds)) {
    stop("'LoBounds'<='UpBounds' violated")
  }
  if (length(LoBounds) != p)
    stop("length(LoBounds) != dim(Cov)")
  if (length(UpBounds) != p)
    stop("length(UpBounds)!=dim(Cov)")
  Var <- diag(Cov)
  opt = rep(1, p)
  C <- list(diag(Var) - Cov, -UpBounds, LoBounds)
  A <- vector("list", p)
  for (i in 1:p) {
    b <- rep(0, p)
    b[i] <- 1
    A[[i]] <- list(diag(b), -b, b)
  }
  K <- list(type = c("s", "l", "l"), size = rep(p, 3))
  result <- Rcsdp::csdp(C, A, opt, K, control = Rcsdp::csdp.control(printlevel = 0))
  if (result$status >= 4 || result$status == 2) {
    warning("Failure of csdp, status of solution=", result$status)
    lb <- list(glb = NA, solution = NA, status = result$status,
               Call = cl)
  }
  else {
    if (result$status != 0) {
      warning("status of solution=", result$status)
    }
    item.diag <- result$y
    names(item.diag) <- colnames(Cov)
    glb <- (sum(Cov) - sum(Var) + sum(result$y))/sum(Cov)
  }
  return(glb)
}




# adjusted code from psych Package to make it execute faster:
# Revelle, W. (2018) psych: Procedures for Personality and Psychological Research,
# Northwestern University, Evanston, Illinois, USA, https://CRAN.R-project.org/package=psych Version = 1.8.4.
glbOnArray <- function(Cov) {

  d <- dim(Cov)
  if (length(d) == 2L) { # turn it into an array if it is a matrix
    d <- c(1L, d)
    dim(Cov) <- d
  }

  nSamples <- d[1L]
  p <- d[2L]

  opt <- rep.int(1L, p)
  A <- vector("list", p)
  for (i in seq_len(p)) {
    b <- rep(0, p)
    b[i] <- 1
    A[[i]] <- list(diag(b), -b, b)
  }
  K <- list(type = c("s", "l", "l"), size = rep(p, 3))

  control <- Rcsdp::csdp.control(printlevel = 0)
  write.control.file2(control)
  on.exit(unlink("param.csdp"))

  prob.info <- get.prob.info2(K, length(b))
  LoBounds <- rep(0, p)

  cv <- Cov[1L, , ]
  Var <- diag(cv)
  C <- list(diag(Var) - cv, -Var, LoBounds)

  # make the Rcsdp object once instead of each iteration
  prob.data <- list(
    C = Rcsdp:::blkmatrix_R2csdp(C, prob.info),
    A = Rcsdp:::constraints_R2csdp(A, prob.info),
    b = as.double(c(0, opt))
  )

  glbs <- numeric(nSamples)
  for (i in seq_len(nSamples)) {

    cv <- Cov[i, , ]
    Var <- diag(cv)

    prob.data$C$blocks[[1L]]$data <- as.double(diag(Var) - cv)
    prob.data$C$blocks[[2L]]$data <- as.double(c(0, -Var))
    # The two lines above are equivalent to recreating the Rcsdp object:
    # C <- list(diag(Var) - cv, -Var, LoBounds)
    # prob.data0 <- Rcsdp:::prepare.data(C, A, opt, prob.info)

    ret <- .Call(
      "csdp",
      as.integer(sum(prob.info$block.sizes)),
      as.integer(prob.info$nconstraints),
      as.integer(prob.info$nblocks),
      as.integer(c(0, prob.info$block.types)),
      as.integer(c(0, prob.info$block.sizes)),
      prob.data$C,
      prob.data$A,
      prob.data$b,
      PACKAGE = "Rcsdp"
    )

    scv <- sum(cv)

    # ret[[3L]][-1L] is equivalent to:
    # ret[1:3] <- Rcsdp:::get.solution(ret[[1L]], ret[[2L]], ret[[3L]], prob.info)
    # result <- structure(ret, names = c("X", "Z", "y", "pobj", "dobj", "status"))
    # or
    # y <- vector_csdp2R(ret[[3L]])

    glbs[i] <- (scv - sum(Var) + sum(ret[[3L]][-1L])) / scv

  }
  return(glbs)
}

get.prob.info2 <- function(K, m) {
  # block.types <- ifelse(K$type == "s", 1, 2)
  block.types <- (K$type != "s") + 1L
  nblocks <- length(K$type)
  block.sizes <- K$size
  nconstraints <- m
  ret <- list(nblocks = nblocks, nconstraints = nconstraints, block.types = block.types, block.sizes = block.sizes)
  return(ret)
}

write.control.file2 <- function(control) {
  fileptr <- file("param.csdp", "w")
  for (i in 1:length(control)) cat(names(control)[i], "=", control[[i]], "\n", sep = "", file = fileptr)
  close(fileptr)
}
