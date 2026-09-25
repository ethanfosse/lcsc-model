# 02_mortality_france.R -- French male mortality, 1816-2020
# Fosse and Winship, "Varieties of Cross-Cohort Differentiation" (Sociological Science).
# Input: Data/mortality.RData (not included; see README.md). Packages: mgcv, plot3D, RColorBrewer, ggplot2, ggridges, dplyr, tidyr.
# Run from this folder: Rscript --vanilla 02_mortality_france.R
# Writes the paper's figures to Figures/ and tables to Output/.

CRITERION_FIGURES <- "REML" # smoothing criterion of the figures model (paper: REML)
CRITERION_TABLES  <- "ML"   # smoothing criterion of the fit-table models (paper: ML)
LEGACY_UBRE       <- FALSE
NTHREADS          <- 4      # threads passed to mgcv::gam.control()
REFIT             <- FALSE

# war-shock ages 15-50, as in the paper's footnote
WAR_AGES <- 15:50
ROBUSTNESS_1941_42 <- TRUE  # also fit the 1941-1942 robustness model

#### 1. Preliminaries ####

# ---- BEGIN EMBEDDED MODELING HELPERS -------------------------------------------
if (!requireNamespace("mgcv", quietly = TRUE)) {
  stop("Embedded LC-SC helpers requires the 'mgcv' package.")
}
suppressPackageStartupMessages(library(mgcv))

# 1. CONFIGURATION

# lcsc_config(): the constants of an LC-SC analysis, checked against the data
lcsc_config <- function(data,
                        age = "a.index", period = "p.index", cohort = "c.index",
                        age_center = NULL, period_center = NULL, cohort_center = NULL,
                        period_full = NULL,
                        knots = c(a = 20, p = 20, c = 20),
                        quiet = FALSE) {

  for (v in c(age, period, cohort)) {
    if (!v %in% names(data)) {
      stop("lcsc_config(): column '", v, "' not found in the data. ",
           "The scripts expect natural-scale index columns named a.index, p.index ",
           "and c.index; pass age=/period=/cohort= if yours are named differently.")
    }
  }

  a.i <- as.numeric(data[[age]])
  p.i <- as.numeric(data[[period]])
  c.i <- as.numeric(data[[cohort]])

  if (any(is.na(a.i)) || any(is.na(p.i)) || any(is.na(c.i))) {
    stop("lcsc_config(): the age, period or cohort index columns contain NA. ",
         "Drop incomplete rows (na.omit) before building the configuration.")
  }

  # cohort must be period - age
  if (!isTRUE(all.equal(c.i, p.i - a.i))) {
    stop("lcsc_config(): cohort is not equal to period - age for every row. ",
         "The LC-SC model requires that identity; recompute the cohort column.")
  }

  derived <- c(a = round(mean(a.i)), p = round(mean(p.i)))
  age_center    <- if (is.null(age_center))    derived[["a"]] else age_center
  period_center <- if (is.null(period_center)) derived[["p"]] else period_center
  derived[["c"]] <- period_center - age_center
  if (is.null(cohort_center)) {
    cohort_center <- derived[["c"]]
  } else if (cohort_center != period_center - age_center) {
    stop("lcsc_config(): cohort_center must equal period_center - age_center (",
         period_center, " - ", age_center, " = ", period_center - age_center,
         "), so that c = p - a holds on the centered scale; got ", cohort_center, ".")
  }

  # The full period grid the spline basis is built on
  period_full <- if (is.null(period_full)) min(p.i):max(p.i) else sort(unique(period_full))
  missing_years <- setdiff(sort(unique(p.i)), period_full)
  if (length(missing_years)) {
    stop("lcsc_config(): period_full does not cover every observed year. ",
         "Missing: ", paste(utils::head(missing_years, 10), collapse = ", "),
         if (length(missing_years) > 10) ", ..." else "",
         "\n  period_full spans ", min(period_full), "-", max(period_full),
         " but the data span ", min(p.i), "-", max(p.i), ".")
  }

  if (!all(c("a", "p", "c") %in% names(knots))) {
    stop("lcsc_config(): `knots` must be a named vector, e.g. c(a = 20, p = 20, c = 20).")
  }
  n_unique <- c(a = length(unique(a.i)), p = length(period_full), c = length(unique(c.i)))
  for (d in c("a", "p", "c")) {
    if (knots[[d]] > n_unique[[d]] + 1) {
      warning("lcsc_config(): ", knots[[d]], " knots requested for '", d,
              "' but only ", n_unique[[d]], " distinct values are available. ",
              "mgcv will complain; reduce the knot count.")
    }
  }

  cfg <- list(
    age = age, period = period, cohort = cohort,
    age_center = age_center, period_center = period_center, cohort_center = cohort_center,
    period_full = period_full,
    grids = list(p = period_full - period_center),   # centered period grid for fit_lcsc()
    knots = knots,
    derived_centers = derived,
    n_unique = n_unique,
    ranges = list(age = range(a.i), period = range(p.i), cohort = range(c.i))
  )
  class(cfg) <- "lcsc_config"
  if (!quiet) print(cfg)
  cfg
}

print.lcsc_config <- function(x, ...) {
  cat("\n---- LC-SC configuration ---------------------------------------------\n")
  cat(sprintf("  age    : %4d-%-4d  (%d distinct)  centered at %d  [mean %d]\n",
              x$ranges$age[1], x$ranges$age[2], x$n_unique[["a"]],
              x$age_center, x$derived_centers[["a"]]))
  cat(sprintf("  period : %4d-%-4d  (%d distinct)  centered at %d  [mean %d]\n",
              x$ranges$period[1], x$ranges$period[2],
              length(unique(x$period_full)), x$period_center, x$derived_centers[["p"]]))
  cat(sprintf("  cohort : %4d-%-4d  (%d distinct)  centered at %d  [= period - age center]\n",
              x$ranges$cohort[1], x$ranges$cohort[2], x$n_unique[["c"]], x$cohort_center))
  cat(sprintf("  period spline built on the full grid %d:%d (%d years)\n",
              min(x$period_full), max(x$period_full), length(x$period_full)))
  cat(sprintf("  knots  : age %d, period %d, cohort %d\n",
              x$knots[["a"]], x$knots[["p"]], x$knots[["c"]]))
  if (any(abs(c(x$age_center, x$period_center) -
              c(x$derived_centers[["a"]], x$derived_centers[["p"]])) > 0)) {
    cat("  note   : a centering constant differs from round(mean(.)) of these data.\n")
    cat("           That is fine -- centering does not change the fit -- but it does\n")
    cat("           change what the intercept means. Deliberate for the published\n")
    cat("           analyses; check it is deliberate for yours.\n")
  }
  cat("----------------------------------------------------------------------\n\n")
  invisible(x)
}

# lcsc_center(): adds the centered modeling columns a, p, c to a data frame
lcsc_center <- function(data, cfg) {
  data$a <- as.numeric(data[[cfg$age]])    - cfg$age_center
  data$p <- as.numeric(data[[cfg$period]]) - cfg$period_center
  data$c <- as.numeric(data[[cfg$cohort]]) - cfg$cohort_center
  data
}

# 2. ORTHOGONALIZATION

# orthogonalize_basis(): the reparameterization at the heart of the model
orthogonalize_basis <- function(x_obs, k, grid = NULL, tol = 1e-8) {
  values <- if (is.null(grid)) sort(unique(x_obs)) else sort(unique(grid))
  if (!all(x_obs %in% values)) {
    bad <- sort(unique(x_obs[!x_obs %in% values]))
    stop("orthogonalize_basis(): observed value(s) ", paste(utils::head(bad, 10), collapse = ", "),
         if (length(bad) > 10) ", ..." else "", " are not in the grid.",
         "\n  For period this almost always means period_full does not cover the data;",
         " see lcsc_config().")
  }
  gdat <- data.frame(x = values)
  # unconstrained basis (k + 1 columns) on the grid
  sm <- mgcv::smoothCon(s(x, bs = "cr", k = k + 1), data = gdat, absorb.cons = FALSE)[[1]]
  Xuni <- sm$X
  smS  <- sm$S[[1]]
  # M = [1, x] on the grid
  M <- cbind(1, values)
  C <- t(M) %*% Xuni
  qrc <- qr(t(C))
  Z <- qr.Q(qrc, complete = TRUE)[, (nrow(C) + 1):ncol(C), drop = FALSE]
  curUni <- Xuni %*% Z                 # orthogonal to 1 and x over the grid
  S <- t(Z) %*% smS %*% Z              # reparameterized penalty (full rank k - 1)
  S <- (S + t(S)) / 2                  # enforce exact symmetry
  worst <- max(abs(t(M) %*% curUni))
  if (!is.finite(worst) || worst > tol) {
    stop("orthogonalize_basis(): the basis is not orthogonal to {intercept, linear}; ",
         "largest inner product ", format(worst, digits = 3),
         ". This should not happen -- check for duplicate or non-finite covariate values.")
  }
  cur <- curUni[match(x_obs, values), , drop = FALSE]
  list(cur = cur, curUni = curUni, values = values, Z = Z, S = S, sm = sm, k = k, orth_error = worst)
}

# basis_at(): evaluate an orthogonalized basis
basis_at <- function(ob, xnew) {
  Xraw <- mgcv::PredictMat(ob$sm, data.frame(x = xnew))
  Xraw %*% ob$Z
}

# orthogonalize_block(): residualize a block of design columns on other columns
orthogonalize_block <- function(H, X0, w = NULL) {
  if (is.null(w)) w <- rep(1, nrow(X0))
  if (length(w) != nrow(X0) || any(w < 0) || !all(is.finite(w))) stop("orthogonalize_block(): invalid weights")
  sw <- sqrt(w)
  qr0 <- qr(X0 * sw)
  if (qr0$rank < ncol(X0)) stop("orthogonalize_block(): the common columns are rank deficient")
  B <- qr.coef(qr0, H * sw)
  B[is.na(B)] <- 0
  Hp <- H - X0 %*% B
  stopifnot(max(abs(crossprod(X0 * w, Hp))) < 1e-6 * max(1, max(abs(H))) * max(1, max(w)))
  list(Hperp = Hp, B = B)
}

# 3. FITTING

# fit_lcsc(): fit an LC-SC GAM from an mgcv prefit object
fit_lcsc <- function(G, data, cfg = NULL, method = c("REML", "ML"), grids = NULL,
                     fs_exclude = NULL, fs_weights = NULL, legacy_ubre = FALSE,
                     nthreads = NULL, in.out = NULL, fs_project = TRUE) {
  method <- match.arg(method)
  if (is.null(G$X) || is.null(G$smooth)) stop("fit_lcsc(): G must be an unfitted mgcv model, gam(..., fit = FALSE)")
  if (nrow(G$X) != nrow(data)) stop("fit_lcsc(): `data` has ", nrow(data), " rows but the prefit has ", nrow(G$X))
  if (is.null(grids)) grids <- if (is.null(cfg)) list() else cfg$grids
  n.smooth <- length(G$smooth)
  if (n.smooth < 1) stop("fit_lcsc(): no smooth terms in the prefit")
  n.pen <- sapply(G$smooth, function(s) length(s$S))
  pen.index <- cumsum(n.pen) - n.pen + 1
  is.fs <- sapply(G$smooth, inherits, "fs.interaction")
  ortho <- vector("list", n.smooth)

  ## ---- (i) univariate cubic regression splines: orthogonalized to {1, x} on the grid --- #
  for (i in which(!is.fs)) {
    s <- G$smooth[[i]]
    if (length(s$term) != 1) stop("fit_lcsc(): smooth ", s$label, " is not univariate")
    if (!inherits(s, "cr.smooth")) stop("fit_lcsc(): smooth ", s$label, " is not a cubic regression spline (bs = 'cr')")
    if (n.pen[i] != 1) stop("fit_lcsc(): smooth ", s$label, " must carry exactly one penalty")
    term <- s$term
    if (!term %in% names(data)) stop("fit_lcsc(): smooth ", s$label, " is over '", term, "', which is not a column of `data`")
    k <- s$bs.dim
    cols <- s$first.para:s$last.para
    x <- data[[term]]
    if (is.null(grids[[term]])) {
      # no grid supplied
      u <- sort(unique(x))
      if (length(u) > 2 && any(abs(diff(u) - min(diff(u))) > 1e-8))
        warning("fit_lcsc(): the observed values of '", term, "' are not evenly spaced and no grid ",
                "was supplied for it; the basis is built on the observed values only. Pass cfg = ",
                "lcsc_config(...) (or grids = list(", term, " = <full grid>)).", call. = FALSE)
    }
    ob <- orthogonalize_basis(x_obs = x, k = k, grid = grids[[term]])
    if (ncol(ob$cur) != length(cols))
      stop("fit_lcsc(): column mismatch for ", s$label, " (", ncol(ob$cur), " vs ", length(cols), ")")
    # write the orthogonalized basis and penalty into the prefit object
    G$X[, cols] <- ob$cur
    G$S[[pen.index[i]]] <- ob$S
    G$smooth[[i]]$S[[1]] <- ob$S
    if (!legacy_ubre) {
      # correct the penalty-rank bookkeeping
      rk <- ncol(ob$S)
      G$rank[pen.index[i]] <- rk
      G$smooth[[i]]$rank <- rk
      G$smooth[[i]]$null.space.dim <- 0
    }
    names(ortho)[i] <- s$label
    ortho[[i]] <- c(ob[c("values", "Z", "S", "sm", "k")],
                    list(type = "cr", term = term, cols = cols, label = s$label))
  }

  ## ---- (ii) factor-smooth interactions: orthogonalized to the common terms ------------- #
  if (any(is.fs) && !fs_project) {
    for (i in which(is.fs)) {
      s <- G$smooth[[i]]
      names(ortho)[i] <- s$label
      ortho[[i]] <- list(type = "fs", term = s$term, cols = s$first.para:s$last.para, label = s$label, sm = s,
                         B = NULL, common.cols = integer(0), common.names = character(0), weights = NULL)
    }
  }
  if (any(is.fs) && fs_project) {
    cn <- G$term.names
    fs.cols <- unlist(lapply(G$smooth[is.fs], function(s) s$first.para:s$last.para))
    cr.cols <- unlist(lapply(G$smooth[!is.fs], function(s) s$first.para:s$last.para))
    para.cols <- setdiff(seq_len(ncol(G$X)), c(fs.cols, cr.cols))
    excl <- if (is.null(fs_exclude)) integer() else
      para.cols[cn[para.cols] %in% fs_exclude | grepl(paste(fs_exclude, collapse = "|"), cn[para.cols])]
    common.cols <- sort(c(setdiff(para.cols, excl), cr.cols))
    for (i in which(is.fs)) {
      s <- G$smooth[[i]]
      cols <- s$first.para:s$last.para
      ob <- orthogonalize_block(G$X[, cols, drop = FALSE], G$X[, common.cols, drop = FALSE], w = fs_weights)
      G$X[, cols] <- ob$Hperp
      names(ortho)[i] <- s$label
      ortho[[i]] <- list(type = "fs", term = s$term, cols = cols, label = s$label, sm = s,
                         B = ob$B, common.cols = common.cols, common.names = cn[common.cols],
                         weights = fs_weights)
    }
  }

  ## ---- (iii) the fit ------------------------------------------------------------------ #
  if (legacy_ubre) {
    mod <- mgcv::gam(G = G)                       # original behaviour: default GCV.Cp/UBRE
  } else {
    ctrl <- if (is.null(nthreads)) mgcv::gam.control() else mgcv::gam.control(nthreads = nthreads)
    mod <- mgcv::gam(G = G, method = method, control = ctrl, in.out = in.out)
  }
  list(mod = mod, X = G$X, S = G$S, ortho = ortho, data = data, cfg = cfg, grids = grids,
       method = if (legacy_ubre) mod$method else method, legacy_ubre = legacy_ubre)
}

# fit_or_load(): fit_lcsc() with the result cached as <dir>/<name>.<criterion>.rds
lcsc_first_run_notice <- function(prefix, total_min, dir = "Models") {
  cached <- list.files(dir, pattern = paste0("^", prefix, ".*\\.rds$"))
  if (length(cached)) {
    cat(sprintf("Cached fits found in %s/ (%d files); models with unchanged settings are reloaded, not re-estimated.\n\n",
                dir, length(cached)))
  } else {
    bar <- strrep("-", 78)
    cat(bar, "\n", sep = "")
    cat(sprintf("FIRST RUN: no cached fits in %s/. The models below are estimated from scratch,\n", dir))
    cat(sprintf("which takes about %s minutes in total on a workstation (each fit prints its own\n", format(total_min)))
    cat("expected and actual time). The fitted models are then cached in ", dir, "/, so that\n", sep = "")
    cat("later runs of this script take about a minute.\n")
    cat(bar, "\n\n", sep = "")
  }
  invisible(length(cached))
}
# Hash serialized values without adding a package dependency
lcsc_fingerprint <- function(x) {
  path <- tempfile("lcsc-fingerprint-", fileext = ".rds")
  on.exit(unlink(path), add = TRUE)
  saveRDS(x, path, compress = FALSE, version = 3)
  unname(tools::md5sum(path))
}

lcsc_cache_inputs <- function(G, data, cfg, grids) {
  fields <- c("X", "y", "w", "offset", "S", "off", "rank", "C", "L", "lsp0",
              "sp", "min.sp", "H", "nsdf", "n.true", "scale", "smooth")
  family_code <- lapply(G$family, function(x) {
    if (is.function(x)) list(formals = formals(x), body = body(x)) else x
  })
  fitting_code <- lapply(list(fit_lcsc, orthogonalize_basis, orthogonalize_block,
                              basis_at, lcsc_predict), function(f) {
    list(formals = formals(f), body = body(f))
  })
  lcsc_fingerprint(list(data = data, cfg = cfg, grids = grids,
                        prefit = G[intersect(fields, names(G))],
                        formula = paste(deparse(G$formula, width.cutoff = 500), collapse = " "),
                        family = family_code, implementation = fitting_code))
}

fit_or_load <- function(name, G, data, cfg = NULL, method = c("REML", "ML"), refit = FALSE,
                        dir = "Models", grids = NULL, fs_exclude = NULL, fs_weights = NULL,
                        in.out = NULL, legacy_ubre = FALSE, nthreads = NULL, formula = NULL,
                        expected_min = NULL, quiet = FALSE, fs_project = TRUE) {
  method <- match.arg(method)
  dir.create(dir, showWarnings = FALSE, recursive = TRUE)
  path <- file.path(dir, paste0(name, ".", if (legacy_ubre) "UBRE" else method, ".rds"))
  settings <- list(cache_schema = 2L,
                   R = as.character(getRversion()), platform = R.version$platform,
                   mgcv = as.character(utils::packageVersion("mgcv")),
                   inputs = lcsc_cache_inputs(G, data, cfg, grids),
                   method = method, legacy_ubre = legacy_ubre, fs_exclude = fs_exclude,
                   fs_weights = if (is.null(fs_weights)) "equal" else fs_weights,
                   in.out = in.out, nthreads = nthreads, fs_project = fs_project,
                   formula = if (is.null(formula)) "" else paste(deparse(formula, width.cutoff = 500), collapse = " "),
                   n = nrow(data), p = ncol(G$X))
  if (!refit && file.exists(path)) {
    rs <- tryCatch(readRDS(path), error = function(e) NULL)
    if (!is.null(rs) && identical(rs$settings, settings) && isTRUE(rs$mod$converged)) {
      if (!quiet) cat(sprintf("[%s, %s] loaded cached fit from %s (mgcv method = %s; fitted in %.1f min)\n",
                              name, method, path, rs$mod$method, rs$elapsed.min))
      return(rs)
    }
    if (!quiet) cat(sprintf("[%s, %s] cache in %s is absent, invalid, or has different inputs/settings -- re-estimating\n",
                            name, method, path))
  }
  if (!quiet) cat(sprintf("[%s, %s] no usable cached fit in %s/ -- fitting with mgcv::gam(method = '%s'%s)%s ...\n",
                          name, method, dir, if (legacy_ubre) "<default: UBRE>" else method,
                          if (is.null(nthreads)) "" else sprintf(", nthreads = %d", nthreads),
                          if (is.null(expected_min)) "" else sprintf("; expect about %s min", format(expected_min))))
  t0 <- Sys.time()
  rs <- fit_lcsc(G, data = data, cfg = cfg, method = method, grids = grids, fs_exclude = fs_exclude,
                 fs_weights = fs_weights, legacy_ubre = legacy_ubre, nthreads = nthreads, in.out = in.out,
                 fs_project = fs_project)
  el <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
  if (!quiet) cat(sprintf("[%s, %s] fitted in %.1f min (mgcv method = %s; converged = %s; outer iterations = %s; %d coefficients)\n",
                          name, method, el, rs$mod$method, rs$mod$converged,
                          paste(rs$mod$outer.info$iter, collapse = ""), length(coef(rs$mod))))
  if (!isTRUE(rs$mod$converged)) stop(sprintf("[%s, %s] mgcv reports non-convergence; fit was not cached", name, method))
  rs$elapsed.min <- el; rs$settings <- settings
  saveRDS(rs, path)
  if (!quiet) cat(sprintf("[%s, %s] cached to %s (%.0f MB)\n", name, method, path, file.size(path) / 1e6))
  rs
}

# 4. PREDICTION

# lcsc_predict(): prediction for an LC-SC model fitted with fit_lcsc()
lcsc_predict <- function(fit, newdata) {
  mod <- fit$mod
  beta <- mod$coefficients
  X <- matrix(0, nrow = nrow(newdata), ncol = length(beta))
  colnames(X) <- names(beta)
  # parametric columns (intercept + parametric terms)
  Terms <- delete.response(mod$pterms)
  term.labels <- attr(mod$pterms, "term.labels")
  need <- setdiff(c(all.vars(Terms), unlist(lapply(fit$ortho, `[[`, "term"))), names(newdata))
  if (length(need)) stop("lcsc_predict(): newdata lacks the column(s) ", paste(need, collapse = ", "))
  Xp <- model.matrix(Terms, newdata)
  X[, seq_len(ncol(Xp))] <- Xp
  n.pterms <- length(term.labels)
  # smooth columns: cr smooths first (an fs block is projected off them)
  types <- sapply(fit$ortho, function(ob) if (is.null(ob$type)) "cr" else ob$type)
  for (i in which(types == "cr")) {
    ob <- fit$ortho[[i]]
    X[, ob$cols] <- basis_at(ob, newdata[[ob$term]])
  }
  for (i in which(types == "fs")) {
    ob <- fit$ortho[[i]]
    H <- mgcv::PredictMat(ob$sm, newdata)
    X[, ob$cols] <- if (is.null(ob$B)) H else H - X[, ob$common.cols, drop = FALSE] %*% ob$B
  }
  fam <- mod$family
  linkinv <- fam$linkinv; dmu.deta <- fam$mu.eta
  Vp <- mod$Vp
  lp <- as.vector(X %*% beta)
  se.lp <- sqrt(pmax(0, rowSums((X %*% Vp) * X)))
  Yhat <- data.frame(est = linkinv(lp), se = se.lp * abs(dmu.deta(lp)))
  # term-wise contributions
  termNames <- c(term.labels, unname(sapply(fit$ortho, `[[`, "label")))
  est <- matrix(0, nrow = nrow(newdata), ncol = length(termNames))
  colnames(est) <- termNames
  se <- est
  assign <- attr(Xp, "assign")
  for (j in seq_len(n.pterms)) {
    jc <- which(assign == j)
    Xj <- X[, jc, drop = FALSE]
    e <- as.vector(Xj %*% beta[jc])
    s <- sqrt(pmax(0, rowSums((Xj %*% Vp[jc, jc, drop = FALSE]) * Xj)))
    est[, j] <- e; se[, j] <- s * abs(dmu_eta_safe(dmu.deta, e))
  }
  for (i in seq_along(fit$ortho)) {
    jc <- fit$ortho[[i]]$cols
    Xj <- X[, jc, drop = FALSE]
    e <- as.vector(Xj %*% beta[jc])
    s <- sqrt(pmax(0, rowSums((Xj %*% Vp[jc, jc, drop = FALSE]) * Xj)))
    est[, n.pterms + i] <- e; se[, n.pterms + i] <- s * abs(dmu_eta_safe(dmu.deta, e))
  }
  list(est = est, se = se, X = X, lp = lp, Yhat = Yhat, intercept = unname(beta[1]))
}
dmu_eta_safe <- function(f, eta) { r <- f(eta); r[!is.finite(r)] <- 0; r }

# lcsc_component(): linear predictor and standard error of a set of design columns
lcsc_component <- function(fit, X, cols, level = 0.95) {
  beta <- fit$mod$coefficients; Vp <- fit$mod$Vp
  Xj <- X[, cols, drop = FALSE]
  lp <- as.vector(Xj %*% beta[cols])
  se <- sqrt(pmax(0, rowSums((Xj %*% Vp[cols, cols, drop = FALSE]) * Xj)))
  z <- qnorm(1 - (1 - level) / 2)
  data.frame(lp = lp, se = se, lo = lp - z * se, hi = lp + z * se)
}

# lcsc_coef_index(): the positions of each block of coefficients, found by name
lcsc_coef_index <- function(mod) {
  if (!inherits(mod, "gam") && !is.null(mod[["mod"]])) mod <- mod[["mod"]]   # accept a fit_lcsc() result as well
  nm <- names(stats::coef(mod))
  out <- list()
  out[["(Intercept)"]] <- which(nm == "(Intercept)")
  for (t in attr(mod$pterms, "term.labels")) {
    hit <- which(nm == t)
    if (length(hit)) out[[t]] <- hit
  }
  for (k in seq_along(mod$smooth)) {
    sm <- mod$smooth[[k]]
    out[[sm$label]] <- sm$first.para:sm$last.para
  }
  out
}

# 5. PERIOD-COHORT SHOCK INDICATORS

# create_dummy_vars(): one indicator per (col1, col2) cell
create_dummy_vars <- function(df, col1_name, col2_name, col1_range, col2_range) {
  dummy_var_names <- c()
  for (val1 in col1_range) {
    for (val2 in col2_range) {
      var_name <- paste0(col1_name, val1, "_", col2_name, val2)
      df[[var_name]] <- ifelse(df[[col1_name]] == val1 & df[[col2_name]] == val2, 1, 0)
      if (sum(df[[var_name]]) > 0) {
        dummy_var_names <- c(dummy_var_names, var_name)
      } else {
        df[[var_name]] <- NULL   # cell not observed: drop the all-zero column
      }
    }
  }
  list(df = df, dummy_var_names = dummy_var_names)
}

# 6. LEXIS-TABLE MATRICES FOR THE FIGURES

mean_by_ap <- function(data_frame) {
  tapply(data_frame[, 4], list(data_frame$a.index, data_frame$p.index), mean)
}
mean_by_pc <- function(data_frame) {
  tapply(data_frame[, 4], list(data_frame$p.index, data_frame$c.index), mean)
}
mean_by_ac <- function(data_frame) {
  tapply(data_frame[, 4], list(data_frame$a.index, data_frame$c.index), mean)
}

# lcsc_axis_offsets(): heat-map label positions, derived from the data
lcsc_axis_offsets <- function(index_values, by = 5, origin_adj = 1) {
  labs <- seq(from = min(index_values), to = max(index_values), by = by)
  list(labs = labs, loc = labs - (min(index_values) - origin_adj))
}

# 7. CONCURVITY ON THE ORTHOGONALIZED DESIGN

# concurvity_lcsc(): mgcv::concurvity() computed on a supplied design matrix
concurvity_lcsc <- function(b, X, full = TRUE) {
  m <- length(b$smooth)
  X <- X[rowSums(is.na(X)) == 0, , drop = FALSE]
  X <- qr.R(qr(X, tol = 0, LAPACK = FALSE))
  stop <- start <- rep(1, m); lab <- rep("", m)
  for (i in 1:m) {
    start[i] <- b$smooth[[i]]$first.para; stop[i] <- b$smooth[[i]]$last.para
    lab[i] <- b$smooth[[i]]$label
  }
  if (min(start) > 1) {
    start <- c(1, start); stop <- c(min(start) - 1, stop); lab <- c("para", lab); m <- m + 1
  }
  measure.names <- c("worst", "observed", "estimate")
  if (full) {
    conc <- matrix(0, 3, m)
    for (i in 1:m) {
      Xi <- X[, -(start[i]:stop[i]), drop = FALSE]
      Xj <- X[, start[i]:stop[i], drop = FALSE]
      r <- ncol(Xi)
      R <- qr.R(qr(cbind(Xi, Xj), LAPACK = FALSE, tol = 0))[, -(1:r), drop = FALSE]
      Rt <- qr.R(qr(R, tol = 0))
      conc[1, i] <- svd(forwardsolve(t(Rt), t(R[1:r, , drop = FALSE])))$d[1]^2
      beta <- b$coef[start[i]:stop[i]]
      conc[2, i] <- sum((R[1:r, , drop = FALSE] %*% beta)^2) / sum((Rt %*% beta)^2)
      conc[3, i] <- sum(R[1:r, ]^2) / sum(R^2)
    }
    colnames(conc) <- lab; rownames(conc) <- measure.names
  } else {
    conc <- list()
    for (i in 1:3) conc[[i]] <- matrix(1, m, m)
    for (i in 1:m) {
      Xi <- X[, start[i]:stop[i], drop = FALSE]; r <- ncol(Xi)
      for (j in 1:m) if (i != j) {
        Xj <- X[, start[j]:stop[j], drop = FALSE]
        R <- qr.R(qr(cbind(Xi, Xj), LAPACK = FALSE, tol = 0))[, -(1:r), drop = FALSE]
        Rt <- qr.R(qr(R, tol = 0))
        conc[[1]][i, j] <- svd(forwardsolve(t(Rt), t(R[1:r, , drop = FALSE])))$d[1]^2
        beta <- b$coef[start[j]:stop[j]]
        conc[[2]][i, j] <- sum((R[1:r, , drop = FALSE] %*% beta)^2) / sum((Rt %*% beta)^2)
        conc[[3]][i, j] <- sum(R[1:r, ]^2) / sum(R^2)
      }
    }
    for (i in 1:3) rownames(conc[[i]]) <- colnames(conc[[i]]) <- lab
    names(conc) <- measure.names
  }
  conc
}

# concurvity_blocks(): the same three measures on an explicit partition of the columns
concurvity_blocks <- function(X, beta, blocks, targets = names(blocks), tol = 1e-8,
                              w = NULL, full = TRUE) {
  if (is.null(names(blocks)) || any(names(blocks) == "")) stop("concurvity_blocks(): `blocks` must be named")
  if (!identical(sort(unlist(blocks, use.names = FALSE)), seq_len(ncol(X))))
    stop("concurvity_blocks(): `blocks` must partition all ", ncol(X), " columns of X ",
         "(every column belongs to exactly one block)")
  if (!all(targets %in% names(blocks)))
    stop("concurvity_blocks(): unknown target(s): ", paste(setdiff(targets, names(blocks)), collapse = ", "))
  if (length(beta) != ncol(X))
    stop("concurvity_blocks(): beta has ", length(beta), " entries but X has ", ncol(X), " columns")
  if (!is.null(w)) {
    if (length(w) != nrow(X) || any(w < 0)) stop("concurvity_blocks(): `w` must be nrow(X) non-negative weights")
    X <- X * sqrt(w)
  }

  ## ---- (i) rank-truncate each block, preserving its fitted contribution --------------- #
  U <- vector("list", length(blocks)); g <- vector("list", length(blocks))
  M <- vector("list", length(blocks))
  nc <- rk <- integer(length(blocks))
  for (i in seq_along(blocks)) {
    cj <- blocks[[i]]; Xj <- X[, cj, drop = FALSE]; bj <- beta[cj]
    sv <- svd(Xj); r <- sum(sv$d > max(sv$d) * tol)
    nc[i] <- ncol(Xj); rk[i] <- r
    if (r == ncol(Xj)) {
      U[[i]] <- Xj; g[[i]] <- bj                      # untouched: M stays NULL, as initialized
    } else {
      U[[i]] <- sv$u[, 1:r, drop = FALSE]             # orthonormal basis of the same span
      M[[i]] <- sv$d[1:r] * t(sv$v[, 1:r, drop = FALSE])       # X_j = U M exactly
      g[[i]] <- M[[i]] %*% bj                         # so that U g = X_j b_j
    }
  }
  Xs <- do.call(cbind, U); bs <- unlist(lapply(g, as.numeric))
  start <- cumsum(c(1, head(rk, -1))); stop <- cumsum(rk)

  ## ---- (ii) mgcv's measures on the well-posed design ---------------------------------- #
  Xq <- qr.R(qr(Xs, tol = 0, LAPACK = FALSE))
  measure <- function(ci, cj, Mi) {     # ci: target columns; cj: comparison columns
    Xi <- Xq[, cj, drop = FALSE]; Xj <- Xq[, ci, drop = FALSE]
    # the comparison columns can be jointly rank deficient
    if (qr(Xi, tol = 1e-10)$rank < ncol(Xi)) {
      sv <- svd(Xi, nv = 0); Xi <- sv$u[, sv$d > max(sv$d) * tol, drop = FALSE]
    }
    r <- ncol(Xi)
    R  <- qr.R(qr(cbind(Xi, Xj), LAPACK = FALSE, tol = 0))[, -(1:r), drop = FALSE]
    Rt <- qr.R(qr(R, tol = 0)); bb <- bs[ci]
    # "estimate" averages over the block's own basis functions
    RM <- if (is.null(Mi)) R[1:r, , drop = FALSE] else R[1:r, , drop = FALSE] %*% Mi
    RtM <- if (is.null(Mi)) Rt else Rt %*% Mi
    c(worst    = svd(forwardsolve(t(Rt), t(R[1:r, , drop = FALSE])))$d[1]^2,
      observed = sum((R[1:r, , drop = FALSE] %*% bb)^2) / sum((Rt %*% bb)^2),
      estimate = sum(RM^2) / sum(RtM^2))
  }
  ti <- match(targets, names(blocks))
  if (full) {
    conc <- matrix(NA_real_, 3, length(ti),
                   dimnames = list(c("worst", "observed", "estimate"), targets))
    for (k in seq_along(ti)) {
      ci <- start[ti[k]]:stop[ti[k]]
      conc[, k] <- measure(ci, setdiff(seq_len(ncol(Xq)), ci), M[[ti[k]]])
    }
  } else {
    conc <- lapply(1:3, function(z) matrix(1, length(ti), length(blocks),
                   dimnames = list(targets, names(blocks))))
    for (k in seq_along(ti)) for (j in seq_along(blocks)) if (ti[k] != j) {
      v <- measure(start[ti[k]]:stop[ti[k]], start[j]:stop[j], M[[ti[k]]])
      for (z in 1:3) conc[[z]][k, j] <- v[z]
    }
    names(conc) <- c("worst", "observed", "estimate")
  }
  attr(conc, "ncol") <- setNames(nc, names(blocks))
  attr(conc, "rank") <- setNames(rk, names(blocks))
  attr(conc, "deficient") <- names(blocks)[rk < nc]
  conc
}

# lcsc_blocks(): the column partition of a fit_lcsc() model, with the paper's term names
lcsc_blocks <- function(b, X, extra = list()) {
  first <- sapply(b$smooth, `[[`, "first.para"); last <- sapply(b$smooth, `[[`, "last.para")
  nm <- if (!is.null(colnames(X))) colnames(X) else names(coef(b))
  para <- seq_len(min(first) - 1L)
  blocks <- list(); blocks[["intercept"]] <- para[1]; used <- para[1]
  for (e in names(extra)) {
    cj <- if (is.character(extra[[e]])) which(grepl(extra[[e]], nm)) else extra[[e]]
    cj <- setdiff(intersect(cj, para), used)
    if (length(cj)) { blocks[[e]] <- cj; used <- c(used, cj) }
  }
  rest <- setdiff(para, used); if (length(rest)) blocks[["linear"]] <- rest
  for (i in seq_along(b$smooth)) blocks[[b$smooth[[i]]$label]] <- first[i]:last[i]
  tail.cols <- setdiff(seq_len(ncol(X)), unlist(blocks, use.names = FALSE))
  if (length(tail.cols)) blocks[["other"]] <- tail.cols
  blocks
}

# 8. FIT STATISTICS (the paper's definitions)

# deviance-based R2
deviance_r2 <- function(residual_deviance, null_deviance) 1 - residual_deviance / null_deviance

# deviance-based adjusted R2
adj_deviance_r2 <- function(residual_deviance, null_deviance, n, edf) {
  r2 <- 1 - residual_deviance / null_deviance
  1 - (1 - r2) * ((n - 1) / (n - edf))
}

# Richards
richards_chat <- function(LL_model, LL_sat, edf_model, edf_sat) {
  2 * (LL_sat - LL_model) / (edf_sat - edf_model)
}

# QAIC and QBIC for a Poisson model given c_hat, with one extra degree of freedom
qaic_qbic <- function(LL, edf, c_hat, n) {
  dof <- edf + 1
  c(QAIC = -2 * LL / c_hat + 2 * dof, QBIC = -2 * LL / c_hat + dof * log(n))
}

# analytic log-likelihood of the saturated Poisson model
poisson_saturated_loglik <- function(y) sum(dpois(y, lambda = y, log = TRUE))

# Poisson deviance computed safely
safe_deviance <- function(y, mu) sum(pmax(poisson()$dev.resids(y, mu, rep(1, length(y))), 0))

# Pearson dispersion of a fitted Poisson GAM: Pearson chi-square / (n - EDF)
pearson_dispersion <- function(model, mu = NULL, edf = NULL) {
  if (is.null(mu)) { y <- model$y; mu <- fitted(model); edf <- sum(model$edf) } else y <- model
  sum((y - mu)^2 / mu) / (length(y) - edf)
}

# EDF of the i-th smooth of a fitted mgcv model
smooth_edf <- function(m, i) { s <- m$smooth[[i]]; sum(m$edf[s$first.para:s$last.para]) }

# period_from_residuals(): average of the cell-specific residuals for each period
period_from_residuals <- function(y, mu, exposure, period) {
  per <- sort(unique(period))
  O <- tapply(y, period, sum)[as.character(per)]
  E <- tapply(mu, period, sum)[as.character(per)]
  rr <- tapply((y - mu) / exposure * 1000, period, mean)[as.character(per)]
  data.frame(period = per, observed = as.numeric(O), expected = as.numeric(E),
             ratio = as.numeric(O / E), log_ratio = as.numeric(log(O / E)),
             mean_resid_per1000 = as.numeric(rr), row.names = NULL)
}

# period_leftovers(): the common period fluctuations of the models with varying LC curves
period_leftovers <- function(y, mu, period, exclude = rep(FALSE, length(y))) {
  stopifnot(length(y) == length(mu), length(period) == length(y), length(exclude) == length(y))
  per <- sort(unique(period))
  keep <- !exclude
  O <- as.numeric(tapply(y[keep], factor(period[keep], levels = per), sum))
  E <- as.numeric(tapply(mu[keep], factor(period[keep], levels = per), sum))
  if (any(!is.finite(O / E)) || any(O <= 0)) stop("period_leftovers(): a period has no usable cells or no events")
  g <- log(O / E)
  g_cell <- g[match(period, per)]
  list(table = data.frame(period = per, observed = O, expected = E, g = g, se = 1 / sqrt(O)),
       g_cell = g_cell, mu = ifelse(exclude, mu, mu * exp(g_cell)))
}

# Format a number as it appears in the manuscript
fmt_tex <- function(x, digits = 2, big.mark = ",") {
  out <- formatC(round(x, digits), format = "f", digits = digits, big.mark = big.mark)
  gsub(" ", "", out)
}

# significance stars as used in the manuscript's model-summary table
sig_stars <- function(p) ifelse(p < 0.001, "^{***}", ifelse(p < 0.01, "^{**}", ifelse(p < 0.05, "^{*}", "")))

# write LaTeX lines to a file (UTF-8, LF line endings), creating the directory if needed
write_tex <- function(lines, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  con <- file(path, open = "wb"); on.exit(close(con))
  writeLines(lines, con, sep = "\n", useBytes = TRUE)
  invisible(path)
}

# registry of the numbers quoted in the manuscript text (one CSV per example)
.mn <- new.env()
.mn$tab <- data.frame(key = character(), value = character(), numeric = numeric(),
                      where = character(), stringsAsFactors = FALSE)

# mn_add(): register a number quoted in the manuscript text
mn_add <- function(key, value, digits = 2, where = "") {
  v <- if (is.numeric(value)) fmt_tex(value, digits) else as.character(value)
  num <- if (is.numeric(value)) as.numeric(value) else suppressWarnings(as.numeric(gsub(",", "", value)))
  .mn$tab <- rbind(.mn$tab, data.frame(key = key, value = v, numeric = num, where = where,
                                      stringsAsFactors = FALSE))
  invisible(v)
}

# 10. GUARD RAILS FOR THE ANALYSIS SCRIPTS

# lcsc_require_data(): the first thing every analysis script calls
lcsc_require_data <- function(file, what = basename(file)) {
  if (!file.exists(file)) {
    stop('Cannot find ', file, '.\n',
         'Set the working directory to the folder containing this script and Data/.\n',
         'This distribution does not include the dataset; see README.md for its source and required columns.',
         call. = FALSE)
  }
  invisible(TRUE)
}

# lcsc_session_header(): a few lines at the top of every run saying
lcsc_session_header <- function(script) {
  cat("\n", strrep("=", 78), "\n", sep = "")
  cat(script, "\n")
  cat("  working directory : ", getwd(), "\n", sep = "")
  cat("  R version         : ", R.version.string, "\n", sep = "")
  cat("  mgcv version      : ", as.character(utils::packageVersion("mgcv")), "\n", sep = "")
  cat("  started           : ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n", sep = "")
  cat(strrep("=", 78), "\n", sep = "")
  invisible(NULL)
}

# timing wrapper
timed <- function(expr, label = "") {
  t0 <- Sys.time()
  r <- expr
  cat(sprintf("[%s] elapsed: %.2f min\n", label, as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  r
}

# print a named numeric vector / matrix rounded
pr <- function(x, digits = 4) print(round(x, digits))

message("Embedded LC-SC helpers loaded: lcsc_config(), fit_lcsc(), lcsc_predict() and helpers.")
# ---- END EMBEDDED HELPERS; EXAMPLE ANALYSIS CONTINUES ---------------------------

lcsc_session_header("02_mortality_france.R")
lcsc_require_data("Data/mortality.RData", "HMD France files")

# loading libraries
suppressPackageStartupMessages({
  library(mgcv)         # penalized-spline GAMs (Simon Wood)
  library(plot3D)       # 3D plotting of surfaces (hist3D)
  library(RColorBrewer) # color palettes
  library(ggplot2)      # ggplot graphics
  library(ggridges)     # ridgeline plot
  library(dplyr)        # %>% and mutate (used in the ridgeline plot)
  library(tidyr)        # gather (used in the ridgeline plot)
})

cat("Settings: CRITERION_FIGURES =", CRITERION_FIGURES, "| CRITERION_TABLES =", CRITERION_TABLES,
    "| LEGACY_UBRE =", LEGACY_UBRE, "| NTHREADS =", NTHREADS, "| REFIT =", REFIT,
    "| ROBUSTNESS_1941_42 =", ROBUSTNESS_1941_42, "\n\n")

# Figures below place some labels with jitter()
set.seed(20240617)

# output directories
fig.dir <- file.path("Figures", "Mortality")
dir.create(fig.dir, recursive = TRUE, showWarnings = FALSE)
dir.create("Output", showWarnings = FALSE)
dir.create("Models", showWarnings = FALSE)

#### 2. Data Preparation ####

# loading the dataset (data frame `df`: HMD France 1x1 male deaths and exposures)
load("Data/mortality.RData")

# age 0 is unusual, and extreme old age is sparse in early years, so we keep ages 1-90
df <- subset(df, age >= 1)
df <- subset(df, age <= 90)

# listwise deletion of missing values
df <- na.omit(df)

# creating a, p, and c variables on their natural scales (cohort = period - age)
df$p <- df$period
df$a <- df$age
df$c <- df$p - df$a

# outcome and exposure (offset)
df$y <- as.integer(round(df$Male))
df$pop <- df$MaleExposure
# thus: we are modeling the age-specific mortality rate in an age-period cell

# unique cell identifier
df$id <- 1:nrow(df)

# war-shock dummies (create_dummy_vars in the embedded helpers)

# Crimean War (1853-1856)
rs <- create_dummy_vars(df, "a", "p", WAR_AGES, 1853:1856)
df <- rs$df; crimean <- rs$dummy_var_names # 144 dummies

# Franco-Prussian War (1870-1871)
rs <- create_dummy_vars(df, "a", "p", WAR_AGES, 1870:1871)
df <- rs$df; prussia <- rs$dummy_var_names # 72 dummies

# World War 1 (1914-1918)
rs <- create_dummy_vars(df, "a", "p", WAR_AGES, 1914:1918)
df <- rs$df; ww1 <- rs$dummy_var_names # 180 dummies

# World War 2, 1940 (fall of France)
rs <- create_dummy_vars(df, "a", "p", WAR_AGES, 1940)
df <- rs$df; ww2.1940 <- rs$dummy_var_names # 36 dummies

# World War 2, 1943-1944 (occupation and liberation)
rs <- create_dummy_vars(df, "a", "p", WAR_AGES, 1943:1944)
df <- rs$df; ww2.1943 <- rs$dummy_var_names # 72 dummies

# all 504 war-shock dummies
war.dummies <- c(crimean, prussia, ww1, ww2.1940, ww2.1943)
stopifnot(length(war.dummies) == 504)

# a.index, p.index, and c.index variables in the data (original scales)
df$a.index <- df$a; df$p.index <- df$p; df$c.index <- df$c

# index vectors of the unique age, period, and cohort values (for plot labels)
a.index <- as.numeric(as.character(sort(unique(df$a)))) # ages 1-90
p.index <- as.numeric(as.character(sort(unique(df$p)))) # periods 1816-2020
c.index <- as.numeric(as.character(sort(unique(df$c)))) # cohorts 1726-2019

# observed age-specific mortality rate (ASMR) per 1,000 males
df$yhat.raw <- (df$Male / df$MaleExposure) * 1000

# ---- CONFIGURATION: the only place the analysis constants are written --------------- #
cfg <- lcsc_config(df,
                   age_center    = 46,
                   period_center = 1918,
                   cohort_center = 1872,
                   period_full   = 1816:2020,
                   knots         = c(a = 90, p = 205, c = 20))

# centered modeling columns a, p, c (overwriting the natural-scale ones)
df <- lcsc_center(df, cfg)

cat("Analysis sample:", nrow(df), "age-period cells;", length(a.index), "ages,",
    length(p.index), "periods,", length(c.index), "cohorts;", length(war.dummies), "war-shock dummies\n\n")

#### 3. Functions for the Analysis ####

# convert the output of concurvity_lcsc() (full and pairwise) into a long data frame
concurvity_to_df <- function(conc.full, conc.pair, model_name) {
  d1 <- data.frame(model = model_name, type = "full",
                   measure = rep(rownames(conc.full), times = ncol(conc.full)),
                   term = rep(colnames(conc.full), each = nrow(conc.full)),
                   term_given = NA_character_,
                   value = as.vector(conc.full), stringsAsFactors = FALSE)
  d2 <- do.call(rbind, lapply(names(conc.pair), function(m) {
    M <- conc.pair[[m]]
    data.frame(model = model_name, type = "pairwise", measure = m,
               term = rep(colnames(M), each = nrow(M)),       # column j: term whose dependence is measured
               term_given = rep(rownames(M), times = ncol(M)), # row i: term it is compared with
               value = as.vector(M), stringsAsFactors = FALSE)
  }))
  rbind(d1, d2)
}

#### 4. Observed Mortality-Rate Matrices ####

# observed (raw) ASMRs arranged as age-period, period-cohort, and age-cohort matrices
yhat.raw <- df[, c("a.index", "p.index", "c.index", "yhat.raw", "pop")]
matAP.raw <- mean_by_ap(yhat.raw); matPC.raw <- mean_by_pc(yhat.raw); matAC.raw <- mean_by_ac(yhat.raw)

#### 5. Estimating the Models (Poisson GAMs, REML and ML) ####

## PRELIMINARIES: SETTING UP THE DATA AND FORMULAS

# the estimation dataset
data <- df[, c("a", "p", "c", "y", "a.index", "p.index", "c.index", "pop", war.dummies)]

# number of knots for age, period, and cohort (from the configuration above)
ak <- cfg$knots[["a"]]  # 90 knots for age
pk <- cfg$knots[["p"]]  # 205 knots for period
ck <- cfg$knots[["c"]]  # 20 knots for cohort

# (a) LC-SC model formula
formula.lcsc <- as.formula(y ~ a + c + s(a, bs = "cr", k = ak - 1) + s(p, bs = "cr", k = pk - 1) + s(c, bs = "cr", k = ck - 1))

# (b) LC-SC model with the 504 period-cohort war-shock dummies on the right-hand side
rhs1 <- c("a + c + s(a, bs = 'cr', k = ak-1) + s(p, bs = 'cr', k = pk-1) + s(c, bs = 'cr', k = ck-1)")
rhs2 <- paste(war.dummies, collapse = " + ")
formula.dum <- as.formula(paste("y ~", rhs1, "+", rhs2))

# unfitted Poisson prefits (fit = FALSE)
G.lcsc <- mgcv::gam(formula.lcsc, offset = log(pop), family = 'poisson', data = data, fit = FALSE)
G.dum  <- mgcv::gam(formula.dum,  offset = log(pop), family = 'poisson', data = data, fit = FALSE)
cat("Coefficients: LC-SC model =", ncol(G.lcsc$X), "; LC-SC model with war dummies =", ncol(G.dum$X), "\n\n")

# each fit is cached in Models/
lcsc_first_run_notice("mort_", total_min = if (ROBUSTNESS_1941_42) 15 else 10)
fit_cached <- function(name, G, method, formula, data. = data, expected_min = NULL) {
  fit_or_load(name, G, data = data., cfg = cfg, method = method, refit = REFIT,
              legacy_ubre = LEGACY_UBRE, nthreads = NTHREADS, formula = formula,
              expected_min = expected_min)
}

# (1) LC-SC MODEL WITH WAR DUMMIES, REML
rs.dum <- fit_cached("mort_m.lcsc.dum", G.dum, CRITERION_FIGURES, formula.dum, expected_min = 4)
m.lcsc.dum <- rs.dum$mod

# (2) LC-SC MODEL, ML (m.lcsc.ML) -- fit-statistics table
rs.lcsc.ML <- fit_cached("mort_m.lcsc", G.lcsc, CRITERION_TABLES, formula.lcsc, expected_min = 1)
m.lcsc.ML <- rs.lcsc.ML$mod

# (3) LC-SC MODEL WITH WAR DUMMIES, ML (m.lcsc.dum.ML) -- fit-statistics table
rs.dum.ML <- fit_cached("mort_m.lcsc.dum", G.dum, CRITERION_TABLES, formula.dum, expected_min = 4)
m.lcsc.dum.ML <- rs.dum.ML$mod

# the prefits are no longer needed
rm(G.lcsc, G.dum); invisible(gc())

## SUMMARY OF THE FITS
cat("\nSmoothing criterion actually used by mgcv (m$method):\n")
cat("  m.lcsc.dum    :", m.lcsc.dum$method, "\n")
cat("  m.lcsc.ML     :", m.lcsc.ML$method, "\n")
cat("  m.lcsc.dum.ML :", m.lcsc.dum.ML$method, "\n")

# per-smooth EDFs of the three models (s(a), s(p), s(c)), plus the total EDF
smooth_edfs <- function(m) sapply(seq_along(m$smooth), function(i) smooth_edf(m, i))
edf.tab <- rbind("m.lcsc.dum (REML)"  = c(smooth_edfs(m.lcsc.dum),    total = sum(m.lcsc.dum$edf)),
                 "m.lcsc.ML"          = c(smooth_edfs(m.lcsc.ML),     total = sum(m.lcsc.ML$edf)),
                 "m.lcsc.dum.ML"      = c(smooth_edfs(m.lcsc.dum.ML), total = sum(m.lcsc.dum.ML$edf)))
colnames(edf.tab)[1:3] <- sapply(m.lcsc.dum$smooth, `[[`, "label")
cat("\nEffective degrees of freedom (smooth terms and total):\n"); pr(edf.tab, 2)
cat("\nSmoothing parameters (log10 lambda): m.lcsc.dum:", round(log10(m.lcsc.dum$sp), 3),
    "| m.lcsc.ML:", round(log10(m.lcsc.ML$sp), 3), "| m.lcsc.dum.ML:", round(log10(m.lcsc.dum.ML$sp), 3), "\n")
cat("\nSmooth-term table of the figures model (m.lcsc.dum):\n")
print(summary(m.lcsc.dum)$s.table)

# EDF of each smooth next to its maximum
edf.max <- sapply(m.lcsc.ML$smooth, function(s) s$last.para - s$first.para + 1)
names(edf.max) <- sapply(m.lcsc.ML$smooth, `[[`, "label")
cat("\nMaximum EDF of each smooth (orthogonalized basis dimension):", paste(names(edf.max), edf.max, collapse = ", "), "\n")
for (i in seq_along(edf.max)) {
  lab <- gsub("[()]", "", names(edf.max)[i])
  mn_add(paste0("mort_edf_", lab, "_lcsc"), smooth_edf(m.lcsc.ML, i), 1,
         paste("Section III.3 footnote: EDF of", names(edf.max)[i], "in the LC-SC model (ML)"))
  mn_add(paste0("mort_edf_", lab, "_dum"), smooth_edf(m.lcsc.dum.ML, i), 1,
         paste("Section III.3 footnote (optional): EDF of", names(edf.max)[i], "in the model with varying period shocks (ML)"))
  mn_add(paste0("mort_edfmax_", lab), edf.max[[i]], 0, paste("Section III.3 footnote: maximum EDF of", names(edf.max)[i]))
}

#### 6. Model Fit Statistics (French Mortality Fit-Statistics Table) ####

## (6a) SATURATED POISSON MODEL (ANALYTIC)
n <- nrow(data)
LL.sat <- poisson_saturated_loglik(data$y)
EDF.sat <- n

# (6b) OVERDISPERSION ESTIMATE c_hat
LL.ML <- as.numeric(logLik(m.lcsc.ML));         EDF.ML <- sum(m.lcsc.ML$edf)
LL.dum.ML <- as.numeric(logLik(m.lcsc.dum.ML)); EDF.dum.ML <- sum(m.lcsc.dum.ML$edf)
c_hat <- richards_chat(LL_model = LL.dum.ML, LL_sat = LL.sat, edf_model = EDF.dum.ML, edf_sat = EDF.sat)

cat("\nOverdispersion:\n")
cat("  variance-to-mean ratio of y:", round(var(data$y) / mean(data$y), 2), "\n")
cat("  c_hat (Richards 2008, m.lcsc.dum.ML vs saturated):", round(c_hat, 4), "\n")
cat("  Pearson dispersion, m.lcsc.ML    :", round(pearson_dispersion(m.lcsc.ML), 4), "\n")
cat("  Pearson dispersion, m.lcsc.dum.ML:", round(pearson_dispersion(m.lcsc.dum.ML), 4), "\n")

## (6c) FIT-STATISTICS TABLE
ic.sat <- c(QAIC=NA_real_, QBIC=NA_real_)  # Saturated likelihood is a benchmark
ic.ML     <- qaic_qbic(LL.ML,     EDF.ML,     c_hat, n)
ic.dum.ML <- qaic_qbic(LL.dum.ML, EDF.dum.ML, c_hat, n)

tab <- data.frame(
  Model    = c("Saturated model", "LC-SC model", "LC-SC model with varying period shocks"),
  Criterion = c("--", m.lcsc.ML$method, m.lcsc.dum.ML$method),
  LLV      = c(LL.sat, LL.ML, LL.dum.ML),
  R2_D     = c(1,
               deviance_r2(m.lcsc.ML$deviance, m.lcsc.ML$null.deviance),
               deviance_r2(m.lcsc.dum.ML$deviance, m.lcsc.dum.ML$null.deviance)),
  Adj_R2_D = c(NA, # undefined for the saturated model (EDF = n)
               adj_deviance_r2(m.lcsc.ML$deviance, m.lcsc.ML$null.deviance, n = n, edf = EDF.ML),
               adj_deviance_r2(m.lcsc.dum.ML$deviance, m.lcsc.dum.ML$null.deviance, n = n, edf = EDF.dum.ML)),
  QAIC     = c(ic.sat["QAIC"], ic.ML["QAIC"], ic.dum.ML["QAIC"]),
  QBIC     = c(ic.sat["QBIC"], ic.ML["QBIC"], ic.dum.ML["QBIC"]),
  EDF      = c(EDF.sat, EDF.ML, EDF.dum.ML),
  c_hat    = c_hat,
  Pearson_dispersion = c(NA, pearson_dispersion(m.lcsc.ML), pearson_dispersion(m.lcsc.dum.ML)),
  stringsAsFactors = FALSE)

# console version of the fit-statistics table
tab.print <- tab
tab.print$LLV  <- formatC(tab$LLV, format = "f", digits = 2, big.mark = ",")
tab.print$R2_D <- formatC(tab$R2_D, format = "f", digits = 5)
tab.print$Adj_R2_D <- ifelse(is.na(tab$Adj_R2_D), "---", formatC(tab$Adj_R2_D, format = "f", digits = 5))
tab.print$QAIC <- formatC(tab$QAIC, format = "f", digits = 2, big.mark = ",")
tab.print$QBIC <- formatC(tab$QBIC, format = "f", digits = 2, big.mark = ",")
tab.print$EDF  <- formatC(tab$EDF, format = "f", digits = 2, big.mark = ",")
tab.print$c_hat <- formatC(tab$c_hat, format = "f", digits = 4)
tab.print$Pearson_dispersion <- ifelse(is.na(tab$Pearson_dispersion), "---", formatC(tab$Pearson_dispersion, format = "f", digits = 4))
cat("\nFrench Mortality: Fit Statistics of LC-SC and Saturated Models\n")
cat("(paper table: R2_D 1.000 / 0.966 / 0.990; Adj R2_D --- / 0.966 / 0.989; EDF 18,450.00 / 309.93 / 811.88)\n")
print(tab.print, row.names = FALSE)

# writing the table
write.csv(tab, file.path("Output", "mortality_fit_table.csv"), row.names = FALSE)

# LaTeX rows of the table exactly as typeset in the manuscript
mort_rows <- function(t) {
  r2 <- function(x) if (is.na(x)) "\\multicolumn{1}{c}{\\ \\ \\ ---}" else sprintf("$%s$", fmt_tex(x, 3))
  cell <- function(i) sprintf("$%s$ & $%s$ & %s & $%s$ & $%s$ & $%s$",
                              fmt_tex(t$LLV[i], 2), fmt_tex(t$R2_D[i], 3), r2(t$Adj_R2_D[i]),
                              if(is.na(t$QAIC[i])) "---" else fmt_tex(t$QAIC[i], 2), if(is.na(t$QBIC[i])) "---" else fmt_tex(t$QBIC[i], 2), fmt_tex(t$EDF[i], 2))
  c(paste0("Saturated model & ", cell(1), " \\\\[1.5ex]"),
    paste0("LC-SC model & ", cell(2), " \\\\[1.5ex]"),
    paste0("LC-SC model with \\\\ varying period shocks & ", cell(3), " \\\\[2.5ex]"))
}
write_tex(mort_rows(tab), file.path("Output", "tex", "mortality_fit_table.tex"))
cat("LaTeX rows of the fit-statistics table written to Output/tex/mortality_fit_table.tex\n")

# numbers of this example quoted in the manuscript text
mn_add("mort_n_cells", n, 0, "S x T = 18,450 (everywhere)")
mn_add("mort_n_dummies", length(war.dummies), 0, "footnote: 504 war dummies")
mn_add("mort_war_age_lo", min(WAR_AGES), 0, "footnote: war dummies for ages 15 to 50")
mn_add("mort_war_age_hi", max(WAR_AGES), 0, "footnote: war dummies for ages 15 to 50")
mn_add("mort_c_hat", c_hat, 2, "footnote on QAIC/QBIC: estimated nu (c_hat)")
mn_add("mort_edf_lcsc", EDF.ML, 2, "table: EDF of the LC-SC model (ML)")
mn_add("mort_edf_dum", EDF.dum.ML, 2, "table: EDF of the model with varying period shocks (ML)")
mn_add("mort_pearson_dum", pearson_dispersion(m.lcsc.dum.ML), 0, "footnote: Pearson dispersion of the model with varying period shocks")

## (6c') ROBUSTNESS: ADDITIONAL PERIOD-COHORT DUMMIES FOR 1941-1942 (ML)
if (ROBUSTNESS_1941_42) {
  # NOTE: df$a and df$p are centered at this point
  rs41 <- create_dummy_vars(df, "a.index", "p.index", WAR_AGES, 1941:1942)
  df41 <- rs41$df; vichy <- rs41$dummy_var_names            # 72 dummies
  stopifnot(length(vichy) == 72, all(colSums(df41[, vichy]) == 1))  # each dummy marks exactly one cell
  data41 <- df41[, c("a", "p", "c", "y", "a.index", "p.index", "c.index", "pop", war.dummies, vichy)]
  formula.dum41 <- as.formula(paste("y ~", rhs1, "+", paste(c(war.dummies, vichy), collapse = " + ")))
  G.dum41 <- mgcv::gam(formula.dum41, offset = log(pop), family = 'poisson', data = data41, fit = FALSE)
  rs.dum41.ML <- fit_cached("mort_m.lcsc.dum4142", G.dum41, CRITERION_TABLES, formula.dum41, data. = data41, expected_min = 5)
  m.dum41 <- rs.dum41.ML$mod
  LL.dum41 <- as.numeric(logLik(m.dum41)); EDF.dum41 <- sum(m.dum41$edf)
  c_hat41 <- richards_chat(LL_model = LL.dum41, LL_sat = LL.sat, edf_model = EDF.dum41, edf_sat = EDF.sat)
  ic.dum41.tab <- qaic_qbic(LL.dum41, EDF.dum41, c_hat, n)     # table's nu-hat (504-dummy model)
  ic.dum.own   <- qaic_qbic(LL.dum.ML, EDF.dum.ML, c_hat41, n) # 504-dummy model under the 576-dummy nu-hat
  ic.dum41.own <- qaic_qbic(LL.dum41, EDF.dum41, c_hat41, n)
  rob <- data.frame(
    Model = c("LC-SC model with varying period shocks (504 dummies, ML)",
              "LC-SC model with varying period shocks + 1941-1942 dummies (576 dummies, ML)"),
    LLV = c(LL.dum.ML, LL.dum41), EDF = c(EDF.dum.ML, EDF.dum41),
    deviance = c(m.lcsc.dum.ML$deviance, m.dum41$deviance),
    R2_D = c(deviance_r2(m.lcsc.dum.ML$deviance, m.lcsc.dum.ML$null.deviance), deviance_r2(m.dum41$deviance, m.dum41$null.deviance)),
    Adj_R2_D = c(adj_deviance_r2(m.lcsc.dum.ML$deviance, m.lcsc.dum.ML$null.deviance, n, EDF.dum.ML),
                 adj_deviance_r2(m.dum41$deviance, m.dum41$null.deviance, n, EDF.dum41)),
    QAIC_table_c_hat = c(ic.dum.ML["QAIC"], ic.dum41.tab["QAIC"]), QBIC_table_c_hat = c(ic.dum.ML["QBIC"], ic.dum41.tab["QBIC"]),
    c_hat_own = c(c_hat, c_hat41),
    QAIC_own_c_hat = c(ic.dum.own["QAIC"], ic.dum41.own["QAIC"]), QBIC_own_c_hat = c(ic.dum.own["QBIC"], ic.dum41.own["QBIC"]),
    stringsAsFactors = FALSE)
  cat("\nRobustness: adding 72 period-cohort dummies for 1941-1942 (ages 15-50), ML fits:\n")
  print(format(rob, digits = 8), row.names = FALSE)
  mn_add("mort_rob_n_dummies_4142", length(vichy), 0, "footnote: 72 additional dummies for 1941-1942")
  mn_add("mort_rob_dQAIC_4142", ic.dum41.tab["QAIC"] - ic.dum.ML["QAIC"], 1, "footnote: change in QAIC (table nu-hat) from adding the 1941-42 dummies")
  mn_add("mort_rob_dQBIC_4142", ic.dum41.tab["QBIC"] - ic.dum.ML["QBIC"], 1, "footnote: change in QBIC (table nu-hat) from adding the 1941-42 dummies")
  rm(rs.dum41.ML, G.dum41, df41, data41); invisible(gc())
}

## (6d) CONCURVITY (footnote in Section V.4)
conc.ML.full  <- concurvity_lcsc(m.lcsc.ML, rs.lcsc.ML$X, full = TRUE)
conc.ML.pair  <- concurvity_lcsc(m.lcsc.ML, rs.lcsc.ML$X, full = FALSE)
conc.dum.full <- concurvity_lcsc(m.lcsc.dum.ML, rs.dum.ML$X, full = TRUE)
conc.dum.pair <- concurvity_lcsc(m.lcsc.dum.ML, rs.dum.ML$X, full = FALSE)

cat("\nConcurvity, LC-SC model (ML) -- full (multivariate):\n"); pr(conc.ML.full, 4)
cat("\nConcurvity, LC-SC model (ML) -- pairwise:\n"); print(lapply(conc.ML.pair, round, 4))
cat("\nConcurvity, LC-SC model with varying period shocks (ML) -- full (multivariate):\n"); pr(conc.dum.full, 4)
cat("\nConcurvity, LC-SC model with varying period shocks (ML) -- pairwise:\n"); print(lapply(conc.dum.pair, round, 4))
cat("\nMaximum observed/estimate multivariate concurvity over the smooth terms:",
    round(max(conc.ML.full[c("observed", "estimate"), -1]), 4), "(LC-SC),",
    round(max(conc.dum.full[c("observed", "estimate"), -1]), 4), "(LC-SC + war dummies)\n")

for (msr in c("worst", "observed", "estimate")) for (tm in c("s(a)", "s(p)", "s(c)")) {
  mn_add(paste0("mort_conc_lcsc_", msr, "_", gsub("[()]", "", tm)), conc.ML.full[msr, tm], 3,
         "Appendix F: multivariate concurvity, LC-SC model (ML), fitted design")
  mn_add(paste0("mort_conc_dum_", msr, "_", gsub("[()]", "", tm)), conc.dum.full[msr, tm], 3,
         "Appendix F: multivariate concurvity, model with varying period shocks (ML)")
}

# the same measures on an explicit partition of the columns
conc_blocks_df <- function(rs, model) {
  b <- rs$mod
  bl <- lcsc_blocks(b, rs$X, extra = list(lc_slope = "^a$", sc_slope = "^c$", war = "^a[0-9]+_p[0-9]+$"))
  cb <- concurvity_blocks(rs$X, coef(b), bl, full = TRUE)
  cat("\nConcurvity (explicit blocks),", model, ":\n"); pr(cb, 4)
  do.call(rbind, lapply(colnames(cb), function(tm) data.frame(
    model = model, metric = "unweighted", type = "full", measure = rownames(cb), term = tm,
    term_given = NA_character_, value = cb[, tm], ncol = attr(cb, "ncol")[[tm]],
    rank = attr(cb, "rank")[[tm]], stringsAsFactors = FALSE)))
}
cb.df <- rbind(conc_blocks_df(rs.lcsc.ML, "LC-SC (ML)"), conc_blocks_df(rs.dum.ML, "LC-SC + varying period shocks (ML)"))
stopifnot(max(abs(cb.df$value[cb.df$model == "LC-SC (ML)" & cb.df$term == "s(c)"] - conc.ML.full[, "s(c)"])) < 1e-8)
write.csv(cb.df, file.path("Output", "mortality_concurvity_blocks.csv"), row.names = FALSE)

# the ML fits are not used below
rm(rs.lcsc.ML, rs.dum.ML); invisible(gc())

#### 7. Model-Based Predicted-Rate Matrices (LC-SC Model with War Dummies, REML) ####

## (7a) FULL PREDICTIONS FROM THE LC-SC MODEL WITH DUMMIES

newdata <- data
H <- lcsc_predict(rs.dum, newdata = newdata)
H_est <- H[["est"]] # term-wise contributions (predict = "terms")

# summing the term contributions, adding the intercept and the log-population offset
stopifnot(isTRUE(all.equal(H$intercept + rowSums(H_est), H$lp)))
newdata$lin_pred <- H$intercept + rowSums(H_est) + log(newdata$pop)

# exponentiating to get the predicted counts, then converting to rates per 1,000
newdata$yhat <- exp(newdata$lin_pred)
newdata$yhat.rate <- (newdata$yhat / newdata$pop) * 1000

# sanity check: the predictions reproduce mgcv's fitted values
stopifnot(isTRUE(all.equal(newdata$yhat, as.vector(fitted(m.lcsc.dum)), tolerance = 1e-6)))

# predicted-rate matrices
yhat.lcsc.dum <- newdata[, c("a.index", "p.index", "c.index", "yhat.rate")]
matAP.lcsc.dum <- mean_by_ap(yhat.lcsc.dum); matPC.lcsc.dum <- mean_by_pc(yhat.lcsc.dum); matAC.lcsc.dum <- mean_by_ac(yhat.lcsc.dum)

# coefficient blocks of the figures model, located by name (the embedded helpers)
ix   <- lcsc_coef_index(m.lcsc.dum)
i.mu <- ix[["(Intercept)"]]; i.a <- ix[["a"]]; i.c <- ix[["c"]]
i.sa <- ix[["s(a)"]];        i.sp <- ix[["s(p)"]]; i.sc <- ix[["s(c)"]]
i.war <- match(war.dummies, names(coef(m.lcsc.dum)))
stopifnot(!anyNA(i.war), length(i.war) == 504)

# mean exposure, used to convert the linear predictor of single components into rates
log.mean.pop <- log(mean(data$pop))

# predicted rate per 1,000 from a set of columns of the design matrix
component_rate <- function(cols, sort_col) {
  relevant_X <- H$X[, cols, drop = FALSE]
  relevant_X <- unique(relevant_X)
  relevant_X <- relevant_X[order(relevant_X[, sort_col]), , drop = FALSE]
  linear_predictor <- relevant_X %*% coef(m.lcsc.dum)[cols] + log.mean.pop
  as.vector((exp(linear_predictor) / mean(data$pop)) * 1000)
}

## (7b) PERIOD NONLINEARITIES + WAR-SHOCK (PERIOD-COHORT INTERACTION) TERMS ONLY

relevant_indices <- c(i.mu, i.war, i.sp)
relevant_X <- H$X[, relevant_indices]
relevant_beta <- coef(m.lcsc.dum)[relevant_indices]

# computing predicted rates on a copy of the dataset
newdata <- data
newdata$linear_predictor <- relevant_X %*% relevant_beta + log(newdata$pop)
newdata$predicted_rate <- (exp(newdata$linear_predictor) / newdata$pop) * 1000

# creating the relevant matrices
yhat.interact <- newdata[, c("a.index", "p.index", "c.index", "predicted_rate")]
matAP.interact <- mean_by_ap(yhat.interact)
matAC.interact <- mean_by_ac(yhat.interact)
matPC.interact <- mean_by_pc(yhat.interact)

## (7c) SMOOTHED PREDICTIONS WITHOUT SHOCKS

relevant_indices <- c(i.mu, i.a, i.c, i.sa, i.sc)
relevant_X <- H$X[, relevant_indices]
relevant_beta <- coef(m.lcsc.dum)[relevant_indices]

# computing predicted rates on a copy of the dataset
newdata <- data
newdata$linear_predictor <- relevant_X %*% relevant_beta + log(newdata$pop)
newdata$predicted_rate <- (exp(newdata$linear_predictor) / newdata$pop) * 1000

# creating the relevant matrices
yhat.smooth <- newdata[, c("a.index", "p.index", "c.index", "predicted_rate")]
matAP.smooth <- mean_by_ap(yhat.smooth)
matAC.smooth <- mean_by_ac(yhat.smooth)
matPC.smooth <- mean_by_pc(yhat.smooth)

#### 9. Figure 2: Comparative Cohort Careers, Raw ASMRs (Main Text) ####

# specifying name and size of the plot
pdf(file.path(fig.dir, "CohortCareers_CellRawASMR.pdf"), width = 14.75, height = 11.75) # dimensions in inches

# selecting the color palette (greys)
mypal <- colorRampPalette(brewer.pal(9, "Greys"), alpha = FALSE, bias = 1.75)(300)

# single plot
par(mfrow = c(1, 1))

# age-cohort matrix of raw cell means (cohorts in rows after transposing)
d <- matAC.raw
d <- t(d) # transposing the matrix

# setting up a blank plot
plot(x = a.index, y = a.index, col = "black",
     type = "n", ylim = c(0, 500), xlim = c(0, 100),
     xlab = "Age", ylab = "ASMR per 1,000 Males", yaxt = "n", xaxt = "n",
     main = " ", cex.main = 0.9)

# ticks for the graph
axis(side = 2, at = seq(from = 0, to = 500, by = 50), las = 1, cex.axis = 0.8)
axis(side = 1, at = seq(from = 0, to = 100, by = 20), las = 0, cex.axis = 0.8)

# adding lines for every 5th cohort (excluding the first and last few sparse cohorts)
for (i in 5:(length(c.index) - 2)) {
  if ((i - 5) %% 5 == 0) {
    # for lines and text of different colors
    lines(x = a.index, d[i, ], col = mypal[i],
          lty = 1, lwd = 1)
    x.val <- as.numeric(names(na.omit(d[i, ]))[length(names(na.omit(d[i, ])))])
    if (as.numeric(rownames(d)[i]) < 1945) {
      text(x = x.val + 1,
           y = na.omit(d[i, ])[length(na.omit(d[i, ]))], labels = rownames(d)[i],
           cex = 0.5, col = mypal[i])
    }
  }
}

# red coloring for the 1895 cohort
temp <- d[rownames(d) == 1895, ]
lines(x = a.index, temp, lty = 1, lwd = 1.15, col = "red")
text(x = max(a.index) + 1, y = temp[length(temp)], labels = "1895", cex = 0.5, col = "red")

# closing graphics device
invisible(dev.off())

#### 10. Figure 3: LC Slope and SC Slope (Appendix E) ####

pdf(file.path(fig.dir, "TwoPanel_LCSCSlopes.pdf"), width = 13, height = 6) # dimensions in inches
par(mfrow = c(1, 2))

### PLOT 1: AGE SLOPE (intercept and age linear component)
new_data <- data.frame(a.index = a.index, pred = component_rate(c(i.mu, i.a), sort_col = 2))

plot(x = new_data$a.index, y = new_data$pred, type = "l",
     ylim = c(0, 70), xlab = "Age", ylab = "Age-Specific Mortality Rate per 1,000 Males",
     yaxt = "n", xaxt = "n", main = "(a) LC Slope", cex.main = 1.4)
axis(side = 2, at = seq(from = 0, to = 2000, by = 10), las = 1, cex.axis = 0.8)
axis(side = 1, at = seq(from = 0, to = 110, by = 10), las = 0, cex.axis = 0.8)

### PLOT 2: COHORT SLOPE (intercept and cohort linear component)
new_data <- data.frame(c.index = c.index, pred = component_rate(c(i.mu, i.c), sort_col = 2))

plot(x = new_data$c.index, y = new_data$pred, type = "l", xlim = c(1700, 2030),
     ylim = c(0, 70), xlab = "Cohort", ylab = " ",
     yaxt = "n", xaxt = "n", main = "(b) SC Slope", cex.main = 1.4)
axis(side = 2, at = seq(from = 0, to = 5000, by = 10), las = 1, cex.axis = 0.8)
axis(side = 1, at = seq(from = 1720, to = 2020, by = 30), las = 0, cex.axis = 0.8)

# closing graphics device
invisible(dev.off())

#### 11. Figure 4: LC Curve and SC Curve (Appendix E) ####

pdf(file.path(fig.dir, "TwoPanel_LCSCCurves.pdf"), width = 13, height = 6) # dimensions in inches
par(mfrow = c(1, 2))

### PLOT 1: AGE CURVE (intercept, age linear component, s(a))
new_data <- data.frame(a.index = a.index, pred = component_rate(c(i.mu, i.a, i.sa), sort_col = 2))

plot(x = new_data$a.index, y = new_data$pred, type = "l",
     ylim = c(0, 150), xlab = "Age", ylab = "Age-Specific Mortality Rate per 1,000 Males",
     yaxt = "n", xaxt = "n", main = "(a) LC Curve", cex.main = 1.4)
axis(side = 2, at = seq(from = 0, to = 2000, by = 20), las = 1, cex.axis = 0.8)
axis(side = 1, at = seq(from = 0, to = 110, by = 10), las = 0, cex.axis = 0.8)

### PLOT 2: COHORT CURVE (intercept, cohort linear component, s(c))
new_data <- data.frame(c.index = c.index, pred = component_rate(c(i.mu, i.c, i.sc), sort_col = 2))

plot(x = new_data$c.index, y = new_data$pred, type = "l", xlim = c(1700, 2030),
     ylim = c(0, 40), xlab = "Cohort", ylab = " ",
     yaxt = "n", xaxt = "n", main = "(b) SC Curve", cex.main = 1.4)
axis(side = 2, at = seq(from = 0, to = 5000, by = 5), las = 1, cex.axis = 0.8)
axis(side = 1, at = seq(from = 1720, to = 2020, by = 30), las = 0, cex.axis = 0.8)

# closing graphics device
invisible(dev.off())

#### 12. Figure 5: Period Nonlinearities, Overall Average (Appendix E) ####

# specifying name and size of the plot
pdf(file.path(fig.dir, "Period_NonlinearitiesOverall_Average.pdf"), width = 15.5, height = 8.25) # dimensions in inches
par(mfrow = c(1, 1))

# organizing into a miniature dataset
new_data <- data.frame(p.index = p.index, pred = colMeans(matAP.interact, na.rm = TRUE))

# Retrieve the red color from the Set1 palette
set1_red <- brewer.pal(8, "Set1")[1]

# Create the plot
plot(x = new_data$p.index, y = new_data$pred, type = "l",
     ylim = c(0, 35), xlab = "Period", ylab = "Age-Specific Mortality Rate per 1,000 Males",
     yaxt = "n", xaxt = "n", main = " ", cex.main = 0.9)

# Add shaded regions (war years) using rect() function
rect(xleft = 1852.5, ybottom = par("usr")[3], xright = 1856.5, ytop = par("usr")[4],
     border = NA, col = adjustcolor(set1_red, alpha.f = 0.1))
rect(xleft = 1869.5, ybottom = par("usr")[3], xright = 1871.5, ytop = par("usr")[4],
     border = NA, col = adjustcolor(set1_red, alpha.f = 0.1))
rect(xleft = 1914, ybottom = par("usr")[3], xright = 1918.5, ytop = par("usr")[4],
     border = NA, col = adjustcolor(set1_red, alpha.f = 0.1))
rect(xleft = 1938.5, ybottom = par("usr")[3], xright = 1945.5, ytop = par("usr")[4],
     border = NA, col = adjustcolor(set1_red, alpha.f = 0.1))

# Add axis ticks
axis(side = 2, at = seq(from = 0, to = 300, by = 5), las = 1, cex.axis = 0.8)
axis(side = 1, at = seq(from = 1810, to = 2025, by = 10), las = 0, cex.axis = 0.8)

# Create an array of historical events
events <- c(1821, 1830, 1832.5, 1848, 1889, 1919, 1938.75, 1958, 1962, 1968, 1981, 1992, 2017,
            1944.75, 2003, 2020, 1894)
event_labels <- c("Death of \nNapoleon", "July \nRevolution", "June \nRebellion",
                  "French \nRevolution \nof 1848",
                  "Eiffel Tower Completed",
                  "Treaty of Versailles",
                  "Vichy France \n (1940-1944)", "Establishment \n of the 5th Republic",
                  "End of the Algerian War \n (1954-1962)", "'68 Protests",
                  "Mitterrand Elected", "Treaty of Maastricht", "Macron Elected",
                  "Invasion of \nNormandy", "2003 Heat Wave", "COVID-19 \n Pandemic",
                  "Start of the \n Dreyfus Affair")

# Calculate y-positions for labels and vertical lines
event_y_positions <- sapply(events, function(x) {
  # Find the closest year in the data
  closest_year_index <- which.min(abs(new_data$p.index - x))

  # Return the corresponding mortality rate
  new_data$pred[closest_year_index]
})

# Indicate which events are labeled above (rather than below) the line
events_above <- c("French \nRevolution \nof 1848", "Crimean War \n(1853-1856)",
                  "2003 Heat Wave", 'COVID-19 \n Pandemic',
                  "End of the Algerian War \n (1954-1962)",
                  "June \nRebellion", "Start of the \n Dreyfus Affair")

# Distance of labels from the line
distance_from_line <- 4

# Add labels and vertical lines
for (i in seq_along(events)) {
  # Determine label position (above or below the line)
  if (event_labels[i] %in% events_above) {
    position <- 1
  } else {
    position <- -1
  }

  # Add a label above or below the line
  text(x = events[i], y = event_y_positions[i] + distance_from_line * position,
       labels = event_labels[i], cex = 0.4, adj = c(0.5, 0))

  # Add a vertical dashed line from the label to the line
  segments(x0 = events[i], y0 = event_y_positions[i], x1 = events[i],
           y1 = event_y_positions[i] + distance_from_line * position,
           lty = "dashed", col = adjustcolor("gray", alpha.f = 0.9))
}

# add labels for the major wars
text(x = 1854.5, y = 30, labels = "Crimean War \n (1853-1856)", cex = 0.4, adj = c(0.5, 0))
text(x = 1870.5, y = 30, labels = "Franco-Prussian War \n (1870-1871)", cex = 0.4, adj = c(0.5, 0))
text(x = 1916, y = 30, labels = "WW1 \n(1914-1918)", cex = 0.4, adj = c(0.5, 0))
text(x = 1942, y = 30, labels = "WW2 \n(1939-1945)", cex = 0.4, adj = c(0.5, 0))

# closing the graphics device and reverting to defaults
invisible(dev.off())

#### 13. Figures 6 and 7: Period Nonlinearities at Ages 20 and 40 (Appendix E) ####

period_nonlinearities_at_age <- function(age, file) {
  pdf(file.path(fig.dir, file), width = 15.5, height = 8.25) # dimensions in inches
  par(mfrow = c(1, 1))

  # organizing into a miniature dataset (row of the age-period matrix for this age)
  new_data <- data.frame(p.index = p.index, pred = matAP.interact[as.character(age), ])

  # Retrieve the red color from the Set1 palette
  set1_red <- brewer.pal(8, "Set1")[1]

  # Create the plot
  plot(x = new_data$p.index, y = new_data$pred, type = "l",
       ylim = c(0, 350), xlab = "Period", ylab = "Age-Specific Mortality Rate per 1,000 Males",
       yaxt = "n", xaxt = "n", main = " ", cex.main = 0.9)

  # Add shaded regions (war years) using rect() function
  rect(xleft = 1852.5, ybottom = par("usr")[3], xright = 1856.5, ytop = par("usr")[4],
       border = NA, col = adjustcolor(set1_red, alpha.f = 0.1))
  rect(xleft = 1869.5, ybottom = par("usr")[3], xright = 1871.5, ytop = par("usr")[4],
       border = NA, col = adjustcolor(set1_red, alpha.f = 0.1))
  rect(xleft = 1914, ybottom = par("usr")[3], xright = 1918.5, ytop = par("usr")[4],
       border = NA, col = adjustcolor(set1_red, alpha.f = 0.1))
  rect(xleft = 1938.5, ybottom = par("usr")[3], xright = 1945.5, ytop = par("usr")[4],
       border = NA, col = adjustcolor(set1_red, alpha.f = 0.1))

  # Add axis ticks
  axis(side = 2, at = seq(from = 0, to = 400, by = 50), las = 1, cex.axis = 0.8)
  axis(side = 1, at = seq(from = 1810, to = 2025, by = 10), las = 0, cex.axis = 0.8)

  # add labels for the major wars
  text(x = 1854.5, y = 325, labels = "Crimean War \n (1853-1856)", cex = 0.4, adj = c(0.5, 0))
  text(x = 1870.5, y = 325, labels = "Franco-Prussian War \n (1870-1871)", cex = 0.4, adj = c(0.5, 0))
  text(x = 1916, y = 325, labels = "WW1 \n(1914-1918)", cex = 0.4, adj = c(0.5, 0))
  text(x = 1942, y = 325, labels = "WW2 \n(1939-1945)", cex = 0.4, adj = c(0.5, 0))

  # closing the graphics device
  invisible(dev.off())
}
period_nonlinearities_at_age(20, "Period_Nonlinearities20.pdf")
period_nonlinearities_at_age(40, "Period_Nonlinearities40.pdf")

#### 14. Figure 8: Age-Cohort Ridgeline Plot of War Effects (Appendix E) ####

# age-cohort matrix of period nonlinearities + war effects
d <- matAC.interact
cohorts <- as.character(1870:1930) # convert cohorts to character for matching column names
d_subset <- d[0:90, cohorts]

# Convert the matrix to a data frame
ridge_df <- as.data.frame(d_subset)

# Add a row ID (age) to the data frame
ridge_df$row_id <- as.numeric(as.character(rownames(ridge_df)))

# Gather the columns into key-value pairs
df_long <- ridge_df %>%
  gather(year, value, -row_id) %>%
  mutate(year = as.numeric(year))

# Keep every 5th cohort
df_filtered <- df_long[df_long$year %% 5 == 0, ]

# Output to PDF
pdf(file.path(fig.dir, "AC_Ridges_WarEffects.pdf"), width = 8, height = 5.5) # dimensions in inches

# Create the ggridges plot
print(
  ggplot(df_filtered, aes(x = row_id, y = as.factor(year), height = value)) +
    geom_density_ridges(
      stat = "identity",
      alpha = 0.35, # Transparency
      colour = "darkgray", # Line color
      scale = 3.75
    ) +
    scale_x_continuous(breaks = c(0, 20, 40, 60, 80)) + # Ticks at ages 0, 20, 40, 60, 80
    theme(
      panel.background = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "none" # Remove legend
    ) +
    labs(
      title = " ",
      x = "Age",
      y = "Cohort"
    )
)

# Close the PDF device
invisible(dev.off())

#### 15. Figure 9: Comparative Cohort Careers, Smoothed Yhats (Appendix E) ####

# specifying name and size of the plot
pdf(file.path(fig.dir, "CohortCareers_CellASMR_YhatSmoothed.pdf"), width = 14.75, height = 11.75) # dimensions in inches

# selecting the color palette (greys)
mypal <- colorRampPalette(brewer.pal(9, "Greys"), alpha = FALSE, bias = 1.75)(300)

# single plot
par(mfrow = c(1, 1))

# age-cohort matrix of smoothed predicted rates (cohorts in rows after transposing)
d <- matAC.smooth
d <- t(d) # transposing the matrix

# setting up a blank plot
plot(x = a.index, y = a.index, col = "black",
     type = "n", ylim = c(0, 400), xlim = c(0, 100),
     xlab = "Age", ylab = "ASMR per 1,000 Males", yaxt = "n", xaxt = "n",
     main = " ", cex.main = 0.9)

# ticks for the graph
axis(side = 2, at = seq(from = 0, to = 500, by = 50), las = 1, cex.axis = 0.8)
axis(side = 1, at = seq(from = 0, to = 100, by = 20), las = 0, cex.axis = 0.8)

# adding lines for every 5th cohort (excluding the first and last few sparse cohorts)
for (i in 5:(length(c.index) - 2)) {
  if ((i - 5) %% 5 == 0) {
    # for lines and text of different colors
    lines(x = a.index, d[i, ], col = mypal[i],
          lty = 1, lwd = 1)
    x.val <- as.numeric(names(na.omit(d[i, ]))[length(names(na.omit(d[i, ])))])
    if (as.numeric(rownames(d)[i]) < 1945) {
      text(x = x.val + 1,
           y = na.omit(d[i, ])[length(na.omit(d[i, ]))], labels = rownames(d)[i],
           cex = 0.5, col = mypal[i])
    }
  }
}

# closing graphics device
invisible(dev.off())

#### 16. Figure 10: Cohort Careers for Selected WW1 Cohorts, Grid (Main Text) ####

# specifying name and size of the plot
pdf(file.path(fig.dir, "CohortCareers_CellASMR_Yhat_WW1Grid.pdf"), width = 13.75, height = 4.75) # dimensions in inches

# 1x3 grid of plots
par(mfrow = c(1, 3))

# red color for the war bands
set1_red <- brewer.pal(8, "Set1")[1]

# age-cohort matrix of predicted rates (cohorts in rows after transposing)
d <- matAC.lcsc.dum
d <- t(d) # transposing the matrix

# Define cohorts of interest
cohorts_of_interest <- c("1880", "1895", "1905")

# Subset the data to include only cohorts of interest
d <- d[rownames(d) %in% cohorts_of_interest, ]

# Loop over each cohort
for (i in 1:nrow(d)) {

  # setting up a blank plot
  plot(x = a.index, y = a.index, col = "black",
       type = "n", ylim = c(0, 350), xlim = c(0, 90),
       xlab = "Age", ylab = "ASMR per 1,000 Males", yaxt = "n", xaxt = "n",
       main = paste0("(", letters[i], ") Cohort ", rownames(d)[i]), cex.main = 0.9)

  # ticks for the graph
  axis(side = 2, at = seq(from = 0, to = 500, by = 50), las = 1, cex.axis = 0.8)
  axis(side = 1, at = seq(from = 0, to = 100, by = 10), las = 0, cex.axis = 0.8)

  # Plot lines for each cohort group
  lines(x = a.index, d[i, ], col = "black", lty = 1, lwd = 1.2)

  # WW1 band (at the ages this cohort passed through 1914-1918)
  rect(xleft = 1914 - as.numeric(cohorts_of_interest[i]), ybottom = par("usr")[3],
       xright = 1918.5 - as.numeric(cohorts_of_interest[i]), ytop = par("usr")[4],
       border = NA, col = adjustcolor(set1_red, alpha.f = 0.1))

  # WW2 band (at the ages this cohort passed through 1939-1945)
  rect(xleft = 1938.5 - as.numeric(cohorts_of_interest[i]), ybottom = par("usr")[3],
       xright = 1945.5 - as.numeric(cohorts_of_interest[i]), ytop = par("usr")[4],
       border = NA, col = adjustcolor(set1_red, alpha.f = 0.1))
}

# closing graphics device
invisible(dev.off())

#### 17. Wrap-Up ####

cat("\nFigures written to", fig.dir, ":\n")
print(list.files(fig.dir, pattern = "\\.pdf$"))
cat("Tables written to Output/: mortality_fit_table.csv, mortality_concurvity_blocks.csv, tex/mortality_fit_table.tex\n")
cat("02_mortality_france.R finished:", format(Sys.time()), "\n")

## END OF R CODE
