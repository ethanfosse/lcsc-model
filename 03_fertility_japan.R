# 03_fertility_japan.R -- Japanese fertility, 1947-2020
# Fosse and Winship, "Varieties of Cross-Cohort Differentiation" (Sociological Science).
# Input: Data/fertility.RData (not included; see README.md). Packages: mgcv, plot3D, RColorBrewer.
# Run from this folder: Rscript --vanilla 03_fertility_japan.R
# Writes the paper's figures to Figures/ and tables to Output/.

CRITERION_FIGURES <- "REML"  # smoothing criterion of the figure models (paper: REML)
CRITERION_TABLES  <- "ML"    # smoothing criterion of the table models (paper: ML)
LEGACY_UBRE       <- FALSE
NTHREADS          <- 8       # threads for mgcv::gam.control()
REFIT             <- FALSE
K_COHORT          <- 15      # basis dimension of each cohort-specific curve delta_u(a)
SP_START_FS       <- 0.01
CONCURVITY_METRICS <- c("unweighted", "weighted")

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

lcsc_session_header("03_fertility_japan.R")
lcsc_require_data("Data/fertility.RData", "HFD Japan files")

# loading libraries
suppressPackageStartupMessages({
  library(mgcv)         # penalized-spline GAMs (Simon Wood)
  library(plot3D)       # 3-D histogram and 2-D image/contour plots
  library(RColorBrewer) # color palettes
})

options(warn = 1)  # print warnings where they occur
cat("Settings: CRITERION_FIGURES =", CRITERION_FIGURES, "| CRITERION_TABLES =", CRITERION_TABLES,
    "| LEGACY_UBRE =", LEGACY_UBRE, "| NTHREADS =", NTHREADS, "| REFIT =", REFIT,
    "| K_COHORT =", K_COHORT, "| SP_START_FS =", if (is.null(SP_START_FS)) "NULL" else SP_START_FS, "\n\n")
t.start <- Sys.time()

# Some figures place labels with jitter()
set.seed(20240617)

# output directories
fig.dir <- file.path("Figures", "Fertility")
dir.create(fig.dir, recursive = TRUE, showWarnings = FALSE)
dir.create("Output", showWarnings = FALSE)
dir.create("Models", showWarnings = FALSE)

#### 2. Data Preparation ####

# Human Fertility Database, Japan
load("Data/fertility.RData") # object: df

# listwise deletion (open-ended/incomplete cells)
df <- na.omit(df)

# age, period, and cohort on their natural scales
df$p <- df$period; df$a <- df$age
df$c <- df$p - df$a

# index (label) copies of the temporal variables
df$a.index <- df$a; df$p.index <- df$p; df$c.index <- df$c

# outcome: births in an age-period cell, rounded to an integer for the Poisson models
df$y <- as.integer(round(df$Total))
df$pop <- df$Exposure

# observed ASFR per 1,000 women (births / exposure * 1000)
df$yhat.raw <- df$Total / df$Exposure * 1000

# index vectors (used for plotting labels and prediction frames)
a.index <- as.numeric(as.character(sort(unique(df$a.index)))) # ages 12-55
p.index <- as.numeric(as.character(sort(unique(df$p.index)))) # periods 1947-2020
c.index <- as.numeric(as.character(sort(unique(df$c.index)))) # cohorts 1892-2008

# ---- CONFIGURATION: the only place the analysis constants are written --------------- #
cfg <- lcsc_config(df,
                   age_center    = 34,
                   period_center = 1984,
                   cohort_center = 1950,
                   period_full   = 1947:2020,
                   knots         = c(a = 44, p = 74, c = 117))

# centered modeling columns a, p, c (overwriting the natural-scale ones)
df <- lcsc_center(df, cfg)

# cohort as a factor (levels = the 117 birth cohorts) for the factor-smooth interaction
df$c.cat <- factor(df$c.index)

# fire-horse dummies (create_dummy_vars in the embedded helpers)
rs <- create_dummy_vars(df, "a.index", "p.index", 18:35, 1966:1966)
df <- rs$df; firehorse <- rs$dummy_var_names
stopifnot(length(firehorse) == 18, all(colSums(df[, firehorse]) == 1))

# estimation data set (all models)
data <- df[, c("a", "p", "c", "c.cat", "y", "pop", "a.index", "p.index", "c.index", firehorse)]

n <- nrow(df)
cat("Analysis sample:", n, "age-period cells;", length(a.index), "ages,", length(p.index), "periods,",
    length(c.index), "cohorts;", length(firehorse), "fire-horse dummies (ages 18-35 x 1966)\n")
cat("Variance-to-mean ratio of y:", round(var(df$y) / mean(df$y), 2), "\n\n")

#### 3. Helper Functions ####

# axis label positions of the heat map, derived from the index ranges
age.lab <- lcsc_axis_offsets(a.index, by = 5, origin_adj = 1)
per.lab <- lcsc_axis_offsets(p.index, by = 5, origin_adj = 1)
coh.lab <- lcsc_axis_offsets(c.index, by = 5, origin_adj = 1)

# two-dimensional Lexis-surface heat map with contour lines and age, period
plotAPCHeatmap <- function(d, z.min = NULL, z.max = NULL, by.z = NULL, save_pdf = TRUE, pdf_name = "2D_APC_heatmap.pdf") {
  mypal <- rev(colorRampPalette(brewer.pal(9, "RdYlBu"), alpha = 1, bias = 1)(150))
  if (is.null(z.min)) z.min <- round(min(d, na.rm = TRUE), digits = 0)
  if (is.null(z.max)) z.max <- round(max(d, na.rm = TRUE), digits = 0)
  if (is.null(by.z)) by.z <- round(diff(seq(from = z.min, to = z.max, length.out = 15))[1], digits = 1)
  if (save_pdf) pdf(pdf_name, width = 11.25, height = 9.5)
  op <- par(mar = c(5.4, 4.1, 4.1, 5.1))
  d <- d[nrow(d):1, ]; d <- t(d); d.contour <- d
  age.plot.labs <- age.lab$labs;    age.loc <- age.lab$loc
  period.plot.labs <- per.lab$labs; period.loc <- per.lab$loc
  cohort.plot.labs <- coh.lab$labs
  plot3D::image2D(z = d, x = 1:nrow(d), y = 1:ncol(d), shade = 0.01, rasterImage = FALSE,
                  col = mypal, colkey = FALSE, axes = F, ylab = "", xlab = "")
  plot3D::contour2D(z = d.contour, x = 1:nrow(d.contour), y = 1:ncol(d.contour), col = "black", labcex = 0.5, lwd = 1, alpha = 0.8,
                    levels = seq(z.min, z.max, b = by.z), add = TRUE)
  text(x = 0.5, y = age.loc, pos = 2, srt = 0, labels = rev(age.plot.labs), xpd = TRUE, cex = 0.6)
  text(x = period.loc, y = length(a.index) + 0.5, pos = 3, srt = 0, labels = period.plot.labs, xpd = TRUE, cex = 0.6)
  text(x = period.loc + 0.25, y = -1.00, pos = 3, srt = 320,
       labels = cohort.plot.labs[1:length(period.loc)], xpd = TRUE, cex = 0.6)
  text(x = length(p.index) + 0.65, y = age.loc, pos = 4, srt = 320,
       labels = cohort.plot.labs[length(period.loc) + 1:length(cohort.plot.labs)], xpd = TRUE, cex = 0.6)
  axis(side = 2, at = age.loc, tck = -0.01, labels = F)
  axis(side = 1, at = period.loc, tck = -0.015, labels = F)
  axis(side = 3, at = period.loc, tck = -0.01, labels = F)
  axis(side = 4, at = age.loc, tck = -0.015, labels = F)
  mtext(side = 1, "Cohort", line = 2.25, cex = 0.8)
  mtext(side = 2, "Age", line = 2.25, cex = 0.8)
  mtext(side = 3, "Period", line = 2, cex = 0.8)
  text(x = length(p.index) + 4.5, y = mean(seq(1:length(a.index))), xpd = T, labels = "Cohort", cex = 0.8, srt = 270)
  par(op)
  if (save_pdf) invisible(dev.off())
}

# y-axis limits: the defaults, widened to the data if needed
ylim_or_default <- function(x, default) {
  r <- range(x, na.rm = TRUE)
  if (r[1] >= default[1] && r[2] <= default[2]) return(default)
  pad <- 0.05 * diff(r)
  c(floor(min(r[1] - pad, default[1])), ceiling(max(r[2] + pad, default[2])))
}

# each fit is cached in Models/
lcsc_first_run_notice("fert_", total_min = 30)
fit_cached <- function(name, G, method, formula, fs_exclude = NULL, fs_weights = NULL, in.out = NULL,
                       expected_min = NULL, fs_project = TRUE) {
  fit_or_load(name, G, data = data, cfg = cfg, method = method, refit = REFIT,
              fs_exclude = fs_exclude, fs_weights = fs_weights, in.out = in.out,
              legacy_ubre = LEGACY_UBRE, nthreads = NTHREADS, formula = formula,
              expected_min = expected_min, fs_project = fs_project)
}

#### 4. Raw-Data Figures ####

# observed ASFRs as age-period and age-cohort matrices
yhat.raw <- df[, c("a.index", "p.index", "c.index", "yhat.raw", "pop")]
matAP.raw <- mean_by_ap(yhat.raw); matAC.raw <- mean_by_ac(yhat.raw)

# (4a) OBSERVED-ASFR HEAT MAP (APPENDIX E)
d <- matAP.raw
z.min <- round(min(d, na.rm = TRUE), digits = 1)
z.max <- round(max(d, na.rm = TRUE), digits = 1)
plotAPCHeatmap(d, z.min, z.max, by.z = 20, save_pdf = TRUE,
               pdf_name = file.path(fig.dir, "2D_matAP.raw.pdf"))

# (4b) 3-D HISTOGRAM OF THE RAW ASFR SURFACE (APPENDIX E)
mypal <- rev(colorRampPalette(brewer.pal(9, "RdYlBu"), alpha = 1, bias = 1)(295))
expand <- 0.3
alpha.plane <- 0.95 # transparency of surface
pdf(file.path(fig.dir, "3D_raw.histogram.pdf"), width = 12, height = 12) # dimensions in inches
d <- matAP.raw
par(mfrow = c(1, 1))
theta <- 145; phi <- 25
hist3D(x = as.numeric(rownames(d)), y = as.numeric(colnames(d)), expand = expand,
       z = d, phi = phi, theta = theta, col = mypal, resfac = 1,
       add = F, alpha.plane, colkey = F, shade = 0.5, xlab = "Age", ylab = "Period",
       zlab = "ASFR per \n1,000 Females", zlim = c(0, 450), opaque.top = F, contour = F)
invisible(dev.off())

## (4c) RAW COMPARATIVE COHORT CAREERS WITH THE 1940 COHORT IN RED (MAIN TEXT)
pdf(file.path(fig.dir, "Raw_CohortCareersRed.pdf"), width = 10, height = 7.25) # dimensions in inches
par(mfrow = c(1, 1))
d <- matAC.raw
d <- t(d) # transposing the matrix
c.index.number <- 1:length(c.index)
mypal <- colorRampPalette(brewer.pal(9, "Greys"), alpha = FALSE, bias = 3)(length(c.index))
mypal_transp <- paste0(mypal, "40")
mypal_transp <- rev(mypal_transp) # reversing the color scheme
plot(x = a.index, y = a.index, col = "black",
     type = "n", ylim = c(0, 350), xlim = c(10, 60),
     xlab = "Age", ylab = "Age-Specific Fertility Rate per 1,000 Women", yaxt = "n", xaxt = "n",
     main = " ", cex.main = 0.9)
axis(side = 2, at = seq(from = 0, to = 400, by = 50), las = 1, cex.axis = 0.8)
axis(side = 1, at = seq(from = 10, to = 60, by = 10), las = 0, cex.axis = 0.8)
for (i in 1:length(c.index.number)) {
  c.ind <- c.index.number[i]
  lines(x = a.index, d[c.ind, ], col = mypal_transp[i], lty = 1, lwd = 1)
}
# overlay: the 1940 cohort in red
c.index.number <- which(c.index %in% 1940)
lines(x = a.index, d[c.index.number, ], col = "red", lty = 1, lwd = 1.2)
invisible(dev.off())

#### 5. Basic LC-SC Model (m.lcsc: REML for the figure, ML for the table) ####

# knots for age, period, and cohort
ak <- cfg$knots[["a"]]; pk <- cfg$knots[["p"]]; ck <- cfg$knots[["c"]]

# LC-SC model: y ~ theta1*a + theta2*c + f(a) + g(p) + h(c), Poisson with log-exposure
formula.lcsc <- as.formula(y ~ a + c + s(a, bs = "cr", k = ak - 1) + s(p, bs = "cr", k = pk - 1) + s(c, bs = "cr", k = ck - 1))

# unfitted Poisson prefit (fit = FALSE)
G.lcsc <- mgcv::gam(formula.lcsc, offset = log(pop), family = "poisson", data = data, fit = FALSE)
cat("Basic LC-SC model:", ncol(G.lcsc$X), "coefficients\n")

# (5a) REML FIT -- FailureAdditiveModel figure
rs.lcsc <- fit_cached("fert_m.lcsc", G.lcsc, CRITERION_FIGURES, formula.lcsc)
m.lcsc <- rs.lcsc$mod

# (5b) ML FIT -- fit-statistics table
rs.lcsc.ML <- fit_cached("fert_m.lcsc", G.lcsc, CRITERION_TABLES, formula.lcsc)
m.lcsc.ML <- rs.lcsc.ML$mod
rm(G.lcsc); invisible(gc())

cat("\nBasic LC-SC model -- smoothing criterion actually used (m$method):",
    m.lcsc$method, "(figure) /", m.lcsc.ML$method, "(table)\n")
cat("Smoothing parameters:", signif(m.lcsc$sp, 4), "(figure) /", signif(m.lcsc.ML$sp, 4), "(table)\n")
edf.lcsc <- rbind(c(sapply(seq_along(m.lcsc$smooth), function(i) smooth_edf(m.lcsc, i)), total = sum(m.lcsc$edf)),
                  c(sapply(seq_along(m.lcsc.ML$smooth), function(i) smooth_edf(m.lcsc.ML, i)), total = sum(m.lcsc.ML$edf)))
dimnames(edf.lcsc) <- list(c(paste0("m.lcsc (", m.lcsc$method, ")"), paste0("m.lcsc.ML (", m.lcsc.ML$method, ")")),
                           c(sapply(m.lcsc$smooth, `[[`, "label"), "total"))
cat("Effective degrees of freedom:\n"); pr(edf.lcsc, 2)
cat("logLik:", round(as.numeric(logLik(m.lcsc)), 2), "(figure) /", round(as.numeric(logLik(m.lcsc.ML)), 2), "(table)\n")
cat("\nParametric and smooth-term tables of the figure model (m.lcsc):\n")
print(summary(m.lcsc)$p.table); print(summary(m.lcsc)$s.table)

# (5c) PREDICTED ASFRs FROM THE ADDITIVE LC-SC MODEL (REML fit)
H <- lcsc_predict(rs.lcsc, newdata = data)
stopifnot(isTRUE(all.equal(H$intercept + rowSums(H$est), H$lp)))
stopifnot(isTRUE(all.equal(exp(H$lp + log(data$pop)), as.numeric(fitted(m.lcsc)), tolerance = 1e-6)))
yhat.lcsc <- data.frame(a.index = data$a.index, p.index = data$p.index, c.index = data$c.index,
                        yhat.rate = exp(H$lp) * 1000)
matAC.lcsc <- mean_by_ac(yhat.lcsc)

# period fluctuations g~(p) of the basic LC-SC model (REML), one value per year
g.lcsc <- as.numeric(tapply(H$est[, "s(p)"], data$p.index, mean))
stopifnot(length(g.lcsc) == length(p.index))
gl.pct <- 100 * expm1(g.lcsc)
mn_add("fert_lcsc_gpct_1947", gl.pct[p.index == 1947], 1, "Section V.5: period fluctuation (%) of the basic LC-SC model in 1947")
mn_add("fert_lcsc_gpct_min", min(gl.pct), 1, "Section V.5: lowest period fluctuation (%) of the basic LC-SC model")
mn_add("fert_lcsc_gpct_min_year", as.character(p.index[which.min(gl.pct)]), 0, "Section V.5: year of the lowest period fluctuation of the basic LC-SC model")
mn_add("fert_lcsc_gpct_2020", gl.pct[p.index == 2020], 1, "Section V.5: period fluctuation (%) of the basic LC-SC model in 2020")
cat(sprintf("Basic LC-SC model (REML): period fluctuations %.1f%% (1947), %.1f%% (lowest, %d), %.1f%% (2020)\n",
            gl.pct[p.index == 1947], min(gl.pct), p.index[which.min(gl.pct)], gl.pct[p.index == 2020]))

# (5d) FAILURE OF THE ADDITIVE LC-SC MODEL (APPENDIX E)
pdf(file.path(fig.dir, "FailureAdditiveModel.pdf"), width = 10, height = 7.25) # dimensions in inches
d <- t(matAC.lcsc)
c.index.number <- 1:length(c.index)
mypal <- colorRampPalette(brewer.pal(9, "Blues"), alpha = FALSE, bias = 3)(length(c.index))
mypal_transp <- rev(paste0(mypal, "80"))
plot(x = a.index, y = a.index, col = "black",
     type = "n", ylim = c(0, 350), xlim = c(10, 60),
     xlab = "Age", ylab = "Age-Specific Fertility Rate per 1,000 Women", yaxt = "n", xaxt = "n",
     main = " ", cex.main = 0.9)
axis(side = 2, at = seq(from = 0, to = 500, by = 50), las = 1, cex.axis = 0.8)
axis(side = 1, at = seq(from = 10, to = 60, by = 10), las = 0, cex.axis = 0.8)
for (i in 1:length(c.index.number)) lines(x = a.index, d[c.index.number[i], ], col = mypal_transp[i], lty = 1, lwd = 1.1)
d <- t(matAC.raw)
mypal <- colorRampPalette(brewer.pal(9, "Reds"), alpha = FALSE, bias = 3)(length(c.index))
mypal_transp <- rev(paste0(mypal, "40"))
for (i in 1:length(c.index.number)) lines(x = a.index, d[c.index.number[i], ], col = mypal_transp[i], lty = 1, lwd = 1)
invisible(dev.off())
rm(H); invisible(gc())

#### 6. Varying-LC Models (m.vt, m.vt.fh): REML for the figures, ML for the table ####

cat("\n---------------------------------------------------------------------------\n")
cat("Varying-LC models (plain factor smooth): the LC curve (a + s(a), 44 knots) of the basic\n")
cat(sprintf("LC-SC model plus a factor-smooth interaction s(a, c.cat, bs = 'fs', xt = 'cr', k = %d) with\n", K_COHORT))
cat(sprintf("%d x %d = %d basis columns, NOT projected: each cohort curve has its own level and slope\n",
            length(c.index), K_COHORT, length(c.index) * K_COHORT))
cat("(which absorb the differences between cohorts; there is no SC curve), and the penalties divide\n")
cat("the fitted surface between the overall LC curve and the cohort curves (Pedersen et al.'s GS\n")
cat("model). There is NO period term: the common period fluctuations are the yearly leftovers of\n")
cat("the fitted curves (period_leftovers(), 6c). The preferred model adds the 18 fire-horse\n")
cat("dummies. gam() with REML (figures) and ML (table), as for every other model.\n")
cat("---------------------------------------------------------------------------\n\n")

# (6a) FORMULAS: no SC curve and no period term
rhs.vt <- sprintf("a + s(a, bs = 'cr', k = %d) + s(a, c.cat, bs = 'fs', xt = 'cr', k = %d)", ak - 1, K_COHORT)
formula.vt    <- as.formula(paste("y ~", rhs.vt))
formula.vt.fh <- as.formula(paste("y ~", rhs.vt, "+", paste(firehorse, collapse = " + ")))
cat("Formula, LC-SC model with varying LC curves (m.vt):\n"); print(formula.vt, showEnv = FALSE)
cat("Formula, LC-SC model with varying LC curves and varying period shocks (m.vt.fh):\n")
cat("  y ~", rhs.vt, "+ <18 fire-horse dummies>\n\n")

# unfitted Poisson prefits (fit = FALSE)
G.vt    <- suppressWarnings(mgcv::gam(formula.vt,    offset = log(pop), family = "poisson", data = data, fit = FALSE))
G.vt.fh <- suppressWarnings(mgcv::gam(formula.vt.fh, offset = log(pop), family = "poisson", data = data, fit = FALSE))
cat("Coefficients: m.vt =", ncol(G.vt$X), "; m.vt.fh =", ncol(G.vt.fh$X), "\n")

# starting values of the four smoothing parameters (see SP_START_FS in the header)
sp.start <- if (is.null(SP_START_FS)) NULL else c(m.lcsc$sp[1], rep(SP_START_FS, 3))
in.out <- if (is.null(sp.start)) NULL else list(sp = sp.start, scale = 1)
cat("Starting smoothing parameters:", if (is.null(sp.start)) "mgcv default" else signif(sp.start, 4), "\n\n")

# (6b) FITS (cached)
rs.vt       <- fit_cached("fert_m.vt",    G.vt,    CRITERION_FIGURES, formula.vt,    in.out = in.out, expected_min = 11, fs_project = FALSE)
rs.vt.fh    <- fit_cached("fert_m.vt.fh", G.vt.fh, CRITERION_FIGURES, formula.vt.fh, in.out = in.out, expected_min = 11, fs_project = FALSE)
rs.vt.ML    <- fit_cached("fert_m.vt",    G.vt,    CRITERION_TABLES,  formula.vt,    in.out = in.out, expected_min = 11, fs_project = FALSE)
rs.vt.fh.ML <- fit_cached("fert_m.vt.fh", G.vt.fh, CRITERION_TABLES,  formula.vt.fh, in.out = in.out, expected_min = 11, fs_project = FALSE)
m.vt <- rs.vt$mod; m.vt.fh <- rs.vt.fh$mod; m.vt.ML <- rs.vt.ML$mod; m.vt.fh.ML <- rs.vt.fh.ML$mod
rm(G.vt, G.vt.fh); invisible(gc())

## (6c) THE COMMON PERIOD FLUCTUATIONS: THE YEARLY LEFTOVERS OF EACH FIT (Appendix C.3)
is.fh <- unname(rowSums(data[, firehorse]) == 1)
stopifnot(sum(is.fh) == length(firehorse), all(data$p.index[is.fh] == 1966))
pl.vt       <- period_leftovers(data$y, as.numeric(fitted(m.vt)),       data$p.index)
pl.vt.fh    <- period_leftovers(data$y, as.numeric(fitted(m.vt.fh)),    data$p.index, exclude = is.fh)
pl.vt.ML    <- period_leftovers(data$y, as.numeric(fitted(m.vt.ML)),    data$p.index)
pl.vt.fh.ML <- period_leftovers(data$y, as.numeric(fitted(m.vt.fh.ML)), data$p.index, exclude = is.fh)
# the 74 leftovers are counted as 72 effective degrees of freedom
EDF_PERIOD <- length(p.index) - 2
lev_trend <- function(pl) {
  t <- pl$table; w <- t$expected
  c(mean = sum(w * t$g) / sum(w), slope = unname(coef(lm(g ~ period, data = t, weights = w))[2]))
}

## (6d) SUMMARY OF THE FITS
smooth_i <- function(m, lab) which(sapply(m$smooth, `[[`, "label") == lab)
fs_i     <- function(m) which(sapply(m$smooth, inherits, "fs.interaction"))
fs.label <- rs.vt.fh$ortho[[fs_i(m.vt.fh)]]$label
report_vt <- function(rs, pl, label) {
  m <- rs$mod; H <- lcsc_predict(rs, data)
  chk <- all.equal(as.numeric(fitted(m)), as.numeric(exp(H$lp + log(data$pop))), tolerance = 1e-6)
  if (!isTRUE(chk)) warning(label, ": fitted values and term contributions disagree: ", chk)
  stopifnot(is.null(rs$ortho[[fs_i(m)]]$B))   # the Simple specification: the block is not projected
  lca <- H$est[, "a"] + H$est[, "s(a)"]
  du <- H$est[, fs.label]
  ll.gam <- as.numeric(logLik(m)); ll <- sum(dpois(data$y, pl$mu, log = TRUE))
  cat(sprintf("\n%s\n", label))
  cat("  method:", m$method, "| coefficients:", length(coef(m)), "| criterion score:", round(m$gcv.ubre, 2),
      "| converged:", m$converged, "| fitted in", round(rs$elapsed.min, 1), "min\n")
  cat("  smoothing parameters (s(a); fs wiggliness, ridge 1, ridge 2):", signif(m$sp, 6), "\n")
  cat("  EDF of the GAM: total", round(sum(m$edf), 2), "| s(a)", round(smooth_edf(m, smooth_i(m, "s(a)")), 2),
      "|", fs.label, round(smooth_edf(m, fs_i(m)), 2),
      "| + period leftovers", EDF_PERIOD, "=", round(sum(m$edf) + EDF_PERIOD, 2), "\n")
  cat("  logLik: GAM alone", round(ll.gam, 2), "| with the period leftovers", round(ll, 2),
      "| deviance (with leftovers):", round(safe_deviance(data$y, pl$mu), 2), "\n")
  cat("  intercept =", signif(H$intercept, 5), "| theta1 =", signif(coef(m)[["a"]], 5), "\n")
  cat("  LC(a) = theta1 a + f~(a) at ages 12/20/27/35/45/55:",
      round(tapply(lca, data$a.index, mean)[as.character(c(12, 20, 27, 35, 45, 55))], 3), "\n")
  cat("  range of delta_u(a):", round(range(du), 3), "| range of g(p):", round(range(pl$table$g), 3), "\n")
  # births-weighted means of delta_u(a), weighted by the model's own fitted births
  w.fit <- as.numeric(fitted(m))
  wmean <- function(by) tapply(du * w.fit, by, sum) / tapply(w.fit, by, sum)
  cat("  births-weighted mean of delta_u(a): within cohorts, range", round(range(wmean(data$c.index)), 3),
      "| within ages, max |.|", signif(max(abs(wmean(data$a.index))), 3),
      "| within periods, max |.|", signif(max(abs(wmean(data$p.index))), 3), "\n")
  lt <- lev_trend(pl)
  cat("  period leftovers: expected-count-weighted mean", signif(lt[["mean"]], 3),
      "| weighted linear trend per year", signif(lt[["slope"]], 3), "\n")
  invisible(H)
}
H.vt    <- report_vt(rs.vt,       pl.vt,       "LC-SC model with varying LC curves (m.vt; figure model)")
H.vt.fh <- report_vt(rs.vt.fh,    pl.vt.fh,    "LC-SC model with varying LC curves and varying period shocks (m.vt.fh; figure model)")
H.vt.ML <- report_vt(rs.vt.ML,    pl.vt.ML,    "LC-SC model with varying LC curves (m.vt.ML; table model)")
H.fh.ML <- report_vt(rs.vt.fh.ML, pl.vt.fh.ML, "LC-SC model with varying LC curves and varying period shocks (m.vt.fh.ML; table model)")

# criterion comparison (REML figure fits vs ML table fits), printed below
crit_row <- function(rs, pl, H, model, role) {
  m <- rs$mod
  fa <- tapply(H$est[, "a"] + H$est[, "s(a)"], data$a.index, mean)[as.character(c(12, 20, 27, 35, 45, 55))]
  fh <- if (all(firehorse %in% names(coef(m)))) mean(coef(m)[firehorse]) else NA_real_
  # the three fs smoothing parameters
  lt <- lev_trend(pl)
  data.frame(model = model, role = role, method = m$method, minutes = rs$elapsed.min,
             score = m$gcv.ubre, sp_a = m$sp[1],
             sp_d_wiggle = m$sp[2], sp_d_null1 = m$sp[3], sp_d_null2 = m$sp[4],
             EDF_gam = sum(m$edf), EDF = sum(m$edf) + EDF_PERIOD, EDF_f = smooth_edf(m, smooth_i(m, "s(a)")),
             EDF_d = smooth_edf(m, fs_i(m)), EDF_g = EDF_PERIOD,
             intercept = unname(coef(m)[1]), theta1 = coef(m)[["a"]],
             logLik_gam = as.numeric(logLik(m)), logLik = sum(dpois(data$y, pl$mu, log = TRUE)),
             deviance = safe_deviance(data$y, pl$mu),
             leftover_mean = lt[["mean"]], leftover_trend = lt[["slope"]],
             t(setNames(as.numeric(fa), paste0("LC_a", c(12, 20, 27, 35, 45, 55)))),
             mean_firehorse_coef = fh, check.names = FALSE)
}
crit <- rbind(crit_row(rs.vt, pl.vt, H.vt, "m.vt", "figures"), crit_row(rs.vt.ML, pl.vt.ML, H.vt.ML, "m.vt", "table"),
              crit_row(rs.vt.fh, pl.vt.fh, H.vt.fh, "m.vt.fh", "figures"), crit_row(rs.vt.fh.ML, pl.vt.fh.ML, H.fh.ML, "m.vt.fh", "table"))
cat("\nCriterion comparison of the varying-LC fits:\n")
print(format(crit[, c("model", "role", "method", "minutes", "EDF", "EDF_f", "EDF_d", "logLik", "LC_a27", "mean_firehorse_coef")], digits = 6), row.names = FALSE)
mn_add("fert_fig_edf_vt_fh", sum(m.vt.fh$edf) + EDF_PERIOD, 2, "Section V.5 footnote: EDF of the preferred model under REML (figures), GAM + 72 for the period leftovers")
mn_add("fert_fig_llv_vt_fh", sum(dpois(data$y, pl.vt.fh$mu, log = TRUE)), 2, "Section V.5 footnote: log-likelihood of the preferred model under REML (figures)")
mn_add("fert_minutes_max", ceiling(max(crit$minutes)), 0, "Section V.5: fitting time of the varying-LC models (upper bound, minutes; machine-dependent)")
# the leftovers' weighted linear trend
r.fh <- crit[crit$model == "m.vt.fh" & crit$role == "table", ]
mn_add("fert_leftover_trend_total", 100 * abs(r.fh$leftover_trend) * (max(p.index) - min(p.index)), 2,
       "Appendix C.3: absolute weighted linear trend of the leftovers over 1947-2020, in percent, preferred ML model")

# fire-horse dummy coefficients (log scale) of the preferred figure model
b.fh <- coef(m.vt.fh)[firehorse]
se.fh <- sqrt(diag(m.vt.fh$Vp))[match(firehorse, names(coef(m.vt.fh)))]
cat("\nFire-horse (1966) dummy coefficients of m.vt.fh, ages 18-35 (log scale; departure of each cell from the cohort curves):\n")
print(round(cbind(estimate = b.fh, se = se.fh, exp = exp(b.fh)), 4))
if (smooth_edf(m.vt.fh, 1) < 2) warning("The common LC curve of m.vt.fh has EDF < 2 (nearly linear) -- check the specification.")

#### 7. Common Period Fluctuations g(p) and the 1966 Shock ####
pct <- function(g) 100 * (exp(g) - 1)
g.vtfh <- setNames(pl.vt.fh$table$g, pl.vt.fh$table$period)
g66 <- g.vtfh[["1966"]]
add.fh <- b.fh - g66                           # additional (cohort-varying) 1966 shocks, log scale
fh.mean.log <- mean(add.fh)                    # their equal-age (geometric) average
tot.mean.log <- mean(b.fh)                     # total 1966 departure at ages 18-35 = g(1966) + fh.mean.log
ifh <- match(firehorse, names(coef(m.vt.fh)))
se.tot <- sqrt(sum(m.vt.fh$Vp[ifh, ifh])) / length(firehorse)   # SE of the average dummy coefficient
cat("\nFire-horse (1966): common fluctuation g(1966) and the additional shocks by age (log scale):\n")
print(round(cbind(total = b.fh, common = g66, additional = add.fh, additional_pct = pct(add.fh)), 4))
cat(sprintf("g(1966) = %.4f (%.1f%%); average additional shock %.4f (%.1f%%); total %.4f (%.1f%%)\n",
            g66, pct(g66), fh.mean.log, pct(fh.mean.log), tot.mean.log, pct(tot.mean.log)))

period.tab <- data.frame(
  period = p.index,
  lcsc_g = g.lcsc, lcsc_pct = pct(g.lcsc),                                   # basic LC-SC model, s(p)
  vt_g = pl.vt$table$g, vt_se = pl.vt$table$se, vt_pct = pct(pl.vt$table$g),
  vtfh_g = pl.vt.fh$table$g, vtfh_se = pl.vt.fh$table$se, vtfh_pct = pct(pl.vt.fh$table$g),
  shock1966_log = ifelse(p.index == 1966, fh.mean.log, 0))
period.tab$vtfh_g_with_shock   <- period.tab$vtfh_g + period.tab$shock1966_log
period.tab$vtfh_pct_with_shock <- pct(period.tab$vtfh_g_with_shock)
stopifnot(abs(period.tab$vtfh_g_with_shock[period.tab$period == 1966] - tot.mean.log) < 1e-12)
# pointwise 95% bounds for Figure 7
period.tab$combined_se <- ifelse(p.index == 1966, se.tot, period.tab$vtfh_se)
period.tab$combined_lo <- period.tab$vtfh_g_with_shock - qnorm(.975) * period.tab$combined_se
period.tab$combined_hi <- period.tab$vtfh_g_with_shock + qnorm(.975) * period.tab$combined_se

# pooled changes in the fitted births of 1966
ip66 <- data$p.index == 1966; peak66 <- ip66 & data$a.index %in% 18:35
stopifnot(identical(peak66, is.fh))
mu.fh <- pl.vt.fh$mu                                                        # fitted births, preferred model
lp.curves <- H.vt.fh$intercept + H.vt.fh$est[, "a"] + H.vt.fh$est[, "s(a)"] + H.vt.fh$est[, fs.label]
counterfactual <- exp(lp.curves + log(data$pop))                             # cohort curves only
stopifnot(isTRUE(all.equal(mu.fh[!is.fh], (counterfactual * exp(pl.vt.fh$g_cell))[!is.fh])),
          isTRUE(all.equal(mu.fh[is.fh], as.numeric(data$y[is.fh]), tolerance = 1e-6)))
pooled <- data.frame(ages = c("18-35", "12-55"), pct_change = sapply(list(peak66, ip66), function(ii)
  100 * (sum(mu.fh[ii]) / sum(counterfactual[ii]) - 1)))
mn_add("fert_pooled_18_35", -pooled$pct_change[1], 1, "Section V.5: pooled 1966 reduction, ages 18-35")
mn_add("fert_pooled_all", -pooled$pct_change[2], 1, "Section V.5: pooled 1966 reduction, ages 12-55")
mn_add("fert_outside_peak_share", 100 * sum(data$y[ip66 & !peak66]) / sum(data$y[ip66]), 2,
       "Section V.5: share of 1966 births outside ages 18-35")

cat("\nCommon period fluctuations g(p) (percent), 1960-1972, with the basic LC-SC model's for comparison\n")
print(round(period.tab[period.tab$period %in% 1960:1972,
                       c("period", "lcsc_pct", "vt_pct", "vtfh_pct", "vtfh_se", "vtfh_pct_with_shock")], 3), row.names = FALSE)
ex66 <- period.tab$period != 1966
cat(sprintf("Range of g(p) excluding 1966, preferred model: %.1f%% (%d) to %.1f%% (%d); m.vt: %.1f%% to %.1f%%\n",
            min(period.tab$vtfh_pct[ex66]), p.index[ex66][which.min(period.tab$vtfh_pct[ex66])],
            max(period.tab$vtfh_pct[ex66]), p.index[ex66][which.max(period.tab$vtfh_pct[ex66])],
            min(period.tab$vt_pct[ex66]), max(period.tab$vt_pct[ex66])))
cat(sprintf("Largest conditional standard error of a leftover other than 1966: %.4f (log scale)\n", max(period.tab$vtfh_se[ex66])))

# numbers of the period analysis quoted in the manuscript text (Section V.5)
for (yr in 1965:1967) {
  mn_add(paste0("fert_vtfh_pct_", yr), period.tab$vtfh_pct[period.tab$period == yr], 1,
         "Section V.5: common period fluctuation (%) of the preferred model")
}
mn_add("fert_vtfh_pnon_range_lo", min(period.tab$vtfh_pct[ex66]), 1, "Section V.5: range of the period fluctuations other than 1966 (%)")
mn_add("fert_vtfh_pnon_range_hi", max(period.tab$vtfh_pct[ex66]), 1, "Section V.5: range of the period fluctuations other than 1966 (%)")
mn_add("fert_vtfh_pnon_year_lo", as.character(p.index[ex66][which.min(period.tab$vtfh_pct[ex66])]), 0,
       "Section V.5: year of the lowest period fluctuation other than 1966")
mn_add("fert_vtfh_pnon_year_hi", as.character(p.index[ex66][which.max(period.tab$vtfh_pct[ex66])]), 0,
       "Section V.5: year of the highest period fluctuation other than 1966")
mn_add("fert_vtfh_se_max", 100 * max(period.tab$vtfh_se[ex66]), 2,
       "Figure 7 note (optional): largest conditional Poisson standard error of a yearly leftover other than 1966 (x 100)")
g.log <- g.vtfh
mn_add("fert_1966_rel_neighbors_pct", -100 * expm1(g.log[["1966"]] - (g.log[["1965"]] + g.log[["1967"]]) / 2), 1,
       "Section V.5: common 1966 decline (%) relative to the geometric mean of 1965 and 1967 (sign dropped)")
cat("\nCommon period fluctuations (%) of the preferred model, by year:\n")
print(round(setNames(period.tab$vtfh_pct, period.tab$period), 1))
age.of <- function(nm) as.numeric(sub("a.index([0-9]+)_.*", "\\1", nm))
mn_add("fert_firehorse_avg_pct", 100 * (1 - exp(fh.mean.log)), 1,
       "Section V.5: average additional 1966 shock (% reduction in the ASFR at ages 18-35)")
mn_add("fert_firehorse_max_pct", 100 * (1 - exp(min(add.fh))), 1, "Section V.5: largest additional 1966 reduction (%)")
mn_add("fert_firehorse_max_age", age.of(names(which.min(add.fh))), 0, "Section V.5: age of the largest additional 1966 reduction")
mn_add("fert_firehorse_top_pct", 100 * expm1(max(add.fh)), 1, "Section V.5: largest additional 1966 change toward higher fertility (%, signed)")
mn_add("fert_firehorse_top_age", age.of(names(which.max(add.fh))), 0, "Section V.5: age of the largest additional 1966 change toward higher fertility")
mn_add("fert_firehorse_n_positive", sum(add.fh > 0), 0, "Section V.5 (optional): number of ages with a positive additional 1966 shock")
mn_add("fert_total_1966_pct", -100 * expm1(tot.mean.log), 1,
       "Section V.5: combined 1966 reduction (common fluctuation + average additional shock, %), ages 18-35")
mn_add("fert_1966_age18_combined_pct", 100 * expm1(b.fh[["a.index18_p.index1966"]]), 1,
       "Section V.5: combined 1966 change at age 18 (common fluctuation + age-18 shock, %)")
mn_add("fert_total_1966_min_pct", 100 * (1 - exp(min(b.fh))), 1, "Section V.5 (optional): largest combined 1966 reduction at one age (%)")
mn_add("fert_total_1966_min_age", age.of(names(which.min(b.fh))), 0, "Section V.5 (optional): age of the largest combined 1966 reduction")

#### 8. Model Fit Statistics (Japanese Fertility Fit-Statistics Table) ####

# (8a) SATURATED POISSON MODEL (ANALYTIC)
LL.sat <- poisson_saturated_loglik(df$y)
EDF.sat <- n

# (8b) DEVIANCES
null.dev <- glm(y ~ 1 + offset(log(pop)), family = poisson, data = df)$deviance
stopifnot(isTRUE(all.equal(null.dev, m.lcsc.ML$null.deviance, tolerance = 1e-6)))
# fitted births of the varying-LC models = GAM fitted births x exp(g(p))
stopifnot(identical(df$y, data$y))
dev.lcsc.ML <- safe_deviance(df$y, fitted(m.lcsc.ML))
dev.vt      <- safe_deviance(df$y, pl.vt.ML$mu)
dev.vt.fh   <- safe_deviance(df$y, pl.vt.fh.ML$mu)

# (8c) LOG-LIKELIHOODS AND EDFs
LL.lcsc.ML <- as.numeric(logLik(m.lcsc.ML)); EDF.lcsc.ML <- sum(m.lcsc.ML$edf)
LL.vt      <- sum(dpois(df$y, pl.vt.ML$mu, log = TRUE));    EDF.vt     <- sum(m.vt.ML$edf) + EDF_PERIOD
LL.vt.fh   <- sum(dpois(df$y, pl.vt.fh.ML$mu, log = TRUE)); EDF.vt.fh  <- sum(m.vt.fh.ML$edf) + EDF_PERIOD

# (8d) OVERDISPERSION ESTIMATE c_hat (Richards 2008, Equation 7)
c_hat <- richards_chat(LL_model = LL.vt.fh, LL_sat = LL.sat, edf_model = EDF.vt.fh, edf_sat = EDF.sat)
disp.lcsc.ML <- pearson_dispersion(m.lcsc.ML)
disp.vt      <- pearson_dispersion(df$y, mu = pl.vt.ML$mu,    edf = EDF.vt)
disp.vt.fh   <- pearson_dispersion(df$y, mu = pl.vt.fh.ML$mu, edf = EDF.vt.fh)
cat("\nOverdispersion:\n")
cat("  variance-to-mean ratio of y:", round(var(df$y) / mean(df$y), 2), "\n")
cat("  c_hat (Richards 2008, m.vt.fh.ML vs saturated):", round(c_hat, 4), "\n")
cat("  Pearson dispersion, m.lcsc.ML:", round(disp.lcsc.ML, 4), "| m.vt.ML:", round(disp.vt, 4),
    "| m.vt.fh.ML:", round(disp.vt.fh, 4), "\n")

# (8e) FIT-STATISTICS TABLE
ic.sat <- c(QAIC=NA_real_, QBIC=NA_real_)  # Saturated likelihood is a benchmark
ic.lcsc.ML <- qaic_qbic(LL.lcsc.ML, EDF.lcsc.ML, c_hat, n)
ic.vt      <- qaic_qbic(LL.vt,      EDF.vt,      c_hat, n)
ic.vt.fh   <- qaic_qbic(LL.vt.fh,   EDF.vt.fh,   c_hat, n)

tab <- data.frame(
  Model     = c("Saturated model", "LC-SC model", "LC-SC model with varying LC curves",
                "LC-SC model with varying LC curves and varying period shocks"),
  Criterion = c("--", m.lcsc.ML$method, m.vt.ML$method, m.vt.fh.ML$method),
  LLV       = c(LL.sat, LL.lcsc.ML, LL.vt, LL.vt.fh),
  R2_D      = c(1, deviance_r2(dev.lcsc.ML, null.dev), deviance_r2(dev.vt, null.dev), deviance_r2(dev.vt.fh, null.dev)),
  Adj_R2_D  = c(NA, adj_deviance_r2(dev.lcsc.ML, null.dev, n = n, edf = EDF.lcsc.ML),
                adj_deviance_r2(dev.vt, null.dev, n = n, edf = EDF.vt),
                adj_deviance_r2(dev.vt.fh, null.dev, n = n, edf = EDF.vt.fh)),
  QAIC      = c(ic.sat["QAIC"], ic.lcsc.ML["QAIC"], ic.vt["QAIC"], ic.vt.fh["QAIC"]),
  QBIC      = c(ic.sat["QBIC"], ic.lcsc.ML["QBIC"], ic.vt["QBIC"], ic.vt.fh["QBIC"]),
  EDF       = c(EDF.sat, EDF.lcsc.ML, EDF.vt, EDF.vt.fh),
  c_hat     = c_hat,
  Pearson_dispersion = c(NA, disp.lcsc.ML, disp.vt, disp.vt.fh),
  deviance  = c(0, dev.lcsc.ML, dev.vt, dev.vt.fh),
  null_deviance = null.dev,
  stringsAsFactors = FALSE)

tab.print <- tab[, c("Model", "Criterion", "LLV", "R2_D", "Adj_R2_D", "QAIC", "QBIC", "EDF", "Pearson_dispersion")]
tab.print$LLV  <- formatC(tab$LLV, format = "f", digits = 2, big.mark = ",")
tab.print$R2_D <- formatC(tab$R2_D, format = "f", digits = 4)
tab.print$Adj_R2_D <- ifelse(is.na(tab$Adj_R2_D), "---", formatC(tab$Adj_R2_D, format = "f", digits = 4))
tab.print$QAIC <- formatC(tab$QAIC, format = "f", digits = 2, big.mark = ",")
tab.print$QBIC <- formatC(tab$QBIC, format = "f", digits = 2, big.mark = ",")
tab.print$EDF  <- formatC(tab$EDF, format = "f", digits = 2, big.mark = ",")
tab.print$Pearson_dispersion <- ifelse(is.na(tab$Pearson_dispersion), "---", formatC(tab$Pearson_dispersion, format = "f", digits = 4))
cat("\nJapanese Fertility: Fit Statistics of LC-SC and Saturated Models\n")
print(tab.print, row.names = FALSE)
cat("c_hat =", round(c_hat, 4), "\n")
write.csv(tab, file.path("Output", "fertility_fit_table.csv"), row.names = FALSE)

# LaTeX rows of the table exactly as typeset in the manuscript
fert_rows <- function(t) {
  r2 <- function(x) if (is.na(x)) "\\multicolumn{1}{c}{---}" else sprintf("$%s$", fmt_tex(x, 4))
  cell <- function(i) sprintf("$%s$ & $%s$ & %s & $%s$ & $%s$ & $%s$",
                              fmt_tex(t$LLV[i], 2), fmt_tex(t$R2_D[i], 4), r2(t$Adj_R2_D[i]),
                              if(is.na(t$QAIC[i])) "---" else fmt_tex(t$QAIC[i], 2), if(is.na(t$QBIC[i])) "---" else fmt_tex(t$QBIC[i], 2), fmt_tex(t$EDF[i], 2))
  c(paste0("Saturated model & ", cell(1), " \\\\[1.5ex]"),
    paste0("LC-SC model & ", cell(2), " \\\\[1.5ex]"),
    paste0("LC-SC model with varying LC curves & ", cell(3), " \\\\[1.5ex]"),
    paste0("LC-SC model with varying LC curves \\\\ and varying period shocks & ", cell(4), " \\\\[2.5ex]"))
}
write_tex(fert_rows(tab), file.path("Output", "tex", "fertility_fit_table.tex"))
cat("LaTeX rows of the fit-statistics table written to Output/tex/fertility_fit_table.tex\n")

# numbers of this example quoted in the manuscript text
mn_add("fert_n_cells", n, 0, "S x T = 3,256 (everywhere)")
mn_add("fert_c_hat", c_hat, 2, "footnote on QAIC/QBIC and table note: estimated nu")
mn_add("fert_n_dummies", length(firehorse), 0, "Section V.5: 18 cohort-specific dummies, ages 18-35")
mn_add("fert_edf_lcsc", EDF.lcsc.ML, 2, "table: EDF of the basic LC-SC model")
mn_add("fert_edf_vt", EDF.vt, 2, "table: EDF of the model with varying LC curves")
mn_add("fert_edf_vt_fh", EDF.vt.fh, 2, "table: EDF of the model with varying LC curves and shocks")
mn_add("fert_edf_lcsc_round", round(EDF.lcsc.ML), 0, "Section V.5: EDF of the basic LC-SC model, rounded")
mn_add("fert_edf_vt_round", round(EDF.vt), 0, "Section V.5: EDF of the varying-LC model, rounded")
mn_add("fert_edf_vt_fh_round", round(EDF.vt.fh), 0, "Section V.5: EDF of the preferred model, rounded")
# EDF breakdown of the preferred model
mn_add("fert_edf_f_vt_fh", smooth_edf(m.vt.fh.ML, smooth_i(m.vt.fh.ML, "s(a)")) + 1, 1,
       "Section V.5: EDF of the overall LC curve LC(a) (theta_1 + f~(a)), preferred model, ML")
mn_add("fert_edf_g_vt_fh", EDF_PERIOD, 0,
       "Section V.5: EDF counted for the common period fluctuations (the yearly leftovers)")
mn_add("fert_n_periods", length(p.index), 0, "Section V.5 / Appendix C.3: number of years, i.e. of period leftovers")
mn_add("fert_edf_h_vt_fh", smooth_edf(m.vt.fh.ML, fs_i(m.vt.fh.ML)), 0,
       "Section V.5: EDF of the cohort-specific deviations, preferred model, ML")
mn_add("fert_edf_fh_vt_fh", sum(m.vt.fh.ML$edf[match(firehorse, names(coef(m.vt.fh.ML)))]), 1,
       "Section V.5 (optional): EDF of the 18 fire-horse dummies, preferred model, ML")
mn_add("fert_n_coef_vt", length(coef(m.vt.ML)), 0, "Section V.5: number of coefficients of the model with varying LC curves")
mn_add("fert_n_coef_vt_fh", length(coef(m.vt.fh.ML)), 0, "Section V.5: number of coefficients of the preferred model")
mn_add("fert_theta1_lcsc", coef(m.lcsc.ML)[["a"]], 3, "Section V.5 (optional): LC slope theta_1, basic LC-SC model, ML")
mn_add("fert_theta1_vt_fh", coef(m.vt.fh.ML)[["a"]], 3, "Section V.5 (optional): LC slope theta_1, preferred model, ML")
mn_add("fert_theta2_lcsc", coef(m.lcsc.ML)[["c"]], 3, "Section V.5 (optional): SC slope theta_2, basic LC-SC model, ML")
mn_add("fert_pearson_vt_fh", disp.vt.fh, 2, "Section V.5 footnote: Pearson dispersion of the preferred ML fit")
mn_add("fert_dQAIC_fh", ic.vt["QAIC"] - ic.vt.fh["QAIC"], 0, "Section V.5 (optional): QAIC reduction from adding the 1966 indicators")
mn_add("fert_dQBIC_fh", ic.vt["QBIC"] - ic.vt.fh["QBIC"], 0, "Section V.5 (optional): QBIC reduction from adding the 1966 indicators")
mn_add("fert_dEDF_fh", EDF.vt.fh - EDF.vt, 0, "Section V.5: EDF added by the 1966 indicators, rounded")
# EDF of each smooth of the basic LC-SC model
for (i in seq_along(m.lcsc.ML$smooth)) {
  s <- m.lcsc.ML$smooth[[i]]; lab <- gsub("[()]", "", s$label)
  mn_add(paste0("fert_edf_", lab, "_lcsc"), smooth_edf(m.lcsc.ML, i), 1,
         paste("Section III.3 footnote: EDF of", s$label, "in the basic LC-SC model (ML)"))
  mn_add(paste0("fert_edfmax_", lab), s$last.para - s$first.para + 1, 0,
         paste("Section III.3 footnote: maximum EDF of", s$label))
}
mn_add("fert_criterion_vt", m.vt.ML$method, where = "table note: criterion of the varying-LC rows")
mn_add("fert_criterion_fig", m.vt.fh$method, where = "figure notes: criterion of the figure models")

## (8f) CONCURVITY OF THE THREE FERTILITY MODELS (footnote in Section V.5)

conc.models <- list("LC-SC (ML)"                                   = rs.lcsc.ML,
                    "LC-SC + varying LC curves (ML)"               = rs.vt.ML,
                    "LC-SC + varying LC curves + 1966 shocks (ML)" = rs.vt.fh.ML)
conc.rows <- list(); conc.full <- list()

for (cm in names(conc.models)) {
  rs.c <- conc.models[[cm]]; b.c <- rs.c$mod; X.c <- rs.c$X; beta.c <- coef(b.c)
  bl <- lcsc_blocks(b.c, X.c, extra = list(lc_slope = "^a$", sc_slope = "^c$", firehorse = "_p\\.index1966$"))
  cat("\n-----------------------------------------------------------------------\n")
  cat("Concurvity:", cm, "\n")
  cat("blocks:", paste(sprintf("%s[%d]", names(bl), sapply(bl, length)), collapse = "  "), "\n")
  for (metric in unique(c("unweighted", CONCURVITY_METRICS))) {   # "unweighted" always: it supplies the quoted values
    w.c <- if (metric == "weighted") b.c$weights else NULL
    cf <- concurvity_blocks(X.c, beta.c, bl, w = w.c, full = TRUE)
    cp <- concurvity_blocks(X.c, beta.c, bl, w = w.c, full = FALSE)
    conc.full[[paste(cm, metric)]] <- cf
    cat("\n  metric:", metric,
        "| rank/ncol:", paste(sprintf("%s %d/%d", names(bl), attr(cf, "rank"), attr(cf, "ncol")), collapse = "  "), "\n")
    cat("  multivariate (term vs the whole rest of the model):\n"); pr(cf, 4)
    cat("  pairwise, worst:\n");    pr(cp$worst, 4)
    cat("  pairwise, observed:\n"); pr(cp$observed, 4)
    for (tg in colnames(cf)) for (ms in rownames(cf))
      conc.rows[[length(conc.rows) + 1]] <- data.frame(
        model = cm, metric = metric, type = "full", measure = ms, term = tg,
        term_given = NA_character_, value = cf[ms, tg],
        ncol = attr(cf, "ncol")[[tg]], rank = attr(cf, "rank")[[tg]], stringsAsFactors = FALSE)
    for (ms in names(cp)) for (tg in rownames(cp[[ms]])) for (gv in colnames(cp[[ms]]))
      if (tg != gv) conc.rows[[length(conc.rows) + 1]] <- data.frame(
        model = cm, metric = metric, type = "pairwise", measure = ms, term = tg,
        term_given = gv, value = cp[[ms]][tg, gv],
        ncol = attr(cp, "ncol")[[tg]], rank = attr(cp, "rank")[[tg]], stringsAsFactors = FALSE)
  }
}
conc.df <- do.call(rbind, conc.rows)
write.csv(conc.df, file.path("Output", "fertility_concurvity.csv"), row.names = FALSE)
cat("\nConcurvity written to Output/fertility_concurvity.csv (", nrow(conc.df), "rows )\n")

# what mgcv's own algorithm returns on the rank-deficient design (not a result)
cat("\nUnguarded computation on the rank-deficient design (mgcv's algorithm, the preferred\n",
    "model): a 'worst' value above 1 for the cohort-specific block is a numerical artifact\n",
    "of the singular solve, NOT concurvity. This is why concurvity_blocks() is used above.\n")
conc.naive <- concurvity_lcsc(rs.vt.fh.ML$mod, rs.vt.fh.ML$X, full = TRUE)
pr(conc.naive, 4)
mn_add("fert_conc_naive_max", max(conc.naive["worst", ]), 3,
       "Appendix F: largest 'worst' value returned by the unguarded computation")

# --- numbers quoted in the manuscript (Appendix F table and text) --------------------- #
key <- function(tm) switch(tm, "s(a)" = "sa", "s(p)" = "sp", "s(c)" = "sc",
                           "s(a,c.cat)" = "h", "firehorse" = "fh", "intercept" = "int",
                           "lc_slope" = "lcs", "sc_slope" = "scs", "linear" = "lin", tm)
cf.lcsc <- conc.full[["LC-SC (ML) unweighted"]]
cf.vtfh <- conc.full[["LC-SC + varying LC curves + 1966 shocks (ML) unweighted"]]
cf.vt   <- conc.full[["LC-SC + varying LC curves (ML) unweighted"]]
for (msr in c("worst", "observed", "estimate")) {
  for (tm in c("s(a)", "s(p)", "s(c)"))
    mn_add(paste0("fert_conc_lcsc_", msr, "_", key(tm)), cf.lcsc[msr, tm], 3,
           "Appendix F: multivariate concurvity, basic LC-SC model (ML)")
  # the varying-LC models have no SC curve and no period term
  for (tm in c("intercept", "lc_slope", "s(a)", "s(a,c.cat)", "firehorse"))
    mn_add(paste0("fert_conc_vtfh_", msr, "_", key(tm)), cf.vtfh[msr, tm], 3,
           "Appendix F: multivariate concurvity, preferred model (ML)")
  for (tm in c("s(a)", "s(a,c.cat)"))
    mn_add(paste0("fert_conc_vt_", msr, "_", key(tm)), cf.vt[msr, tm], 3,
           "Appendix F (optional): multivariate concurvity, varying-LC model (ML)")
}
# the values in the metric of the final IRLS working weights
if ("weighted" %in% CONCURVITY_METRICS) {
  cw <- conc.full[["LC-SC + varying LC curves + 1966 shocks (ML) weighted"]]
  for (tm in c("s(a)", "s(a,c.cat)"))
    mn_add(paste0("fert_conc_w_vtfh_observed_", key(tm)), cw["observed", tm], 3,
           "Appendix F (optional): observed multivariate concurvity, final preferred-ML PIRLS weights")
  rm(cw)
}
# rank of the cohort-specific block
ob.fs <- rs.vt.fh.ML$ortho[[fs_i(m.vt.fh.ML)]]
H.raw <- mgcv::PredictMat(ob.fs$sm, data)
d.all <- svd(H.raw, nu = 0, nv = 0)$d
thr <- max(d.all) * 1e-8                                   # the tolerance of concurvity_blocks()
rank_abs <- function(M) { d <- svd(M, nu = 0, nv = 0)$d; sum(d > thr) }
per.coh <- do.call(rbind, lapply(sort(unique(data$c.index)), function(cc) {
  i <- data$c.index == cc
  cols <- which(colSums(abs(H.raw[i, , drop = FALSE])) > 0)
  data.frame(cohort = cc, n_ages = sum(i), rank = rank_abs(H.raw[i, cols, drop = FALSE]))
}))
rk.unproj <- sum(d.all > thr); rk.used <- attr(cf.vtfh, "rank")[["s(a,c.cat)"]]
if (rk.unproj != sum(per.coh$rank))  # the block is block-diagonal by cohort
  warning("rank of the block (", rk.unproj, ") differs from the sum of the cohort ranks (", sum(per.coh$rank), ")")
if (rk.unproj != rk.used)             # not projected: the block used in the fit is the raw block
  warning("rank of the block used in the fit (", rk.used, ") differs from that of the raw block (", rk.unproj, ")")
loss.few  <- sum(K_COHORT - per.coh$rank[per.coh$n_ages < K_COHORT])
loss.part <- sum(K_COHORT - per.coh$rank[per.coh$n_ages >= K_COHORT])
# rank of the whole design of the preferred model
d.X <- svd(rs.vt.fh.ML$X, nu = 0, nv = 0)$d
rk.X <- sum(d.X > max(d.X) * 1e-8)
cat(sprintf("\nCohort-specific block: %d columns; rank %d over the observed cells (not projected)\n",
            attr(cf.vtfh, "ncol")[["s(a,c.cat)"]], rk.unproj))
cat(sprintf("  (%d dimensions lost in the %d cohorts observed at fewer than %d ages; %d in the %d cohorts observed\n",
            loss.few, sum(per.coh$n_ages < K_COHORT), K_COHORT, loss.part, sum(per.coh$n_ages >= K_COHORT & per.coh$rank < K_COHORT)))
cat(sprintf("   at %d or more ages but over part of the age range)\n", K_COHORT))
cat(sprintf("Whole design of the preferred model: %d columns, rank %d (sum of the block ranks: %d)\n",
            ncol(rs.vt.fh.ML$X), rk.X, sum(attr(cf.vtfh, "rank"))))
mn_add("fert_conc_fs_ncol", attr(cf.vtfh, "ncol")[["s(a,c.cat)"]], 0, "Section V.5 / Appendix F: columns of the cohort-specific block")
mn_add("fert_fs_rank_unproj", rk.unproj, 0, "Appendix F: rank of the cohort-specific block over the observed cells")
mn_add("fert_design_ncol", ncol(rs.vt.fh.ML$X), 0, "Appendix F: columns of the whole design of the preferred model")
mn_add("fert_design_rank", rk.X, 0, "Appendix F: rank of the whole design of the preferred model")
mn_add("fert_design_rank_overlap", sum(attr(cf.vtfh, "rank")) - rk.X, 0,
       "Appendix F (optional): dimensions shared by the blocks of the preferred model (sum of block ranks minus the rank of the design)")
mn_add("fert_fs_loss_fewages", loss.few, 0, "Appendix F (optional): dimensions lost in cohorts observed at fewer than 15 ages")
mn_add("fert_fs_loss_partial", loss.part, 0, "Appendix F (optional): dimensions lost in cohorts observed over part of the age range")
rm(H.raw, per.coh, ob.fs, d.all, d.X)
rm(conc.models, conc.rows, rs.c, b.c, X.c, beta.c); invisible(gc())

rm(rs.lcsc.ML, m.lcsc.ML, rs.vt.ML, rs.vt.fh.ML, H.vt.ML, H.fh.ML); invisible(gc())

#### 9. Figures From the Preferred Model (m.vt.fh, REML) ####

# term contributions of m.vt.fh on the observed cells (Sections 6c and 7)
b0.fh <- H.vt.fh$intercept
lc.a <- as.numeric(H.vt.fh$est[, "a"] + H.vt.fh$est[, "s(a)"])
d.ac <- as.numeric(H.vt.fh$est[, fs.label])
g.p  <- as.numeric(pl.vt.fh$g_cell)
phi.add <- rep(0, nrow(data))
phi.add[is.fh] <- add.fh[paste0("a.index", data$a.index[is.fh], "_p.index1966")]
lp0  <- b0.fh + lc.a + d.ac  # cohort-specific LC curves
pred <- data.frame(a = data$a.index, c = data$c.index, p = data$p.index,
                   lc.a = lc.a, d.ac = d.ac, g.p = g.p, phi = phi.add, lp0 = lp0, asfr0 = exp(lp0) * 1000)
stopifnot(isTRUE(all.equal(exp(lp0 + g.p + phi.add + log(data$pop)), mu.fh, tolerance = 1e-6)))
# peak ages of the cohort-specific LC curves
pk.coh <- 1930:1985
peak_age <- function(v) sapply(pk.coh, function(u) { i <- pred$c == u; pred$a[i][which.max(v[i])] })
pk.net <- peak_age(pred$lp0); pk.fit <- peak_age(pred$lp0 + pred$g.p + pred$phi)
cat(sprintf("Peak ages, cohorts 1930-1985: cohort-specific LC curves %d-%d, fitted careers %d-%d; 1950: %d/%d, 1980: %d/%d\n",
            min(pk.net), max(pk.net), min(pk.fit), max(pk.fit), pk.net[pk.coh == 1950], pk.fit[pk.coh == 1950],
            pk.net[pk.coh == 1980], pk.fit[pk.coh == 1980]))
mn_add("fert_peak_net_1950", pk.net[pk.coh == 1950], 0, "Section V.5: age at which the 1950 cohort's LC curve peaks")
mn_add("fert_peak_net_1980", pk.net[pk.coh == 1980], 0, "Section V.5: age at which the 1980 cohort's LC curve peaks")
mn_add("fert_peak_net_1930", pk.net[pk.coh == 1930], 0, "Section V.5 (optional): age at which the 1930 cohort's LC curve peaks")
mn_add("fert_peak_net_lo", min(pk.net), 0, "Section V.5 (optional): earliest peak age of the cohort LC curves, 1930-1985")
mn_add("fert_peak_net_hi", max(pk.net), 0, "Section V.5 (optional): latest peak age of the cohort LC curves, 1930-1985")

# how the fitted surface is divided between the overall LC curve and the cohort curves
w.fh <- as.numeric(fitted(m.vt.fh))
dev.coh <- tapply(pred$d.ac * w.fh, pred$c, sum) / tapply(w.fh, pred$c, sum)
dev.age <- tapply(pred$d.ac * w.fh, pred$a, sum) / tapply(w.fh, pred$a, sum)
cat(sprintf("Births-weighted mean deviation: within cohorts %.2f to %.2f (1950: %.2f, 1980: %.2f); within ages %.2f to %.2f\n",
            min(dev.coh), max(dev.coh), dev.coh[["1950"]], dev.coh[["1980"]], min(dev.age), max(dev.age)))
mn_add("fert_devmean_age_maxabs", max(abs(dev.age)), 2,
       "Appendix C.3 (optional): largest births-weighted mean of the deviations within an age (linear predictor)")

# reference ASFR for the period figure
ref.row <- which(pred$a == 27 & pred$c == 1950)
stopifnot(length(ref.row) == 1)
asfr.ref <- exp(pred$lp0[ref.row]) * 1000
cat(sprintf("\nReference ASFR (age 27, cohort 1950, period terms set to zero): %.2f per 1,000\n", asfr.ref))
cat(sprintf("  intercept = %.4f, LC(27) = %.4f, delta_1950(27) = %.4f\n",
            b0.fh, pred$lc.a[ref.row], pred$d.ac[ref.row]))

# colors: one per cohort on a spectral scale, violet/blue = oldest, red = newest
sorted_c <- sort(unique(pred$c), decreasing = TRUE)
original_palette <- brewer.pal(11, "Spectral")
no_yellow_palette <- original_palette[-c(5:6)] # exclude yellowish colors
color_palette <- colorRampPalette(no_yellow_palette)(length(sorted_c))
named_palette <- setNames(adjustcolor(color_palette, alpha.f = 0.9), sorted_c)
selected_c <- c(1892, 1920, 1950, 1980, 2008)

## (9a) OVERALL LC CURVE AND COHORT-SPECIFIC DEVIATIONS ON THE LINEAR-PREDICTOR SCALE
pdf(file.path(fig.dir, "Vertical_TwoPanel.pdf"), width = 6.8, height = 13.5) # dimensions in inches (nearly square panels)
par(mfrow = c(2, 1))

# panel (a): overall LC curve LC(a)
newdata <- unique(pred[, c("a", "lc.a")]); newdata <- newdata[order(newdata$a), ]
stopifnot(nrow(newdata) == length(a.index)) # one value per age
ylim.a <- ylim_or_default(newdata$lc.a, c(-12, 7))
plot(x = newdata$a, y = newdata$lc.a, type = "l", xlim = c(10, 57),
     ylim = ylim.a, xlab = "Age", ylab = "Linear Predictor", lwd = 1.2, cex.lab = 1.1,
     yaxt = "n", xaxt = "n", main = "(a) Overall LC Curve", cex.main = 1.5)
axis(side = 2, at = seq(from = -40, to = 40, by = 2), las = 1, cex.axis = 0.9)
axis(side = 1, at = seq(from = 15, to = 55, by = 5), las = 0, cex.axis = 0.9)

# panel (b): cohort-specific deviations delta_u(a), each including its cohort's level
newdata <- pred[, c("a", "c", "d.ac")]
ylim.d <- ylim_or_default(newdata$d.ac, c(-3, 3))
plot(x = newdata$a, y = newdata$d.ac, type = "n", xlim = c(10, 57),
     ylim = ylim.d, xlab = "Age", ylab = "Linear Predictor", cex.lab = 1.1,
     yaxt = "n", xaxt = "n", main = "(b) Cohort-Specific Deviations", cex.main = 1.5)
axis(side = 2, at = seq(from = -40, to = 40, by = 1), las = 1, cex.axis = 0.9)
axis(side = 1, at = seq(from = 15, to = 55, by = 5), las = 0, cex.axis = 0.9)
abline(h = 0, lty = "dashed", col = "darkgray")   # zero deviation = the overall LC curve
for (cc in sorted_c) {
  sub <- newdata[newdata$c == cc, ]; sub <- sub[order(sub$a), ]
  lines(sub$a, sub$d.ac, col = named_palette[as.character(cc)], lwd = 1)
}
legend("topleft", inset = c(0.02, 0.01), legend = selected_c, col = named_palette[as.character(selected_c)],
       lwd = 1.45, title = "Cohort", bty = "n")   # the upper left of the panel is empty
invisible(dev.off())
cat(sprintf("Vertical_TwoPanel: LC(a) range %.2f to %.2f; delta_u(a) range %.2f to %.2f\n",
            min(pred$lc.a), max(pred$lc.a), min(pred$d.ac), max(pred$d.ac)))
mn_add("fert_d_range_lo", min(pred$d.ac), 2, "Appendix figure (optional): smallest cohort-specific deviation (linear predictor)")
mn_add("fert_d_range_hi", max(pred$d.ac), 2, "Appendix figure (optional): largest cohort-specific deviation (linear predictor)")

# (9a') COHORT-SPECIFIC LC CURVES ON THE ASFR SCALE (APPENDIX E)
pred$asfr.common <- exp(b0.fh + pred$lc.a) * 1000                    # overall LC curve as an ASFR
pred$asfr.dev    <- pred$asfr0 - pred$asfr.common                    # contribution of delta_u(a) in ASFR units
common.curve <- unique(pred[, c("a", "asfr.common")]); common.curve <- common.curve[order(common.curve$a), ]
pdf(file.path(fig.dir, "Vertical_TwoPanel_ASFR.pdf"), width = 8.6, height = 13.5) # dimensions in inches
par(mfrow = c(2, 1))
# panel (a): cohort-specific LC curves (ASFR per 1,000) and the overall LC curve
ylim.a2 <- ylim_or_default(c(pred$asfr0, common.curve$asfr.common), c(0, 300))
plot(x = pred$a, y = pred$asfr0, type = "n", xlim = c(10, 70), ylim = ylim.a2,
     xlab = "Age", ylab = "Age-Specific Fertility Rate per 1,000 Women",
     yaxt = "n", xaxt = "n", main = "(a) Cohort-Specific LC Curves", cex.main = 1.5)
axis(side = 2, at = seq(from = 0, to = 2000, by = 50), las = 1, cex.axis = 0.8)
axis(side = 1, at = seq(from = 15, to = 55, by = 5), las = 0, cex.axis = 0.8)
for (cc in sorted_c) {
  sub <- pred[pred$c == cc, ]; sub <- sub[order(sub$a), ]
  lines(sub$a, sub$asfr0, col = named_palette[as.character(cc)], lwd = 1)
}
lines(common.curve$a, common.curve$asfr.common, col = "black", lwd = 2.5)
legend(x = 58, y = ylim.a2[2] - 0.05 * diff(ylim.a2), legend = c(selected_c, "Overall"),
       col = c(named_palette[as.character(selected_c)], "black"), lwd = c(rep(1.45, length(selected_c)), 2.5),
       title = "Cohort", bty = "n")
# panel (b): contribution of the deviations in births per 1,000 women
ylim.b2 <- ylim_or_default(pred$asfr.dev, c(-100, 100))
plot(x = pred$a, y = pred$asfr.dev, type = "n", xlim = c(10, 70), ylim = ylim.b2,
     xlab = "Age", ylab = "Births per 1,000 Women",
     yaxt = "n", xaxt = "n", main = "(b) Deviations from the Overall LC Curve", cex.main = 1.5)
axis(side = 2, at = seq(from = -500, to = 500, by = 25), las = 1, cex.axis = 0.8)
axis(side = 1, at = seq(from = 15, to = 55, by = 5), las = 0, cex.axis = 0.8)
abline(h = 0, lty = "dashed", col = "darkgray")
for (cc in sorted_c) {
  sub <- pred[pred$c == cc, ]; sub <- sub[order(sub$a), ]
  lines(sub$a, sub$asfr.dev, col = named_palette[as.character(cc)], lwd = 1)
}
legend(x = 58, y = ylim.b2[2] - 0.05 * diff(ylim.b2), legend = selected_c, col = named_palette[as.character(selected_c)],
       lwd = 1.45, title = "Cohort", bty = "n")
invisible(dev.off())
cat(sprintf("Vertical_TwoPanel_ASFR: cohort-specific curves peak between %.0f and %.0f per 1,000; deviation contributions range %.0f to %.0f per 1,000\n",
            min(tapply(pred$asfr0, pred$c, max)), max(tapply(pred$asfr0, pred$c, max)), min(pred$asfr.dev), max(pred$asfr.dev)))
mn_add("fert_asfr_dev_min", min(pred$asfr.dev), 0, "Appendix figure (optional): largest negative contribution of the deviations, per 1,000")
mn_add("fert_asfr_dev_max", max(pred$asfr.dev), 0, "Appendix figure (optional): largest positive contribution of the deviations, per 1,000")

## (9a'') HOW THE PIECES ADD UP: BUILD-UP OF THE CAREERS OF SELECTED COHORTS (APPENDIX E)
BUILDUP_COHORTS <- c(1925, 1940, 1980)
phi.cell <- pred$phi                                                  # additional 1966 shock in each cell (0 elsewhere)
pred$asfr.period <- exp(pred$lp0 + pred$g.p) * 1000                  # + common period fluctuation
pred$asfr.fit    <- mu.fh / data$pop * 1000                          # + 1966 shock = fitted career
stopifnot(isTRUE(all.equal(pred$asfr.fit, exp(pred$lp0 + pred$g.p + phi.cell) * 1000, tolerance = 1e-6)))
pred$asfr.obs <- data$y / data$pop * 1000
build <- lapply(BUILDUP_COHORTS, function(cc) { s <- pred[pred$c == cc, ]; s[order(s$a), ] })
names(build) <- BUILDUP_COHORTS

# contributions in births per 1,000 women, in the order of the terms of the equation
for (i in seq_along(build)) {
  s <- build[[i]]
  s$D1 <- s$asfr0 - s$asfr.common; s$D2 <- s$asfr.period - s$asfr0; s$D3 <- s$asfr.fit - s$asfr.period
  stopifnot(isTRUE(all.equal(s$asfr.common + s$D1 + s$D2 + s$D3, s$asfr.fit)))
  s$lp.common <- b0.fh + s$lc.a                                   # log overall LC curve
  s$lp.fit <- s$lp.common + s$d.ac + s$g.p + s$phi                # log fitted career
  stopifnot(isTRUE(all.equal(exp(s$lp.fit) * 1000, s$asfr.fit)),
            isTRUE(all.equal(s$D1, s$asfr.common * (exp(s$d.ac) - 1))),
            isTRUE(all.equal(s$D2, s$asfr.common * exp(s$d.ac) * (exp(s$g.p) - 1))),
            isTRUE(all.equal(s$D3, s$asfr.common * exp(s$d.ac + s$g.p) * (exp(s$phi) - 1))))
  build[[i]] <- s
}
ylim.bu <- ylim_or_default(unlist(lapply(build, function(s) c(s$asfr.fit, s$asfr.common, s$D1, s$D2, s$D3))), c(-150, 300))
col.bu <- c(common = "gray55", cohort = brewer.pal(8, "Set1")[2],
            period = brewer.pal(8, "Set1")[3], shock = brewer.pal(8, "Set1")[1], fitted = "black")
lwd.bu <- c(2, 1.6, 1.6, 1.6, 2.4)

# the contributions as curves
pdf(file.path(fig.dir, "CohortBuildUp_Lines.pdf"), width = 11, height = 8.5) # dimensions in inches
par(mfrow = c(2, 2))
for (i in seq_along(build)) {
  s <- build[[i]]
  plot(x = s$a, y = s$asfr.fit, type = "n", xlim = c(10, 59), ylim = ylim.bu,
       xlab = "Age", ylab = if (i != 2) "Births per 1,000 Women" else " ", cex.lab = 1.1,
       yaxt = "n", xaxt = "n", main = paste0("(", letters[i], ") Cohort ", names(build)[i]), cex.main = 1.4)
  axis(side = 2, at = seq(from = -500, to = 2000, by = 50), las = 1, cex.axis = 0.9)
  axis(side = 1, at = seq(from = 15, to = 55, by = 5), las = 0, cex.axis = 0.9)
  abline(h = 0, lty = "dashed", col = "darkgray")
  lines(s$a, s$asfr.common, col = col.bu["common"], lwd = 2)
  lines(s$a, s$D1, col = col.bu["cohort"], lwd = 1.6)
  lines(s$a, s$D2, col = col.bu["period"], lwd = 1.6)
  if (any(s$D3 != 0)) lines(s$a, s$D3, col = col.bu["shock"], lwd = 1.6)
  lines(s$a, s$asfr.fit, col = col.bu["fitted"], lwd = 2.4)
}
plot.new()                                          # lower right cell: the legend
legend("center", bty = "n", cex = 1.3, lwd = lwd.bu, col = col.bu,
       legend = c("Overall LC curve", "Cohort-specific deviation", "Common period fluctuation",
                  "1966 shock (1940 cohort only)", "Fitted career = sum of the four"))
invisible(dev.off())

# the same three cohorts on the scale of the linear predictor
ylim.log <- ylim_or_default(unlist(lapply(build, function(s) c(s$lp.common, s$lp.fit, s$d.ac, s$g.p, s$phi))), c(-12, 1))
# multiplicative factors on a logarithmic right axis
ylim.fac <- range(unlist(lapply(build, function(s) exp(c(s$d.ac, s$g.p, s$phi))))) * c(1 / 1.25, 1.25)
fac.ticks <- c(0.03125, 0.0625, 0.125, 0.25, 0.5, 1, 2, 4, 8, 16, 32, 64)
fac.ticks <- fac.ticks[fac.ticks >= ylim.fac[1] & fac.ticks <= ylim.fac[2]]
ylim.bu2 <- ylim_or_default(unlist(lapply(build, function(s) c(s$asfr.fit, s$asfr.common))), c(0, 300))
pdf(file.path(fig.dir, "CohortBuildUp_Log.pdf"), width = 10.5, height = 9) # dimensions in inches
layout(matrix(c(1, 4, 2, 5, 3, 6, 7, 7), ncol = 2, byrow = TRUE), heights = c(1, 1, 1, 0.42))
par(cex = 0.83, mar = c(3.6, 5, 2.6, 5.2), mgp = c(2.4, 0.8, 0))   # x-axis title close to the axis (not clipped)
for (i in seq_along(build)) {                       # left column: linear-predictor scale
  s <- build[[i]]
  plot(x = s$a, y = s$lp.fit, type = "n", xlim = c(10, 59), ylim = ylim.log,
       xlab = if (i == 3) "Age" else "", ylab = "",
       cex.lab = 1.05, yaxt = "n", xaxt = "n", main = paste0("(", letters[i], ") Cohort ", names(build)[i]), cex.main = 1.3)
  if (i == 2) title(ylab = "Linear Predictor (log fertility rate)", line = 3.1, cex.lab = 1.05)
  axis(side = 2, at = seq(from = -40, to = 40, by = 2), las = 1, cex.axis = 0.9, mgp = c(3, 1.1, 0))
  axis(side = 1, at = seq(from = 15, to = 55, by = 5), las = 0, cex.axis = 0.9)
  abline(h = 0, lty = "dashed", col = "darkgray")
  lines(s$a, s$lp.common, col = col.bu["common"], lwd = 2)
  lines(s$a, s$d.ac, col = col.bu["cohort"], lwd = 1.6)
  lines(s$a, s$g.p, col = col.bu["period"], lwd = 1.6)
  if (any(s$phi != 0)) lines(s$a, s$phi, col = col.bu["shock"], lwd = 1.6)
  lines(s$a, s$lp.fit, col = col.bu["fitted"], lwd = 2.4)
}
for (i in seq_along(build)) {                       # right column: exponentiated
  s <- build[[i]]
  plot(x = s$a, y = s$asfr.fit, type = "n", xlim = c(10, 59), ylim = ylim.bu2,
       xlab = if (i == 3) "Age" else "", ylab = "",
       cex.lab = 1.05, yaxt = "n", xaxt = "n", main = paste0("(", letters[i + 3], ") Cohort ", names(build)[i]), cex.main = 1.3)
  if (i == 2) title(ylab = "Births per 1,000 Women", line = 3.1, cex.lab = 1.05)
  axis(side = 2, at = seq(from = 0, to = 2000, by = 50), las = 1, cex.axis = 0.9, mgp = c(3, 1.1, 0))
  axis(side = 1, at = seq(from = 15, to = 55, by = 5), las = 0, cex.axis = 0.9)
  lines(s$a, s$asfr.common, col = col.bu["common"], lwd = 2)
  lines(s$a, s$asfr.fit, col = col.bu["fitted"], lwd = 2.4)
  par(new = TRUE)                                    # right axis: multiplicative factors
  plot(x = s$a, y = exp(s$d.ac), type = "n", xlim = c(10, 59), ylim = ylim.fac, axes = FALSE, xlab = "", ylab = "", log = "y")
  axis(side = 4, at = fac.ticks, labels = as.character(fac.ticks), las = 1, cex.axis = 0.9, mgp = c(3, 1.1, 0))
  if (i == 2) mtext("Multiplicative Factor (log scale)", side = 4, line = 4.1, cex = 0.83 * 1.05)
  abline(h = 1, lty = "dashed", col = "darkgray")
  lines(s$a, exp(s$d.ac), col = col.bu["cohort"], lwd = 1.6)
  lines(s$a, exp(s$g.p), col = col.bu["period"], lwd = 1.6)
  if (any(s$phi != 0)) lines(s$a, exp(s$phi), col = col.bu["shock"], lwd = 1.6)
}
par(mar = c(0, 0, 0, 0))                             # bottom: one legend for all six panels
plot.new()
legend("center", bty = "n", ncol = 3, cex = 1.05, lwd = lwd.bu, col = col.bu,
       legend = c("Overall LC curve", "Cohort-specific deviation", "Common period fluctuation",
                  "1966 shock (1940 cohort only)", "Fitted career"))
invisible(dev.off())

# numbers quoted in the notes of the two build-up figures
bu.rows <- list()
for (i in seq_along(build)) {
  s <- build[[i]]
  cat(sprintf("CohortBuildUp: cohort %s observed at ages %d-%d; max |fitted - observed| = %.1f per 1,000\n",
              names(build)[i], min(s$a), max(s$a), max(abs(s$asfr.fit - s$asfr.obs))))
  for (tm in c("deviation", "period", "shock")) {
    D <- switch(tm, deviation = s$D1, period = s$D2, shock = s$D3)
    fac <- exp(switch(tm, deviation = s$d.ac, period = s$g.p, shock = s$phi))
    bu.rows[[length(bu.rows) + 1]] <- data.frame(
      cohort = as.numeric(names(build)[i]), term = tm,
      max_births = max(D), age_max = s$a[which.max(D)], min_births = min(D), age_min = s$a[which.min(D)],
      sum_positive = sum(pmax(D, 0)), sum_negative = sum(pmin(D, 0)),
      max_factor = max(fac), age_max_factor = s$a[which.max(fac)],
      min_factor = min(fac), age_min_factor = s$a[which.min(fac)],
      ages_positive = if (any(D > 0.5)) paste(range(s$a[D > 0.5]), collapse = "-") else NA_character_,
      stringsAsFactors = FALSE)
  }
  bu.rows[[length(bu.rows)]]$max_abs_fit_minus_obs <- max(abs(s$asfr.fit - s$asfr.obs))
}
bu <- do.call(rbind, lapply(bu.rows, function(r) { if (is.null(r$max_abs_fit_minus_obs)) r$max_abs_fit_minus_obs <- NA; r }))
cat("\nBuild-up summary:\n")
print(format(bu, digits = 4), row.names = FALSE)
# the entries quoted in the notes of the two build-up figures are required
bu.quoted <- c("fert_bu_1925_deviation_max", "fert_bu_1925_period_max",
               "fert_bu_1940_period_min", "fert_bu_1940_shock_min",
               "fert_bu_1980_deviation_min", "fert_bu_1980_deviation_max", "fert_bu_1980_period_min",
               "fert_bu_1925_period_facmin", "fert_bu_1925_period_facmax",
               "fert_bu_1980_deviation_facmin", "fert_bu_1980_deviation_facmax",
               "fert_bu_1980_period_facmin", "fert_bu_1980_period_facmax")
bu_add <- function(key, value, digits, what) {
  mn_add(key, value, digits, paste0(what, if (key %in% bu.quoted) "" else " (optional)"))
}
for (r in seq_len(nrow(bu))) {
  kk <- sprintf("fert_bu_%d_%s_", bu$cohort[r], bu$term[r])
  wh <- sprintf("Appendix figure notes (build-up), cohort %d, %s:", bu$cohort[r], bu$term[r])
  bu_add(paste0(kk, "max"), bu$max_births[r], 0, paste(wh, "largest positive contribution, births per 1,000"))
  bu_add(paste0(kk, "min"), -bu$min_births[r], 0, paste(wh, "largest negative contribution, births per 1,000, sign dropped"))
  bu_add(paste0(kk, "sumneg"), -bu$sum_negative[r], 0, paste(wh, "total births removed per 1,000"))
  bu_add(paste0(kk, "sumpos"), bu$sum_positive[r], 0, paste(wh, "total births added per 1,000"))
  bu_add(paste0(kk, "facmax"), bu$max_factor[r], 2, paste(wh, "largest multiplicative factor"))
  bu_add(paste0(kk, "facmin"), bu$min_factor[r], 2, paste(wh, "smallest multiplicative factor"))
}

# (9b) OVERALL LC CURVE AS AN ASFR (APPENDIX E)
newdata <- unique(pred[, c("a", "lc.a")]); newdata <- newdata[order(newdata$a), ]
newdata$yhat_agemain <- exp(b0.fh + newdata$lc.a) * 1000
cat(sprintf("AgeMain_Fertility: overall LC curve peaks at %.1f per 1,000 at age %d\n",
            max(newdata$yhat_agemain), newdata$a[which.max(newdata$yhat_agemain)]))
mn_add("fert_agemain_peak_age", newdata$a[which.max(newdata$yhat_agemain)], 0, "Appendix figure: age at which LC(a) peaks")
mn_add("fert_agemain_peak_asfr", max(newdata$yhat_agemain), 0, "Appendix figure (optional): peak of LC(a) in births per 1,000")
ylim.age <- ylim_or_default(newdata$yhat_agemain, c(0, 250))
pdf(file.path(fig.dir, "AgeMain_Fertility.pdf"), width = 8.25, height = 7.75) # dimensions in inches
plot(x = newdata$a, y = newdata$yhat_agemain, type = "l", xlim = c(10, 59),
     ylim = ylim.age, xlab = "Age", ylab = "Age-Specific Fertility Rate per 1,000 Women",
     yaxt = "n", xaxt = "n", main = " ", cex.main = 0.9)
axis(side = 2, at = seq(from = 0, to = 2000, by = 25), las = 1, cex.axis = 0.8)
axis(side = 1, at = seq(from = 15, to = 55, by = 5), las = 0, cex.axis = 0.8)
invisible(dev.off())

# (9c) COHORT-SPECIFIC LC CURVES AND AGE-SPECIFIC SC CURVES (APPENDIX E)
newdata <- pred[, c("a", "c", "asfr0")]; names(newdata)[3] <- "yhat_agecurves"
newdata <- newdata[order(newdata$c, newdata$a), ]
pdf(file.path(fig.dir, "TwoPanel_LCSC.pdf"), width = 14, height = 7.25) # dimensions in inches
par(mfrow = c(1, 2))
ylim.p1 <- ylim_or_default(newdata$yhat_agecurves[newdata$c %in% c(1950, 1965, 1980)], c(0, 250))
sub_data <- subset(newdata, c == 1950)
plot(x = sub_data$a, y = sub_data$yhat_agecurves, type = "l", xlim = c(10, 59),
     ylim = ylim.p1, xlab = "Age", ylab = "Age-Specific Fertility Rate per 1,000 Women",
     yaxt = "n", xaxt = "n", main = "(a) Cohort-Specific LC Curves", cex.main = 1.4)
axis(side = 2, at = seq(from = 0, to = 2000, by = 25), las = 1, cex.axis = 0.8)
axis(side = 1, at = seq(from = 15, to = 55, by = 5), las = 0, cex.axis = 0.8)
max_point <- sub_data[which.max(sub_data$yhat_agecurves), ]
text(x = max_point$a, y = max_point$yhat_agecurves, labels = max_point$c, pos = 4, cex = 0.5)
for (cc in c(1980, 1965)) {
  sub_data <- subset(newdata, c == cc)
  lines(x = sub_data$a, y = sub_data$yhat_agecurves)
  max_point <- sub_data[which.max(sub_data$yhat_agecurves), ]
  text(x = max_point$a, y = max_point$yhat_agecurves, labels = max_point$c, pos = 4, cex = 0.5)
}
ylim.p2 <- ylim_or_default(newdata$yhat_agecurves[newdata$a %in% c(25, 30, 35, 40)], c(0, 300))
sub_data <- subset(newdata, a == 25)
plot(x = sub_data$c, y = sub_data$yhat_agecurves, type = "l", xlim = c(1900, 2005),
     ylim = ylim.p2, xlab = "Cohort", ylab = " ",
     yaxt = "n", xaxt = "n", main = "(b) Age-Specific SC Curves", cex.main = 1.4)
axis(side = 2, at = seq(from = 0, to = 5000, by = 50), las = 1, cex.axis = 0.8)
axis(side = 1, at = seq(from = 1890, to = 2000, by = 20), las = 0, cex.axis = 0.8)
last_point <- tail(sub_data, n = 1)
text(x = last_point$c, y = last_point$yhat_agecurves, labels = last_point$a, pos = 4, cex = 0.5)
for (aa in c(30, 35, 40)) {
  sub_data <- subset(newdata, a == aa)
  lines(x = sub_data$c, y = sub_data$yhat_agecurves)
  last_point <- tail(sub_data, n = 1)
  text(x = last_point$c, y = last_point$yhat_agecurves, labels = last_point$a, pos = 4, cex = 0.5)
}
invisible(dev.off())
for (cc in c(1950, 1965, 1980)) {
  i <- newdata$c == cc
  cat(sprintf("TwoPanel_LCSC (a): cohort %d peaks at %.1f per 1,000 at age %d\n",
              cc, max(newdata$yhat_agecurves[i]), newdata$a[i][which.max(newdata$yhat_agecurves[i])]))
}
for (cc in c(1925, 1950, 1980)) for (aa in c(25, 40)) {
  v <- newdata$yhat_agecurves[newdata$c == cc & newdata$a == aa]
  cat(sprintf("TwoPanel_LCSC (b): cohort %d at age %d: %.1f per 1,000\n", cc, aa, v))
  mn_add(sprintf("fert_asfr_c%d_a%d", cc, aa), v, 0,
         paste("Section V.5", if (cc == 1925) "(optional)" else "", ": ASFR per 1,000 of the cohort-specific LC curve (TwoPanel_LCSC, panel b)"))
}
for (aa in c(25, 40)) {
  r <- newdata$yhat_agecurves[newdata$c == 1950 & newdata$a == aa] / newdata$yhat_agecurves[newdata$c == 1980 & newdata$a == aa]
  cat(sprintf("TwoPanel_LCSC (b): ratio 1950/1980 at age %d: %.2f\n", aa, r))
  mn_add(sprintf("fert_ratio_1950_1980_a%d", aa), r, 1, "Section V.5: ratio of the 1950 to the 1980 cohort-specific LC curve at this age")
}

## (9d) COMMON PERIOD FLUCTUATIONS WITH THE AVERAGE 1966 FIRE-HORSE SHOCK (MAIN TEXT)
newdata <- data.frame(p = p.index,
                      yhat_pcurves = asfr.ref * exp(period.tab$vtfh_g_with_shock),
                      lo = asfr.ref * exp(period.tab$combined_lo),
                      hi = asfr.ref * exp(period.tab$combined_hi))
cat(sprintf("PeriodRE_Fertility: reference %.2f, range %.2f-%.2f, 1965 = %.2f, 1966 = %.2f, 1967 = %.2f\n",
            asfr.ref, min(newdata$yhat_pcurves), max(newdata$yhat_pcurves),
            newdata$yhat_pcurves[newdata$p == 1965], newdata$yhat_pcurves[newdata$p == 1966], newdata$yhat_pcurves[newdata$p == 1967]))
cat(sprintf("  widest pointwise interval: %.2f per 1,000 (period %d)\n",
            max(newdata$hi - newdata$lo), newdata$p[which.max(newdata$hi - newdata$lo)]))
# logarithmic vertical axis
ylim.per <- range(c(newdata$lo, newdata$hi, asfr.ref)) * c(0.9, 1.1)
# four labeled ticks within the plotted range (about 134-225 per 1,000)
ticks.per <- c(140, 160, 180, 200)
pdf(file.path(fig.dir, "PeriodRE_Fertility.pdf"), width = 8.25, height = 7.75) # dimensions in inches
plot(x = newdata$p, y = newdata$yhat_pcurves, type = "l", ylab = "Reference ASFR x period multiplier (per 1,000, log scale)",
     xlab = "Period", ylim = ylim.per, lwd = 1.25, log = "y", yaxt = "n")
axis(side = 2, at = ticks.per[ticks.per >= ylim.per[1] & ticks.per <= ylim.per[2]], las = 1, cex.axis = 0.9)
abline(h = asfr.ref, lty = "dashed")
lines(x = newdata$p, y = newdata$lo, col = adjustcolor("gray", alpha.f = 0.7))
lines(x = newdata$p, y = newdata$hi, col = adjustcolor("gray", alpha.f = 0.7))
invisible(dev.off())
mn_add("fert_ref_asfr", asfr.ref, 0, "Figure 7 note (optional): reference ASFR per 1,000")

# (9e) PREDICTED-RATES HEAT MAP FROM THE PREFERRED MODEL (APPENDIX E)
yhat.lcsc.vt <- data.frame(a.index = data$a.index, p.index = data$p.index, c.index = data$c.index,
                           yhat.rate = mu.fh / data$pop * 1000)
matAP.lcsc.vt <- mean_by_ap(yhat.lcsc.vt)
d <- matAP.lcsc.vt
z.min <- round(min(d, na.rm = TRUE), digits = 0)
z.max <- round(max(d, na.rm = TRUE), digits = 0)
plotAPCHeatmap(d, z.min, z.max, by.z = 10, save_pdf = TRUE,
               pdf_name = file.path(fig.dir, "2D_matAP.lcsc.vt.pdf"))
cat(sprintf("2D_matAP.lcsc.vt: predicted ASFR range %.1f-%.1f (observed %.1f-%.1f); max |predicted - observed| = %.2f per 1,000\n",
            min(d, na.rm = TRUE), max(d, na.rm = TRUE), min(matAP.raw, na.rm = TRUE), max(matAP.raw, na.rm = TRUE),
            max(abs(d - matAP.raw), na.rm = TRUE)))

#### 10. Wrap-Up ####

cat("\nFigures written to", fig.dir, ":\n")
print(list.files(fig.dir, pattern = "\\.pdf$"))
cat("Tables written to Output/: fertility_fit_table.csv, fertility_concurvity.csv, tex/fertility_fit_table.tex\n")
cat("Cached fits in Models/:", paste(list.files("Models", pattern = "^fert_"), collapse = ", "), "\n")
if (file.exists("Rplots.pdf")) warning("A stray Rplots.pdf was created in the working directory.")
cat(sprintf("Total elapsed: %.1f min\n", as.numeric(difftime(Sys.time(), t.start, units = "mins"))))
cat("03_fertility_japan.R finished:", format(Sys.time()), "\n")

## END OF R CODE
