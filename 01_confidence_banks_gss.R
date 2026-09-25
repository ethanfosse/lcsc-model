# 01_confidence_banks_gss.R -- Confidence in U.S. banks, GSS 1975-2022
# Fosse and Winship, "Varieties of Cross-Cohort Differentiation" (Sociological Science).
# Input: Data/gss2022.RData (not included; see README.md). Packages: mgcv, plot3D, RColorBrewer; speedglm optional.
# Run from this folder: Rscript --vanilla 01_confidence_banks_gss.R
# Writes the paper's figures to Figures/ and tables to Output/.

CRITERION_FIGURES <- "REML"   # smoothing criterion for the model behind all figures
CRITERION_TABLES  <- "ML"     # smoothing criterion for Table 1 and Table 2
LEGACY_UBRE       <- FALSE
NTHREADS          <- 4        # threads for mgcv::gam.control(nthreads = )
USE_SPEEDGLM      <- requireNamespace("speedglm", quietly = TRUE)

# ---- 1. Setup ----------------------------------------------------------------------- #

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

lcsc_session_header("01_confidence_banks_gss.R")
lcsc_require_data("Data/gss2022.RData", "GSS cumulative file")
dir.create(file.path("Figures", "Banks"), recursive = TRUE, showWarnings = FALSE)
dir.create("Output", recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages({
  library(mgcv)         # penalized-spline GAMs
  library(plot3D)       # image2D() / contour2D() for the Lexis heat maps
  library(RColorBrewer) # color palettes
})

cat("Settings: figures =", CRITERION_FIGURES, "| tables =", CRITERION_TABLES,
    "| LEGACY_UBRE =", LEGACY_UBRE, "| NTHREADS =", NTHREADS, "| USE_SPEEDGLM =", USE_SPEEDGLM, "\n")

# Figure CohortCareers.raw places its cohort labels with jitter()
set.seed(20240617)

# ---- 2. Data preparation (individual level, single years) ---------------------------- #

# loading the GSS cumulative file (creates the data frame `df`)
load(file.path("Data", "gss2022.RData"))
gss <- df; rm(df)

# All models below describe the UNWEIGHTED analysis sample

# binary outcome: 1 = "a great deal" or "only some" confidence, 0 = "hardly any"
gss$y <- ifelse(gss$confinan %in% c("a great deal", "only some"), 1,
                ifelse(gss$confinan == "hardly any", 0, NA))

# subsetting to the analysis variables
df <- subset(gss, select = c(y, age, year))

# age, period, and cohort on their natural scales
df$a <- as.numeric(df$age)   # range: 18-89
df$p <- df$year; df$year <- NULL
df$c <- df$p - df$a

# listwise deletion and a unique identifier
df <- na.omit(df)
df$id <- 1:nrow(df)

# raw cell mean of y within each observed (age, period) cell, merged back onto rows
df$yhat.raw <- ave(df$y, df$a, df$p, FUN = function(x) mean(x, na.rm = TRUE))

# index variables on the natural scale (used for graphics)
df$a.index <- df$a; df$p.index <- df$p; df$c.index <- df$p - df$a

# ---- CONFIGURATION: the only place the analysis constants are written --------------- #

cfg <- lcsc_config(df,
                   age_center    = 46,
                   period_center = 1997,
                   cohort_center = 1951,
                   period_full   = 1975:2022,
                   knots         = c(a = 20, p = 20, c = 20))

# centered modeling columns a, p, c, derived from cfg (c = p - a on the centered scale)
df <- lcsc_center(df, cfg)

n.obs <- nrow(df)
cat("Analysis sample: R =", n.obs, "respondents (paper: 44,735);",
    "ages", min(df$a.index), "-", max(df$a.index), ", periods",
    min(df$p.index), "-", max(df$p.index), "(", length(unique(df$p.index)), "survey years )\n")

# ---- 3. Five-year-grouped (categorical) data for the comparison models --------------- #

# subsetting to the analysis variables
df.agg <- subset(gss, select = c(y, age, year))
df.agg$period <- df.agg$year; df.agg$year <- NULL

# 5-year age groups (left endpoints as labels)
age.breaks <- seq(from = 15, to = 90, by = 5)
df.agg$a <- cut(df.agg$age, breaks = age.breaks, right = FALSE,
                labels = age.breaks[-length(age.breaks)])

# 5-year period groups (left endpoints as labels)
period.breaks <- seq(from = 1975, to = 2029, by = 5)
df.agg$p <- cut(df.agg$period, breaks = period.breaks, right = FALSE,
                labels = period.breaks[-length(period.breaks)])

# cohort groups implied by the age and period groups
df.agg$c <- as.numeric(as.character(df.agg$p)) - as.numeric(as.character(df.agg$a))
df.agg$c <- as.factor(df.agg$c)

# listwise deletion and a unique identifier
df.agg <- na.omit(df.agg)
df.agg$id <- 1:nrow(df.agg)

# --- orthogonal polynomial contrasts, orthogonalized to {intercept, linear}
residualize_on_linear <- function(n) {
  M <- contr.poly(n)
  L <- seq_len(n) - 0.5 * (n + 1)  # centered linear component
  M[, 1] <- L
  L <- cbind(rep(1, n), L)                          # intercept + linear
  Mstar <- M - L %*% solve(t(L) %*% L) %*% t(L) %*% M  # residualize on {1, linear}
  Mstar[, 1] <- L[, 2]                              # keep the linear component in column 1
  stopifnot(max(abs(t(L) %*% Mstar[, -1])) < 1e-10) # nonlinear columns orthogonal to {1, linear}
  Mstar
}

# age, period, and cohort contrast matrices
I <- length(levels(df.agg$a)); matAstar <- residualize_on_linear(I)
J <- length(levels(df.agg$p)); matPstar <- residualize_on_linear(J)
K <- length(levels(df.agg$c)); matCstar <- residualize_on_linear(K)

# assigning the contrasts to the grouped factors
contrasts(df.agg$a) <- matAstar
contrasts(df.agg$p) <- matPstar
contrasts(df.agg$c) <- matCstar

# design matrices for the nonlinear components (dropping intercept and linear column)
Ma.nonlin <- data.frame(model.matrix(~ a, data = df.agg))[, -c(1, 2)]
Mp.nonlin <- data.frame(model.matrix(~ p, data = df.agg))[, -c(1, 2)]
Mc.nonlin <- data.frame(model.matrix(~ c, data = df.agg))[, -c(1, 2)]

# adding the nonlinear components to the grouped dataset
df.agg <- cbind(df.agg, Ma.nonlin)
df.agg <- cbind(df.agg, Mp.nonlin)
df.agg <- cbind(df.agg, Mc.nonlin)

# extracting and adding the linear components
a.L <- data.frame(model.matrix(~ a, data = df.agg))[, c(2)]
p.L <- data.frame(model.matrix(~ p, data = df.agg))[, c(2)]
c.L <- data.frame(model.matrix(~ c, data = df.agg))[, c(2)]
df.agg <- cbind(df.agg, a.L)
df.agg <- cbind(df.agg, p.L)
df.agg <- cbind(df.agg, c.L)

# group factor identifying each cell of the age-period Lexis table
temp <- expand.grid(a = levels(df.agg$a), p = levels(df.agg$p))
temp$group <- rownames(temp); temp$group <- as.numeric(temp$group)
df.agg.old <- merge(df.agg, temp, by.x = c("a", "p"), by.y = c("a", "p"))
df.agg.old$group <- as.factor(df.agg.old$group)
contrasts(df.agg.old$group) <- contr.sum(length(levels(df.agg.old$group)))
df.agg <- df.agg.old

stopifnot(nrow(df.agg) == n.obs)

# ---- 4. Categorical models (LC-SC and saturated) ------------------------------------ #

cat("Fitting the categorical (5-year-grouped) LC-SC and saturated models ...\n")

# regression formula for the categorical LC-SC model
LHS <- paste("y ~ ", sep = " ")
a.nonlin <- colnames(Ma.nonlin)
p.nonlin <- colnames(Mp.nonlin)
c.nonlin <- colnames(Mc.nonlin)
RHS <- paste(c("a.L", "c.L", a.nonlin, p.nonlin, c.nonlin), collapse = " + ")
fm.noP <- as.formula(paste(LHS, RHS, sep = ""))

# categorical LC-SC model
m.cat <- glm(fm.noP, data = df.agg, family = "binomial")

# categorical saturated model (one parameter per 5-year Lexis cell)
RHS <- c("1 + group")
fm.LAPCH <- as.formula(paste(LHS, RHS, sep = ""))
m.LAPCH.cat <- glm(fm.LAPCH, data = df.agg, family = "binomial")

# ---- 5. Categorical APC deviations (nonlinearities) --------------------------------- #

# intercept of the categorical LC-SC model, on the response (probability) scale
intercept <- coef(m.cat)["(Intercept)"]
intercept <- exp(intercept) / (1 + exp(intercept))

# index values for the categorical age, period, and cohort groups
a.index <- as.numeric(as.character(sort(unique(df.agg$a))))
p.index <- as.numeric(as.character(sort(unique(df.agg$p))))
c.index <- as.numeric(as.character(sort(unique(df.agg$c))))

# contrast matrices as data frames, with names matching the model coefficients
age.contrasts <- data.frame(matAstar)
period.contrasts <- data.frame(matPstar)
cohort.contrasts <- data.frame(matCstar)
colnames(age.contrasts) <- gsub("\\^", "\\.", paste0("a", colnames(matAstar)))
colnames(period.contrasts) <- gsub("\\^", "\\.", paste0("p", colnames(matPstar)))
colnames(cohort.contrasts) <- gsub("\\^", "\\.", paste0("c", colnames(matCstar)))

# prediction frame for age deviations (only age nonlinearities are non-zero)
predframe.age <- age.contrasts
predframe.age["a.L"] <- 0
predframe.age[, colnames(period.contrasts)] <- 0
predframe.age[, colnames(cohort.contrasts)] <- 0

# prediction frame for period deviations
predframe.period <- period.contrasts
predframe.period["p.L"] <- 0
predframe.period[, colnames(age.contrasts)] <- 0
predframe.period[, colnames(cohort.contrasts)] <- 0

# prediction frame for cohort deviations
predframe.cohort <- cohort.contrasts
predframe.cohort["c.L"] <- 0
predframe.cohort[, colnames(age.contrasts)] <- 0
predframe.cohort[, colnames(period.contrasts)] <- 0

# (1) age deviations
yhat.age <- predict(m.cat, newdata = predframe.age, se.fit = TRUE, type = "response")
yhat.age <- as.data.frame(yhat.age)
b.age <- yhat.age[, 1] - intercept
se.age <- yhat.age[, 2]
low.ci.age <- b.age - 1.96 * se.age
high.ci.age <- b.age + 1.96 * se.age
AgeDeviations <- data.frame(Age = a.index, Deviation = b.age, SE = se.age,
                            `CI Lo` = low.ci.age, `CI Hi` = high.ci.age, check.names = FALSE)

# (2) period deviations
yhat.period <- predict(m.cat, newdata = predframe.period, se.fit = TRUE, type = "response")
yhat.period <- as.data.frame(yhat.period)
b.period <- yhat.period[, 1] - intercept
se.period <- yhat.period[, 2]
low.ci.period <- b.period - 1.96 * se.period
high.ci.period <- b.period + 1.96 * se.period
PeriodDeviations <- data.frame(Period = p.index, Deviation = b.period, SE = se.period,
                               `CI Lo` = low.ci.period, `CI Hi` = high.ci.period, check.names = FALSE)

# (3) cohort deviations
yhat.cohort <- predict(m.cat, newdata = predframe.cohort, se.fit = TRUE, type = "response")
yhat.cohort <- as.data.frame(yhat.cohort)
b.cohort <- yhat.cohort[, 1] - intercept
se.cohort <- yhat.cohort[, 2]
low.ci.cohort <- b.cohort - 1.96 * se.cohort
high.ci.cohort <- b.cohort + 1.96 * se.cohort
CohortDeviations <- data.frame(Cohort = c.index, Deviation = b.cohort, SE = se.cohort,
                               `CI Lo` = low.ci.cohort, `CI Hi` = high.ci.cohort, check.names = FALSE)

# ---- 6. Continuous models: LC-SC (REML for figures, ML for tables) and saturated ----- #

# (1) LC-SC model

# analysis dataset
data <- df[, c("a", "p", "c", "y", "a.index", "p.index", "c.index")]

# number of knots for age, period, and cohort (from the configuration above)
ak <- cfg$knots[["a"]]
pk <- cfg$knots[["p"]]
ck <- cfg$knots[["c"]]

# LC-SC regression formula
formula <- as.formula(y ~ a + c +
                        s(a, bs = "cr", k = ak - 1) +
                        s(p, bs = "cr", k = pk - 1) +
                        s(c, bs = "cr", k = ck - 1))

# the unfitted model (fit = FALSE)
G <- mgcv::gam(formula, family = "binomial", data = data, fit = FALSE)

# continuous LC-SC model for the FIGURES (smoothing parameters by REML)
cat("Fitting the continuous LC-SC GAM for the figures (", CRITERION_FIGURES, ") ...\n", sep = "")
t.fig <- system.time(
  rs <- fit_lcsc(G, data = data, cfg = cfg, method = CRITERION_FIGURES,
                 legacy_ubre = LEGACY_UBRE, nthreads = NTHREADS)
)
m.lcsc <- rs$mod
cat(sprintf("  figures model: method = %s, elapsed = %.1f s, total EDF = %.4f\n",
            m.lcsc$method, t.fig["elapsed"], sum(m.lcsc$edf)))

# continuous LC-SC model for the TABLES (smoothing parameters by ML)
cat("Fitting the continuous LC-SC GAM for the tables (", CRITERION_TABLES, ") ...\n", sep = "")
t.tab <- system.time(
  rs.ML <- fit_lcsc(G, data = data, cfg = cfg, method = CRITERION_TABLES,
                    legacy_ubre = LEGACY_UBRE, nthreads = NTHREADS)
)
m.lcsc.ML <- rs.ML$mod
cat(sprintf("  tables model:  method = %s, elapsed = %.1f s, total EDF = %.4f\n",
            m.lcsc.ML$method, t.tab["elapsed"], sum(m.lcsc.ML$edf)))

cat("\n==== summary(): continuous LC-SC model for the figures (", m.lcsc$method, ") ====\n", sep = "")
print(summary(m.lcsc))
cat("\n==== summary(): continuous LC-SC model for the tables (", m.lcsc.ML$method, ") ====\n", sep = "")
print(summary(m.lcsc.ML))

# (2) saturated model (one parameter per single-year age-period Lexis cell)

# group factor for each observed age-period cell, with sum-to-zero contrasts
df$aF <- as.factor(df$a); df$pF <- as.factor(df$p)
temp <- expand.grid(aF = levels(df$aF), pF = levels(df$pF))
temp$group <- rownames(temp); temp$group <- as.numeric(temp$group)
df.old <- merge(df, temp, by.x = c("aF", "pF"), by.y = c("aF", "pF"))
df.old$group <- as.factor(df.old$group)
contrasts(df.old$group) <- contr.sum(length(levels(df.old$group)))
df <- df.old

# fully saturated individual-level model
fit_saturated <- function(formula, data) {
  if (USE_SPEEDGLM) {
    # speedglm warns that it drops the sum-to-zero contrasts of `group`
    m <- suppressWarnings(speedglm::speedglm(formula, data = data, family = binomial(), fitted = TRUE))
    list(engine = "speedglm", logLik = as.numeric(logLik(m)), deviance = m$deviance,
         null.deviance = m$nulldev, n = m$n, rank = m$rank, coef = coef(m),
         AIC = AIC(m), BIC = BIC(m))
  } else {
    m <- glm(formula, data = data, family = binomial())
    list(engine = "glm", logLik = as.numeric(logLik(m)), deviance = m$deviance,
         null.deviance = m$null.deviance, n = nobs(m), rank = m$rank, coef = coef(m),
         AIC = AIC(m), BIC = BIC(m))
  }
}
cat("Fitting the continuous saturated model (one parameter per cell, ",
    if (USE_SPEEDGLM) "speedglm" else "glm", ") ...\n", sep = "")
t.sat <- system.time(m.LAPCH <- fit_saturated(y ~ group, data = df))
cat(sprintf("  saturated model: %d observed age-period cells, elapsed = %.1f s (%s)\n",
            length(levels(df$group)), t.sat["elapsed"], m.LAPCH$engine))

# ---- 7. Fit statistics: Table 1 of the paper ----------------------------------------- #

# degrees of freedom / number of parameters
p.LAPCH.cat <- ncol(model.matrix(m.LAPCH.cat))     # 5-year cells (paper: 150)
p.cat       <- ncol(model.matrix(m.cat))           # categorical LC-SC (paper: 46)
p.LAPCH     <- length(m.LAPCH$coef)                # single-year cells (= intercept + cells - 1)
edf.gam     <- sum(m.lcsc.ML$edf)                  # continuous LC-SC, ML fit
stopifnot(!anyNA(m.LAPCH$coef), p.LAPCH == m.LAPCH$rank,
          length(m.cat$y) == n.obs, nobs(m.lcsc.ML) == n.obs, m.LAPCH$n == n.obs)

# AIC and BIC from the log-likelihood and the EDF shown in the table
fit_row <- function(LL, dev, nulldev, n, edf) {
  data.frame(LLV = LL,
             R2_D = deviance_r2(dev, nulldev),
             Adj_R2_D = adj_deviance_r2(dev, nulldev, n, edf),
             AIC = -2 * LL + 2 * edf, BIC = -2 * LL + edf * log(n), EDF = edf)
}
tab1 <- rbind(
  fit_row(as.numeric(logLik(m.LAPCH.cat)), m.LAPCH.cat$deviance, m.LAPCH.cat$null.deviance, n.obs, p.LAPCH.cat),
  fit_row(as.numeric(logLik(m.cat)), m.cat$deviance, m.cat$null.deviance, n.obs, p.cat),
  fit_row(m.LAPCH$logLik, m.LAPCH$deviance, m.LAPCH$null.deviance, n.obs, p.LAPCH),
  fit_row(as.numeric(logLik(m.lcsc.ML)), m.lcsc.ML$deviance, m.lcsc.ML$null.deviance, n.obs, edf.gam)
)
stopifnot(isTRUE(all.equal(tab1$AIC[1:3], c(AIC(m.LAPCH.cat), AIC(m.cat), m.LAPCH$AIC))),
          isTRUE(all.equal(tab1$BIC[1:3], c(BIC(m.LAPCH.cat), BIC(m.cat), m.LAPCH$BIC))))
tab1 <- cbind(Data  = c("Categorical", "Categorical", "Continuous", "Continuous"),
              Model = c("Saturated model (m.LAPCH.cat)", "LC-SC model (m.cat)",
                        "Saturated model (m.LAPCH)", paste0("LC-SC model (m.lcsc.ML, ", m.lcsc.ML$method, ")")),
              tab1)
rownames(tab1) <- NULL

cat("\n==== Table 1: fit statistics of the LC-SC and saturated models (n = ", n.obs, ") ====\n", sep = "")
tab1.print <- tab1
tab1.print$LLV <- round(tab1$LLV, 2); tab1.print$R2_D <- round(tab1$R2_D, 4)
tab1.print$Adj_R2_D <- round(tab1$Adj_R2_D, 4); tab1.print$AIC <- round(tab1$AIC, 2)
tab1.print$BIC <- round(tab1$BIC, 2); tab1.print$EDF <- round(tab1$EDF, 2)
print(tab1.print, row.names = FALSE)
cat(sprintf(paste0("Note: AIC = -2 LLV + 2 EDF and BIC = -2 LLV + EDF log(R) with the EDF of the table.\n",
                   "      For the GAM, mgcv's AIC()/BIC() would use the corrected df = %.4f (Wood et al.\n",
                   "      2016) instead of sum(edf) = %.4f, giving AIC = %.2f and BIC = %.2f.\n"),
            attr(logLik(m.lcsc.ML), "df"), edf.gam, AIC(m.lcsc.ML), BIC(m.lcsc.ML)))
write.csv(tab1, file.path("Output", "banks_table1_fit_statistics.csv"), row.names = FALSE)

# Table 1 as LaTeX rows, formatted exactly as in the manuscript
tab1_rows <- function(t) {
  cell <- function(i) sprintf("$%s$ & $%s$ & $%s$ & $%s$ & $%s$ & $%s$",
                              fmt_tex(t$LLV[i], 2), fmt_tex(t$R2_D[i], 3), fmt_tex(t$Adj_R2_D[i], 3),
                              fmt_tex(t$AIC[i], 2), fmt_tex(t$BIC[i], 2), fmt_tex(t$EDF[i], 2))
  brace <- function(lab) sprintf("\\hspace{-0.25em}\\ldelim\\{{2.5}{*}[\\parbox{3cm-\\tabcolsep-\\widthof{$\\Big($}}{\\begin{tabular}{c} %s \\end{tabular} }]", lab)
  c(paste0(brace("Categorical"), " & Saturated model & ", cell(1), " \\\\[1.5ex]"),
    paste0("& LC-SC model & ", cell(2), " \\\\[2.5ex]"),
    paste0(brace("Continuous"), " & Saturated model & ", cell(3), " \\\\[1.5ex]"),
    paste0("& LC-SC model & ", cell(4), " \\\\[2.5ex]"))
}
write_tex(tab1_rows(tab1), file.path("Output", "tex", "banks_table1.tex"))
cat("LaTeX rows of Table 1 written to Output/tex/banks_table1.tex\n")

# numbers quoted in the text of Section IV.1
mn_add("banks_R", n.obs, 0, "sample size R, everywhere")
mn_add("banks_edf_cat_lcsc", p.cat, 0, "Section IV.1: categorical LC-SC model uses 46 of 150 df")
mn_add("banks_edf_cat_sat", p.LAPCH.cat, 0, "Section IV.1: 150 df of the categorical saturated model")
mn_add("banks_pct_cat", 100 * p.cat / p.LAPCH.cat, 1, "Section IV.1: 30.7 percent")
mn_add("banks_edf_gam", edf.gam, 2, "Section IV.1 and Table 1: EDF of the continuous LC-SC model (ML)")
mn_add("banks_edf_cont_sat", p.LAPCH, 0, "Section IV.1 and Table 1: number of single-year cells")
mn_add("banks_pct_gam", 100 * edf.gam / p.LAPCH, 1, "Section IV.1: 1.4 percent")

# deviance R2 and adjusted deviance R2 for the four models (also as a CSV)
r2.tab <- data.frame(
  Model = c("m.LAPCH.cat (categorical saturated)", "m.cat (categorical LC-SC)",
            "m.LAPCH (continuous saturated)", paste0("m.lcsc.ML (continuous LC-SC, ", m.lcsc.ML$method, ")")),
  n = n.obs,
  EDF = c(p.LAPCH.cat, p.cat, p.LAPCH, edf.gam),
  null_deviance = c(m.LAPCH.cat$null.deviance, m.cat$null.deviance, m.LAPCH$null.deviance, m.lcsc.ML$null.deviance),
  residual_deviance = c(m.LAPCH.cat$deviance, m.cat$deviance, m.LAPCH$deviance, m.lcsc.ML$deviance)
)
r2.tab$R2_D <- deviance_r2(r2.tab$residual_deviance, r2.tab$null_deviance)
r2.tab$Adj_R2_D <- adj_deviance_r2(r2.tab$residual_deviance, r2.tab$null_deviance, r2.tab$n, r2.tab$EDF)
cat("\n==== Deviance R2 and adjusted deviance R2 ====\n")
print(data.frame(r2.tab[, c("Model", "EDF")], R2_D = round(r2.tab$R2_D, 5), Adj_R2_D = round(r2.tab$Adj_R2_D, 5)),
      row.names = FALSE)

# ---- 8. LC-SC model summary: Table 2 of the paper ------------------------------------ #

# Fully penalized nonlinear terms require zero-dimensional penalty null spaces
if (!LEGACY_UBRE) {
  stopifnot(all(vapply(m.lcsc.ML$smooth, function(sm) sm$null.space.dim == 0,
                      logical(1))), !is.null(m.lcsc.ML$R))
  XtWX <- crossprod(rs.ML$X * sqrt(m.lcsc.ML$weights))
  stopifnot(max(abs(crossprod(m.lcsc.ML$R) - XtWX)) < 1e-8 * max(abs(XtWX)))
}
sm.ML <- summary(m.lcsc.ML, re.test = TRUE)

# parametric part
p.tab <- sm.ML$p.table
or.tab <- data.frame(
  term     = rownames(p.tab),
  paper    = c("mu", "theta1 (LC slope)", "theta2 (SC slope)"),
  estimate = p.tab[, "Estimate"],
  se       = p.tab[, "Std. Error"],
  OR       = exp(p.tab[, "Estimate"]),
  OR.lo95  = exp(p.tab[, "Estimate"] - 1.96 * p.tab[, "Std. Error"]),
  OR.hi95  = exp(p.tab[, "Estimate"] + 1.96 * p.tab[, "Std. Error"]),
  z        = p.tab[, "z value"],
  p.value  = p.tab[, "Pr(>|z|)"],
  row.names = NULL
)
cat("\n==== Table 2 (parametric terms, ", m.lcsc.ML$method, " fit): odds ratios and 95% CIs ====\n", sep = "")
print(data.frame(or.tab[, 1:2], round(or.tab[, 3:8], 4), p.value = signif(or.tab$p.value, 3)), row.names = FALSE)

# Test each nonlinear function against zero, conditional on the other model terms
s.tab <- sm.ML$s.table
smooth.tab <- data.frame(
  term   = rownames(s.tab),
  paper  = c("f_alpha~(a) age", "g_pi~(p) period", "h_gamma~(c) cohort"),
  EDF    = s.tab[, "edf"],
  Ref.df = s.tab[, "Ref.df"], Chi.sq = s.tab[, "Chi.sq"], p.value = s.tab[, "p-value"],
  row.names = NULL
)
cat("\n==== Table 2 (smooth terms, ", m.lcsc.ML$method, " fit): EDFs and test statistics ====\n", sep = "")
print(data.frame(EDF = round(smooth.tab$EDF, 4),
                 reference.rank = smooth.tab$Ref.df,
                 statistic = round(smooth.tab$Chi.sq, 4),
                 p.value = format.pval(smooth.tab$p.value, digits = 3, eps = 0.001),
                 row.names = smooth.tab$term))

# combined Table 2 as one CSV
na.p <- rep(NA_real_, nrow(or.tab)); na.s <- rep(NA_real_, nrow(smooth.tab))
tab2 <- data.frame(
  part     = c(rep("parametric", nrow(or.tab)), rep("smooth", nrow(smooth.tab))),
  term     = c(or.tab$term, smooth.tab$term),
  paper    = c(or.tab$paper, smooth.tab$paper),
  estimate = c(or.tab$estimate, na.s), se = c(or.tab$se, na.s),
  OR = c(or.tab$OR, na.s), OR.lo95 = c(or.tab$OR.lo95, na.s), OR.hi95 = c(or.tab$OR.hi95, na.s),
  z = c(or.tab$z, na.s), p.value = c(or.tab$p.value, smooth.tab$p.value),
  EDF = c(na.p, smooth.tab$EDF), Ref.df = c(na.p, smooth.tab$Ref.df), Chi.sq = c(na.p, smooth.tab$Chi.sq),
  method = m.lcsc.ML$method
)
write.csv(tab2, file.path("Output", "banks_table2_model_summary.csv"), row.names = FALSE)

# Table 2 as LaTeX rows, formatted exactly as in the manuscript
tab2_rows <- function(or, sm) {
  par.lab <- c("$\\mu$", "$\\theta_{1}$", "$\\theta_{2}$")
  sm.lab  <- c("$f_{\\widetilde{\\alpha}}$", "$g_{\\widetilde{\\pi}}$", "$h_{\\widetilde{\\gamma}}$")
  sapply(1:3, function(i) sprintf("%s & $%s%s$ & $(%s, %s)$ & %s & %s & $%s$ & $%s%s$ \\\\",
                                  par.lab[i], fmt_tex(or$OR[i], 3), sig_stars(or$p.value[i]),
                                  fmt_tex(or$OR.lo95[i], 3), fmt_tex(or$OR.hi95[i], 3),
                                  if (i == 1) "$\\quad$" else "", sm.lab[i],
                                  fmt_tex(sm$EDF[i], 2),
                                  fmt_tex(sm$Chi.sq[i], 2), sig_stars(sm$p.value[i])))
}
write_tex(tab2_rows(or.tab, smooth.tab), file.path("Output", "tex", "banks_table2.tex"))
cat("LaTeX rows of Table 2 written to Output/tex/banks_table2.tex\n")
for (i in 1:3) {
  mn_add(paste0("banks_OR_", c("mu", "theta1", "theta2")[i]), or.tab$OR[i], 3, "Table 2 / Section IV.1")
  mn_add(paste0("banks_edf_s", c("a", "p", "c")[i]), smooth.tab$EDF[i], 2, "Table 2")
}

# the figures model (REML)
edf.fig <- sapply(seq_along(m.lcsc$smooth), function(i) smooth_edf(m.lcsc, i))
names(edf.fig) <- sapply(m.lcsc$smooth, `[[`, "label")
cat("\n==== Figures model (", m.lcsc$method, "): smooth-term EDFs and total EDF ====\n", sep = "")
print(round(c(edf.fig, total = sum(m.lcsc$edf)), 4))
cat("(tables model, ", m.lcsc.ML$method, ": ", paste(sprintf("%s = %.4f", rownames(s.tab), s.tab[, "edf"]), collapse = ", "),
    ", total = ", round(edf.gam, 4), ")\n", sep = "")

# ---- 9. Concurvity of the LC-SC model (Appendix F) ----------------------------------- #

conc.full.ML <- concurvity_lcsc(m.lcsc.ML, rs.ML$X, full = TRUE)
conc.pair.ML <- concurvity_lcsc(m.lcsc.ML, rs.ML$X, full = FALSE)
cat("\n==== Concurvity, orthogonalized design (correct), tables model (", m.lcsc.ML$method,
    "): overall (full = TRUE) ====\n", sep = "")
print(round(conc.full.ML, 4))
cat("\n==== Concurvity, orthogonalized design (correct), tables model (", m.lcsc.ML$method,
    "): pairwise (full = FALSE) ====\n", sep = "")
for (msr in names(conc.pair.ML)) { cat("-- ", msr, " --\n", sep = ""); print(round(conc.pair.ML[[msr]], 4)) }
for (msr in c("worst", "observed", "estimate")) for (tm in c("s(a)", "s(p)", "s(c)"))
  mn_add(paste0("banks_conc_", msr, "_", gsub("[()]", "", tm)), conc.full.ML[msr, tm], 3,
         "Appendix F (multivariate concurvity, ML model, fitted design)")
# pairwise worst-case concurvity of the cohort smooth with the intercept
for (tm in c("para", "s(a)", "s(p)"))
  mn_add(paste0("banks_conc_pair_worst_sc_", gsub("[()]", "", tm)), conc.pair.ML$worst["s(c)", tm], 3,
         paste("Appendix F: pairwise worst-case concurvity of s(c) with", c(para = "the intercept", "s(a)" = "s(a) (optional)", "s(p)" = "s(p) (optional)")[[tm]]))

# the same measures on an explicit partition of the columns
bl.ML <- lcsc_blocks(m.lcsc.ML, rs.ML$X, extra = list(lc_slope = "^a$", sc_slope = "^c$"))
cb.ML <- concurvity_blocks(rs.ML$X, coef(m.lcsc.ML), bl.ML, full = TRUE)
stopifnot(max(abs(cb.ML[, c("s(a)", "s(p)", "s(c)")] - conc.full.ML[, c("s(a)", "s(p)", "s(c)")])) < 1e-8)
cat("\n==== Concurvity, orthogonalized design, tables model: intercept and slopes as separate blocks ====\n")
print(round(cb.ML, 4))
cb.df <- do.call(rbind, lapply(colnames(cb.ML), function(tm) data.frame(
  model = "LC-SC (ML)", metric = "unweighted", type = "full", measure = rownames(cb.ML), term = tm,
  term_given = NA_character_, value = cb.ML[, tm], ncol = attr(cb.ML, "ncol")[[tm]],
  rank = attr(cb.ML, "rank")[[tm]], stringsAsFactors = FALSE)))
write.csv(cb.df, file.path("Output", "banks_concurvity_blocks.csv"), row.names = FALSE)
for (msr in c("worst", "observed", "estimate")) for (tm in c("intercept", "lc_slope", "sc_slope"))
  mn_add(paste0("banks_conc_", msr, "_", c(intercept = "int", lc_slope = "lcs", sc_slope = "scs")[[tm]]), cb.ML[msr, tm], 3,
         "Appendix F (multivariate concurvity, ML model, fitted design)")

conc.full <- concurvity_lcsc(m.lcsc, rs$X, full = TRUE)
conc.pair <- concurvity_lcsc(m.lcsc, rs$X, full = FALSE)
cat("\n==== Concurvity, orthogonalized design (correct), figures model (", m.lcsc$method,
    "): overall (full = TRUE) ====\n", sep = "")
print(round(conc.full, 4))
cat("\n==== Concurvity, orthogonalized design (correct), figures model (", m.lcsc$method,
    "): pairwise (full = FALSE) ====\n", sep = "")
for (msr in names(conc.pair)) { cat("-- ", msr, " --\n", sep = ""); print(round(conc.pair[[msr]], 4)) }

conc.full.mgcv <- mgcv::concurvity(m.lcsc, full = TRUE)
conc.pair.mgcv <- mgcv::concurvity(m.lcsc, full = FALSE)
cat("\n==== Concurvity, mgcv default on the un-orthogonalized design -- NOT appropriate for these\n",
    "     models; shown only for comparison with earlier versions of the analysis: overall ====\n", sep = "")
print(round(conc.full.mgcv, 4))
cat("\n==== (same, mgcv default, un-orthogonalized design): pairwise ====\n")
for (msr in names(conc.pair.mgcv)) { cat("-- ", msr, " --\n", sep = ""); print(round(conc.pair.mgcv[[msr]], 4)) }

# ---- 10. Plot helper functions ------------------------------------------------------- #

# mean_by_ap(): / mean_by_pc() / mean_by_ac()

# full single-year ranges of age, period, and cohort
a.index <- min(df$a.index):max(df$a.index)
p.index <- cfg$period_full
c.index <- min(df$c.index):max(df$c.index)

# axis label positions of the heat maps, derived from the index ranges
age.lab <- lcsc_axis_offsets(a.index, by = 5, origin_adj = 2)
per.lab <- lcsc_axis_offsets(p.index, by = 5, origin_adj = 1)
coh.lab <- lcsc_axis_offsets(c.index, by = 5, origin_adj = 1)

# Lexis heat map with contour lines (age x period surface with cohort diagonals)
plotAPCHeatmap <- function(d, z.min = NULL, z.max = NULL, by.z = NULL, save_pdf = TRUE,
                           breaks = seq(0, 1, length.out = 100),
                           pdf_name = "2D_APC_heatmap.pdf") {

  # color palette
  mypal <- rev(colorRampPalette(brewer.pal(9, "RdYlBu"), alpha = 1, bias = 1)(200))

  if (is.null(by.z)) {
    by.z <- round(diff(seq(from = z.min, to = z.max, length.out = 15))[1], digits = 1)
  }
  if (is.null(z.min)) {
    z.min <- round(min(d, na.rm = TRUE), digits = 0)
  }
  if (is.null(z.max)) {
    z.max <- round(max(d, na.rm = TRUE), digits = 0)
  }

  if (save_pdf) {
    pdf(pdf_name, width = 11.25, height = 9.5)
  }

  # expanding the margins
  par(mar = c(5.4, 4.1, 4.1, 5.1))

  # orienting the matrix for image2D (reverse rows, then transpose)
  d <- d[nrow(d):1, ]
  d <- t(d)
  d.contour <- d

  # age labels (left), period labels (top), cohort labels (bottom and right)
  age.plot.labs <- age.lab$labs;    age.loc <- age.lab$loc
  period.plot.labs <- per.lab$labs; period.loc <- per.lab$loc
  cohort.plot.labs <- coh.lab$labs

  # heat map plus contour lines
  plot3D::image2D(z = d, x = 1:nrow(d), y = 1:ncol(d), shade = 0.01, rasterImage = TRUE,
                  col = mypal, colkey = FALSE, axes = F, ylab = "", xlab = "", breaks = breaks)
  plot3D::contour2D(z = d.contour, x = 1:nrow(d.contour), y = 1:ncol(d.contour), col = "black",
                    labcex = 0.5, lwd = 1, alpha = 0.8,
                    levels = seq(z.min, z.max, b = by.z), add = TRUE)

  # axis labels: age (left), period (top), cohort (bottom and right)
  text(x = 0.5, y = age.loc, pos = 2, srt = 0, labels = rev(age.plot.labs), xpd = TRUE, cex = 0.6)
  text(x = period.loc, y = length(a.index) + 0.5, pos = 3, srt = 0, labels = period.plot.labs, xpd = TRUE, cex = 0.6)
  text(x = period.loc + 0.25, y = -3.00, pos = 3, srt = 0, srt = 320,
       labels = cohort.plot.labs[1:length(period.loc)], xpd = TRUE, cex = 0.6)
  text(x = length(p.index) + 0.65, y = age.loc, pos = 4, srt = 320,
       labels = cohort.plot.labs[length(period.loc):length(cohort.plot.labs)] + 2, xpd = TRUE, cex = 0.6)

  # tick marks
  axis(side = 2, at = age.loc, tck = -0.01, labels = F)
  axis(side = 1, at = period.loc, tck = -0.015, labels = F)
  axis(side = 3, at = period.loc, tck = -0.01, labels = F)
  axis(side = 4, at = age.loc, tck = -0.015, labels = F)

  # axis titles
  mtext(side = 1, "Cohort", line = 2.25, cex = 0.8)
  mtext(side = 2, "Age", line = 2.25, cex = 0.8)
  mtext(side = 3, "Period", line = 2, cex = 0.8)
  text(x = length(p.index) + 4.5, y = mean(seq(1:length(a.index))), xpd = T, labels = "Cohort", cex = 0.8, srt = 270)

  if (save_pdf) {
    invisible(dev.off())
  } else {
    # returning to default margins (only meaningful when drawing on an open device)
    par(mar = c(5.1, 4.1, 4.1, 2.1))
  }
  invisible(NULL)
}

# conventional cohort careers
cohort_careers_plot_full <- function(matAC, a_index, c_index, xlim = c(15, 95), ylim = c(0, 1),
                                     output_filename = "FullCohortCareers.pdf") {
  pdf(output_filename, width = 11.5, height = 10.5)

  par(mfrow = c(1, 1))

  # matrix of cohort careers (rows become cohorts after transposing)
  d <- matAC; xlim <- c(15, 100); ylim <- c(0, 1) # the published figure's axes
  d <- t(d)

  # blank plot
  plot(x = a_index, y = rep(NA, length(a_index)), col = "black",
       type = "n", ylim = ylim, xlim = xlim,
       xlab = "Age", ylab = "Pr(Confidence)", yaxt = "n", xaxt = "n",
       main = " ", cex.main = 0.9)

  # axis ticks
  axis(side = 2, at = seq(from = 0, to = 1, by = 0.2), las = 1, cex.axis = 0.8)
  axis(side = 1, at = seq(from = 15, to = 90, by = 15), las = 0, cex.axis = 0.8)

  # cohorts plotted
  c_index_subset <- c_index[-c(1:4, (length(c_index) - 1):(length(c_index)))]
  selected_cohorts <- c_index_subset[seq(1, length(c_index_subset), 1)] # every cohort
  c_index_number <- which(c_index %in% selected_cohorts)

  # grey palette with transparency
  mypal <- colorRampPalette(rev(brewer.pal(9, "Greys")))(110)
  mypal <- adjustcolor(mypal, alpha.f = 0.6)

  # one line (and jittered end label) per cohort
  for (i in 1:length(c_index_number)) {
    c_ind <- c_index_number[i]
    lines(x = a_index, d[c_ind, ], col = mypal[i],
          lty = 1, lwd = 1)
    x_val <- as.numeric(names(na.omit(d[c_ind, ]))[length(names(na.omit(d[c_ind, ])))])
    text(x = x_val + 1,
         y = jitter(na.omit(d[c_ind, ])[length(na.omit(d[c_ind, ]))], factor = 0.5),
         labels = rownames(d)[c_ind],
         cex = 0.45, col = mypal[i])
  }

  invisible(dev.off())
}

# Lexis heat map without contours, with per-cell value labels (raw cell means)
plotAPCHeatmap_nocontour <- function(d, z.min = NULL, z.max = NULL, by.z = NULL, save_pdf = TRUE,
                                     breaks = NULL, pdf_name = "2D_APC_heatmap.pdf") {

  # color palette
  mypal <- colorRampPalette(rev(brewer.pal(9, "RdYlBu")[-1]), bias = 0.3)(150)

  if (is.null(by.z)) {
    by.z <- round(diff(seq(from = z.min, to = z.max, length.out = 15))[1], digits = 1)
  }
  if (is.null(z.min)) {
    z.min <- round(min(d, na.rm = TRUE), digits = 0)
  }
  if (is.null(z.max)) {
    z.max <- round(max(d, na.rm = TRUE), digits = 0)
  }

  if (save_pdf) {
    pdf(pdf_name, width = 11.5, height = 11.5)
  }

  # expanding the margins
  par(mar = c(5.4, 4.1, 4.1, 5.1))

  # orienting the matrix for image2D (reverse rows, then transpose)
  d <- d[nrow(d):1, ]
  d <- t(d)

  # age labels (left), period labels (top), cohort labels (bottom and right)
  age.plot.labs <- age.lab$labs;    age.loc <- age.lab$loc
  period.plot.labs <- per.lab$labs; period.loc <- per.lab$loc
  cohort.plot.labs <- coh.lab$labs

  # background color for missing cells
  verylightgray <- gray(0.95)

  # padding with NA rows/columns so border cells are drawn fully
  d <- cbind(NA, d, NA)
  d <- rbind(NA, d, NA)

  # heat map
  if (is.null(breaks)) {
    plot3D::image2D(z = d, x = 1:nrow(d), y = 1:ncol(d), shade = 0, rasterImage = F,
                    col = mypal, colkey = FALSE, axes = F, ylab = "", xlab = "", border = "lightgray",
                    facets = T, NAcol = verylightgray)
  } else {
    plot3D::image2D(z = d, x = 1:nrow(d), y = 1:ncol(d), shade = 0, rasterImage = F,
                    col = mypal, colkey = FALSE, axes = F, ylab = "", xlab = "", border = "lightgray",
                    breaks = breaks, facets = T, NAcol = verylightgray)
  }

  # per-cell value labels
  for (i in seq_along(1:nrow(d))) {
    for (j in seq_along(1:ncol(d))) {
      cell_text <- ifelse(is.na(d[i, j]), "", sprintf("%.1f", d[i, j]))
      text(x = i, y = j, labels = cell_text, cex = 0.57, col = "black")
    }
  }

  # axis labels: age (left), period (top), cohort (bottom and right)
  text(x = 0.5, y = age.loc + 1, pos = 2, srt = 0, labels = rev(age.plot.labs), xpd = TRUE, cex = 0.6)
  text(x = period.loc + 1, y = length(a.index) + 2.5, pos = 3, srt = 0, labels = period.plot.labs, xpd = TRUE, cex = 0.6)
  text(x = period.loc + 1.25, y = -3.00, pos = 3, srt = 0, srt = 320,
       labels = cohort.plot.labs[1:length(period.loc)], xpd = TRUE, cex = 0.6)
  text(x = length(p.index) + 2.65, y = age.loc + 1, pos = 4, srt = 320,
       labels = cohort.plot.labs[length(period.loc):length(cohort.plot.labs)] + 2, xpd = TRUE, cex = 0.6)

  # tick marks
  axis(side = 2, at = age.loc + 1, tck = -0.01, labels = F)
  axis(side = 1, at = period.loc + 1, tck = -0.015, labels = F)
  axis(side = 3, at = period.loc + 1, tck = -0.01, labels = F)
  axis(side = 4, at = age.loc + 1, tck = -0.015, labels = F)

  # axis titles
  mtext(side = 1, "Cohort", line = 2.25, cex = 0.8)
  mtext(side = 2, "Age", line = 2.25, cex = 0.8)
  mtext(side = 3, "Period", line = 2, cex = 0.8)
  text(x = length(p.index) + 4.5, y = mean(seq(1:length(a.index))), xpd = T, labels = "Cohort", cex = 0.8, srt = 270)

  if (save_pdf) {
    invisible(dev.off())
  } else {
    # returning to default margins (only meaningful when drawing on an open device)
    par(mar = c(5.1, 4.1, 4.1, 2.1))
  }
  invisible(NULL)
}

# comparative ("Tuftean") cohort careers
cohort_careers_tuftean_plot_gap <- function(matAC_apc, a_index, c_index, yhat_age, xlim = c(15, 100), ylim = c(0, 1), by.y = NULL,
                                            output_filename = "CohortCareersTuftean.pdf") {
  pdf(output_filename, width = 10, height = 7.25)

  par(mfrow = c(1, 1))

  # matrix of cohort careers (rows become cohorts after transposing)
  d <- matAC_apc
  d <- t(d)

  c_index_number <- 1:length(c_index)

  # grey palette with transparency for the background careers
  mypal <- colorRampPalette(brewer.pal(9, "Greys"), alpha = FALSE, bias = 3)(118)
  mypal_transp <- paste0(mypal, "20")

  plot(x = a_index, y = a_index,
       col = "black",
       type = "n", ylim = ylim, xlim = xlim,
       xlab = "Age", ylab = "Pr(Confidence)", yaxt = "n", xaxt = "n",
       main = " ", cex.main = 0.9)

  if (is.null(by.y)) { by.y <- 1 }

  # axis ticks
  axis(side = 2, at = seq(from = 0, to = 10, by = by.y), las = 1, cex.axis = 0.8)
  axis(side = 1, at = seq(from = 15, to = 100, by = 15), las = 0, cex.axis = 0.8)

  # background careers
  for (i in 1:length(c_index_number)) {
    if (i %% 2 != 1) next

    c_ind <- c_index_number[i]

    lines(x = a_index, d[c_ind, ], col = mypal_transp[i],
          lty = 1, lwd = 1)
    x_val <- as.numeric(names(na.omit(d[c_ind, ]))[length(names(na.omit(d[c_ind, ])))])
  }

  # highlighted cohorts overlaid in color
  selected_cohorts <- c(1900, 1915, 1930, 1945, 1960, 1975, 1990)
  c_index_number <- which(c_index %in% selected_cohorts)

  mypal <- c("#9B59B6FF", "#2980B9FF", "#1ABC9CFF", "#E74C3CFF", "#34495EFF", "#E67E22FF", "#9A7D0AFF", "#641E16FF")

  for (i in 1:length(c_index_number)) {
    c_ind <- c_index_number[i]

    lines(x = a_index, d[c_ind, ], col = mypal[i], lty = 1, lwd = 1.3)

    x_val <- as.numeric(names(na.omit(d[c_ind, ]))[length(names(na.omit(d[c_ind, ])))])
    text(x = x_val + 2,
         y = na.omit(d[c_ind, ])[length(na.omit(d[c_ind, ]))],
         labels = rownames(d)[c_ind],
         cex = 0.55, col = mypal[i])
  }

  legend(x = 95, y = 0.65,
         box.lty = 0, legend = rownames(d)[c_index_number],
         col = mypal, lty = rep("solid", 9), cex = 0.6)

  invisible(dev.off())
}

# linear interpolation of missing raw cell means within cohorts (by age)
interpolate_missing_values <- function(data) {

  # unique cohort values
  unique_c_values <- unique(data$c.index)

  # results container
  interpolated_data <- data.frame()

  for (c_val in unique_c_values) {
    # subset for this cohort, sorted by age
    data_subset <- data[data$c.index == c_val, ]
    data_subset <- data_subset[order(data_subset$a.index), ]

    # first and last non-missing values
    first_non_na <- min(which(!is.na(data_subset$yhat.raw)))
    last_non_na <- max(which(!is.na(data_subset$yhat.raw)))

    # interpolate between the first and last non-missing values only
    y_interpolated <- data_subset$yhat.raw
    if (first_non_na < last_non_na) {
      y_interpolated[first_non_na:last_non_na] <- approx(data_subset$a.index[first_non_na:last_non_na],
                                                         data_subset$yhat.raw[first_non_na:last_non_na],
                                                         xout = data_subset$a.index[first_non_na:last_non_na],
                                                         method = "linear", rule = 2)$y
    }

    data_subset$yhat.raw.interpol <- y_interpolated
    interpolated_data <- rbind(interpolated_data, data_subset)
  }

  return(interpolated_data)
}

# panel-drawing versions of the comparative cohort-careers plot
cohort_careers_panel <- function(matAC_apc, a_index, c_index, ylab, xlim = c(15, 100), ylim = c(0.5, 1),
                                 by.y = NULL, panel_label = " ") {

  # matrix of cohort careers (rows become cohorts after transposing)
  d <- matAC_apc
  d <- t(d)

  c_index_number <- 1:length(c_index)

  # grey palette with transparency for the background careers
  mypal <- colorRampPalette(brewer.pal(9, "Greys"), alpha = FALSE, bias = 3)(118)
  mypal_transp <- paste0(mypal, "40")

  plot(x = a_index, y = a_index,
       col = "black",
       type = "n", ylim = ylim, xlim = xlim,
       xlab = "Age", ylab = ylab, yaxt = "n", xaxt = "n",
       main = panel_label, cex.main = 1.3)

  if (is.null(by.y)) { by.y <- 1 }

  # axis ticks
  axis(side = 2, at = seq(from = 0, to = 10, by = by.y), las = 1, cex.axis = 0.8)
  axis(side = 1, at = seq(from = 15, to = 100, by = 15), las = 0, cex.axis = 0.8)

  # background careers
  for (i in 1:length(c_index_number)) {
    c_ind <- c_index_number[i]

    lines(x = a_index, d[c_ind, ], col = mypal_transp[i],
          lty = 1, lwd = 1)
    x_val <- as.numeric(names(na.omit(d[c_ind, ]))[length(names(na.omit(d[c_ind, ])))])
  }

  # highlighted cohorts overlaid in color
  selected_cohorts <- c(1900, 1915, 1930, 1945, 1960, 1975, 1990)
  c_index_number <- which(c_index %in% selected_cohorts)

  mypal <- c("#9B59B6FF", "#2980B9FF", "#1ABC9CFF", "#E74C3CFF", "#34495EFF", "#E67E22FF", "#9A7D0AFF", "#641E16FF")

  for (i in 1:length(c_index_number)) {
    c_ind <- c_index_number[i]

    lines(x = a_index, d[c_ind, ], col = mypal[i], lty = 1, lwd = 1.3)

    x_val <- as.numeric(names(na.omit(d[c_ind, ]))[length(names(na.omit(d[c_ind, ])))])
    text(x = x_val + 2,
         y = na.omit(d[c_ind, ])[length(na.omit(d[c_ind, ]))],
         labels = rownames(d)[c_ind],
         cex = 0.55, col = mypal[i])
  }

  legend(x = 95, y = .75,
         box.lty = 0, legend = rownames(d)[c_index_number],
         col = mypal, lty = rep("solid", 9), cex = 0.6)
}

# prediction grid for the continuous LC-SC model
make_prediction_grid <- function() {
  newdata <- expand.grid(a.index = a.index, p.index = p.index)
  newdata$c.index <- newdata$p.index - newdata$a.index
  newdata <- lcsc_center(newdata, cfg)
  newdata[, c("a.index", "p.index", "c.index", "a", "p", "c")]
}

# coefficient blocks of the figures model, located by NAME rather than by position
ix   <- lcsc_coef_index(m.lcsc)
i.mu <- ix[["(Intercept)"]]; i.a <- ix[["a"]]; i.c <- ix[["c"]]
i.sa <- ix[["s(a)"]];        i.sp <- ix[["s(p)"]]; i.sc <- ix[["s(c)"]]
stopifnot(identical(i.sa, rs$ortho[[1]]$cols), identical(i.sp, rs$ortho[[2]]$cols),
          identical(i.sc, rs$ortho[[3]]$cols))

# ---- 11. Raw cell means on the full Lexis grid --------------------------------------- #

cat("Drawing the figures (Figures/Banks/) ...\n")

# full grid of single-year age-period (and implied cohort) combinations
df1 <- expand.grid(a.index = a.index, p.index = p.index); df1$c.index <- df1$p.index - df1$a.index
df1$a <- df1$a.index - mean(df1$a.index); df1$p <- df1$p.index - mean(df1$p.index)
df1$c <- df1$p - df1$a

# observed raw cell means (with the Lexis-cell group id)
df2 <- df[, c("a.index", "p.index", "group", "yhat.raw")]
df2$group <- as.numeric(as.character(df2$group))
df2_agg <- aggregate(cbind(group, yhat.raw) ~ a.index + p.index, df2, mean)
newdata <- merge(df1, df2_agg, by = c("a.index", "p.index"), all.x = TRUE)
newdata$group <- as.factor(newdata$group)

# linear interpolation of missing cells within cohorts
newdata <- interpolate_missing_values(newdata)

# ---- 12. Figure: CohortCareers.raw.pdf ----------------------------------------------- #

# interpolated cell means (4th column feeds the matrix constructors)
yhat.raw.complete <- newdata[, c("a.index", "p.index", "c.index", "yhat.raw.interpol", "a", "c", "p")]
matAP.apc_gam <- mean_by_ap(yhat.raw.complete); matPC.apc_gam <- mean_by_pc(yhat.raw.complete); matAC.apc_gam <- mean_by_ac(yhat.raw.complete)

cohort_careers_plot_full(matAC = matAC.apc_gam, a_index = a.index, c_index = c.index,
                         output_filename = file.path("Figures", "Banks", "CohortCareers.raw.pdf"))

# ---- 13. Figure: 2D_matAP.raw.pdf ----------------------------------------------------- #

# raw (non-interpolated) cell means
yhat.raw <- newdata[, c("a.index", "p.index", "c.index", "yhat.raw", "a", "c", "p")]
matAP.apc_gam <- mean_by_ap(yhat.raw); matPC.apc_gam <- mean_by_pc(yhat.raw); matAC.apc_gam <- mean_by_ac(yhat.raw)

d <- matAP.apc_gam
plotAPCHeatmap_nocontour(d, z.min = 0.5, z.max = 0.9, by.z = 5, save_pdf = TRUE,
                         pdf_name = file.path("Figures", "Banks", "2D_matAP.raw.pdf"))

# raw age-period matrix saved for the eta figure below
matAP.raw <- matAP.apc_gam

# ---- 14. Figure: 2D_matAP.apc_lcsc.pdf ------------------------------------------------ #

# prediction grid
newdata <- make_prediction_grid()

# LC-SC predictions on the grid
H <- lcsc_predict(rs, newdata)
newdata$yhat <- as.numeric(H[["Yhat"]][["est"]])

# fitted-value matrices
yhat <- newdata[, c("a.index", "p.index", "c.index", "yhat")]
matAP.apc_gam <- mean_by_ap(yhat); matPC.apc_gam <- mean_by_pc(yhat); matAC.apc_gam <- mean_by_ac(yhat)

d <- round(matAP.apc_gam, digits = 3)
plotAPCHeatmap(d, z.min = 0.55, z.max = 0.965, by.z = 0.05, save_pdf = TRUE,
               breaks = seq(0.55, 0.965, length.out = 201),
               pdf_name = file.path("Figures", "Banks", "2D_matAP.apc_lcsc.pdf"))

# ---- 15. Figure: CohortCareersTuftean.apc_lcsc.pdf ------------------------------------ #

cohort_careers_tuftean_plot_gap(matAC_apc = matAC.apc_gam, a_index = a.index, ylim = c(0.5, 1), by.y = 0.05,
                                c_index = c.index, yhat_age = a.index,
                                output_filename = file.path("Figures", "Banks", "CohortCareersTuftean.apc_lcsc.pdf"))

# fitted age-period matrix saved for the eta figure below
matAP.lcsc <- matAP.apc_gam

# ---- 16. Figure: Logit_TwoPanel.pdf --------------------------------------------------- #

# prediction grid (identical to the one above)
newdata <- make_prediction_grid()

# term-by-term predictions (link scale)
H <- lcsc_predict(rs, newdata)
terms <- H[["est"]]
stopifnot(identical(unname(colnames(terms)), c("a", "c", "s(a)", "s(p)", "s(c)")))

# linear predictor excluding s(p), plus the intercept
intercept <- H$intercept
linear_predictor <- intercept + rowSums(terms[, colnames(terms) != "s(p)"])

# probabilities via the logistic function
predictions_without_sp <- exp(linear_predictor) / (1 + exp(linear_predictor))
newdata$yhat <- predictions_without_sp
newdata$yhat_lp <- linear_predictor

# age-cohort matrices (mean_by_ac() reads the 4th column)
yhat <- newdata[, c("a.index", "p.index", "c.index", "yhat", "a", "c", "p", "yhat_lp")]
matAC.pnonzero.prob <- mean_by_ac(yhat)   # probability scale
yhat <- newdata[, c("a.index", "p.index", "c.index", "yhat_lp", "a", "c", "p", "yhat")]
matAC.pnonzero.logit <- mean_by_ac(yhat)  # log-odds scale

# two stacked panels in one pdf
pdf(file.path("Figures", "Banks", "Logit_TwoPanel.pdf"), width = 10, height = 14.5)
par(mfrow = c(2, 1))
cohort_careers_panel(matAC_apc = matAC.pnonzero.logit, a_index = a.index, c_index = c.index,
                     ylab = "Log-odds of Confidence", ylim = c(0, 3.2), by.y = 0.2, panel_label = "(a) Log-Odds Scale")
cohort_careers_panel(matAC_apc = matAC.pnonzero.prob, a_index = a.index, c_index = c.index,
                     ylab = "Pr(Confidence)", ylim = c(0.5, 1), by.y = 0.05, panel_label = "(b) Probability Scale")
invisible(dev.off())

# ---- 17. Figure: 2D_matAP.eta.pdf ----------------------------------------------------- #

matAP.eta <- matAP.raw - matAP.lcsc
matAP.eta <- round(matAP.eta, digits = 1)

d <- matAP.eta
z.min <- round(min(d, na.rm = TRUE), digits = 2)
z.max <- round(max(d, na.rm = TRUE), digits = 2)
plotAPCHeatmap_nocontour(d, z.min, z.max, by.z = 0.05, save_pdf = TRUE,
                         breaks = seq(-0.9, 0.3, length.out = 151),
                         pdf_name = file.path("Figures", "Banks", "2D_matAP.eta.pdf"))

# ---- 18. Figure: ThreePanel_APC_Gam.pdf ----------------------------------------------- #

# prediction grid (identical to the one above)
newdata <- make_prediction_grid()

# design matrix for the grid (H is reused by the three figures that follow)
H <- lcsc_predict(rs, newdata)
stopifnot(ncol(H$X) == length(coef(m.lcsc)))

# one component (a set of design columns) as predicted probabilities with 95% bounds
component_curve <- function(cols, sort_col = NULL, sort_by = NULL) {
  sub <- H$X[, cols, drop = FALSE]
  keep <- which(!duplicated(sub))
  key <- if (is.null(sort_by)) sub[keep, sort_col] else sort_by[keep]
  rows <- keep[order(key)]
  comp <- lcsc_component(rs, H$X[rows, , drop = FALSE], cols)
  data.frame(pred = plogis(comp$lp), lo = plogis(comp$lo), hi = plogis(comp$hi))
}

pdf(file.path("Figures", "Banks", "ThreePanel_APC_Gam.pdf"), width = 12, height = 4.75)
par(mfrow = c(1, 3))

## PANEL (a): THE LC CURVE (INTERCEPT + AGE LINEAR + AGE NONLINEARITIES)

new_data <- cbind(a.index = a.index, component_curve(c(i.mu, i.a, i.sa), sort_col = 2))  # sorted by "a"

plot(x = new_data$a, y = new_data$pred, type = "l", ylim = c(0.4, 1),
     xlab = "Age", ylab = "Pr(Confidence)", xlim = c(15, 95),
     yaxt = "n", xaxt = "n", main = "(a) LC Curve", cex.main = 1.5)
axis(side = 2, at = seq(from = 0, to = 10, by = 0.1), las = 1, cex.axis = 0.8)
axis(side = 1, at = seq(from = 15, to = 100, by = 15), las = 0, cex.axis = 0.8)
lines(x = new_data$a, y = new_data$lo, col = adjustcolor("gray", alpha.f = 0.7))
lines(x = new_data$a, y = new_data$hi, col = adjustcolor("gray", alpha.f = 0.7))

## PANEL (b): THE SC CURVE (INTERCEPT + COHORT LINEAR + COHORT NONLINEARITIES)

new_data <- cbind(c.index = c.index, component_curve(c(i.mu, i.c, i.sc), sort_col = 2))  # sorted by "c"

plot(x = new_data$c, y = new_data$pred, type = "l", ylim = c(0.4, 1),
     xlab = "Cohort", ylab = " ",
     yaxt = "n", xaxt = "n", main = "(b) SC Curve", cex.main = 1.5)
axis(side = 2, at = seq(from = 0.4, to = 1, by = 0.1), las = 1, cex.axis = 0.8)
axis(side = 1, at = seq(from = 1885, to = 2010, by = 15), las = 0, cex.axis = 0.8)
lines(x = new_data$c, y = new_data$lo, col = adjustcolor("gray", alpha.f = 0.7))
lines(x = new_data$c, y = new_data$hi, col = adjustcolor("gray", alpha.f = 0.7))

## PANEL (c): THE PERIOD NONLINEARITIES (SHOCKS)

# intercept plus the period nonlinearities, one row per period
p_linear <- H$X[, i.a] + H$X[, i.c]
new_data <- cbind(p.index = p.index, component_curve(c(i.mu, i.sp), sort_by = p_linear))
stopifnot(nrow(new_data) == length(p.index))

plot(x = new_data$p.index, y = new_data$pred, type = "l", ylim = c(0.5, 1),
     xlab = "Period", ylab = "Period Fluctuations",
     yaxt = "n", xaxt = "n", main = "(c) Period Fluctuations", cex.main = 1.5)
axis(side = 2, at = seq(from = 0.5, to = 1, by = 0.1), las = 1, cex.axis = 0.8)
axis(side = 1, at = seq(from = 1970, to = 2020, by = 5), las = 0, cex.axis = 0.8)

# reference line at the mean of the period fluctuations
abline(a = mean(new_data$pred), b = 0, lty = "dashed")

lines(x = new_data$p.index, y = new_data$lo, col = adjustcolor("gray", alpha.f = 0.7))
lines(x = new_data$p.index, y = new_data$hi, col = adjustcolor("gray", alpha.f = 0.7))

# recession shading (recessions with a >= 1% peak-to-trough decline in GDP)
set1_red <- brewer.pal(8, "Set1")[1]

# 1980, 1981-1982
rect(xleft = 1980, ybottom = par("usr")[3], xright = 1982, ytop = par("usr")[4],
     border = NA, col = adjustcolor(set1_red, alpha.f = 0.1))

# 1990-1991
rect(xleft = 1990, ybottom = par("usr")[3], xright = 1991, ytop = par("usr")[4],
     border = NA, col = adjustcolor(set1_red, alpha.f = 0.1))

# 2007-2009
rect(xleft = 2007, ybottom = par("usr")[3], xright = 2009, ytop = par("usr")[4],
     border = NA, col = adjustcolor(set1_red, alpha.f = 0.1))

# 2020
rect(xleft = 2020, ybottom = par("usr")[3], xright = 2020.5, ytop = par("usr")[4],
     border = NA, col = adjustcolor(set1_red, alpha.f = 0.1))

invisible(dev.off())

# ---- 19. Figure: TwoPanel_Slopes_Gam.pdf ---------------------------------------------- #

pdf(file.path("Figures", "Banks", "TwoPanel_Slopes_Gam.pdf"), width = 14, height = 7)
par(mfrow = c(1, 2))

## PANEL (a): THETA1 (intercept and age linear component only)

new_data <- cbind(a.index = a.index, component_curve(c(i.mu, i.a), sort_col = 2))  # sorted by "a"

plot(x = new_data$a, y = new_data$pred, type = "l", ylim = c(0.4, 1),
     xlab = "Age", ylab = "Pr(Confidence)", xlim = c(15, 95),
     yaxt = "n", xaxt = "n", main = "(a) LC Slope", cex.main = 1.4)
axis(side = 2, at = seq(from = 0, to = 10, by = 0.1), las = 1, cex.axis = 0.8)
axis(side = 1, at = seq(from = 15, to = 100, by = 15), las = 0, cex.axis = 0.8)
lines(x = new_data$a, y = new_data$lo, col = adjustcolor("gray", alpha.f = 0.7))
lines(x = new_data$a, y = new_data$hi, col = adjustcolor("gray", alpha.f = 0.7))

## PANEL (b): THETA2 (intercept and cohort linear component only)

new_data <- cbind(c.index = c.index, component_curve(c(i.mu, i.c), sort_col = 2))  # sorted by "c"

plot(x = new_data$c, y = new_data$pred, type = "l", ylim = c(0.4, 1),
     xlab = "Cohort", ylab = " ",
     yaxt = "n", xaxt = "n", main = "(b) SC Slope", cex.main = 1.4)
axis(side = 2, at = seq(from = 0.4, to = 1, by = 0.1), las = 1, cex.axis = 0.8)
axis(side = 1, at = seq(from = 1885, to = 2010, by = 15), las = 0, cex.axis = 0.8)
lines(x = new_data$c, y = new_data$lo, col = adjustcolor("gray", alpha.f = 0.7))
lines(x = new_data$c, y = new_data$hi, col = adjustcolor("gray", alpha.f = 0.7))

invisible(dev.off())

# ---- 20. Figure: ThreePanel_NonlinearitiesComparison.pdf ------------------------------ #

# AGE: expanding the 5-year categorical deviations to piecewise-constant steps
expanded_data <- data.frame()
for (row in 1:nrow(AgeDeviations)) {
  age_group <- AgeDeviations[row, "Age"]
  repeat_rows <- data.frame(Age = seq(age_group, age_group + 4.75, by = 0.25),
                            Deviation = AgeDeviations[row, "Deviation"],
                            SE = AgeDeviations[row, "SE"],
                            `CI Lo` = AgeDeviations[row, "CI Lo"],
                            `CI Hi` = AgeDeviations[row, "CI Hi"])
  expanded_data <- rbind(expanded_data, repeat_rows)
}
merged_age_filled <- expanded_data

# PERIOD: expanding the 5-year categorical deviations to piecewise-constant steps
expanded_period_data <- data.frame()
for (row in 1:nrow(PeriodDeviations)) {
  period_group <- PeriodDeviations[row, "Period"]
  repeat_rows <- data.frame(Period = seq(period_group, period_group + 4.75, by = 0.25),
                            Deviation = PeriodDeviations[row, "Deviation"],
                            SE = PeriodDeviations[row, "SE"],
                            `CI Lo` = PeriodDeviations[row, "CI Lo"],
                            `CI Hi` = PeriodDeviations[row, "CI Hi"])
  expanded_period_data <- rbind(expanded_period_data, repeat_rows)
}
merged_period_filled <- expanded_period_data

# COHORT: the cohort groups are based on overlapping ten-year intervals
CohortDeviations[, 1] <- CohortDeviations[, 1] + -2.5

expanded_cohort_data <- data.frame()
for (row in 1:nrow(CohortDeviations)) {
  cohort_group <- CohortDeviations[row, "Cohort"]
  repeat_rows <- data.frame(Cohort = seq(cohort_group, cohort_group + 4.9, by = 0.1),
                            Deviation = CohortDeviations[row, "Deviation"],
                            SE = CohortDeviations[row, "SE"],
                            `CI Lo` = CohortDeviations[row, "CI Lo"],
                            `CI Hi` = CohortDeviations[row, "CI Hi"])
  expanded_cohort_data <- rbind(expanded_cohort_data, repeat_rows)
}
merged_cohort_filled <- expanded_cohort_data

# the GAM nonlinearity of one dimension as a response-scale deviation
nonlin_dev <- function(cols, sort_col) {
  relevant_X <- cbind(sort_col, H$X[, cols, drop = FALSE])
  relevant_X <- unique(relevant_X)
  relevant_X <- relevant_X[order(relevant_X[, 1]), ]
  relevant_X <- relevant_X[, -1]
  linear_predictor <- relevant_X %*% coef(m.lcsc)[cols]
  as.vector(plogis(H$intercept + linear_predictor) - plogis(H$intercept))
}

pdf(file.path("Figures", "Banks", "ThreePanel_NonlinearitiesComparison.pdf"), width = 13, height = 5)
par(mfrow = c(1, 3))

### PANEL (a): AGE NONLINEARITIES

# expanding margin for ylab
par(mar = c(5, 6, 4, 2) + 0.1)

new_data <- data.frame(a.index = a.index, pred = nonlin_dev(i.sa, H$X[, i.a]))

plot(x = new_data$a.index, y = new_data$pred - mean(new_data$pred), type = "l", ylim = c(-0.2, 0.15), col = "red", lwd = 1.2,
     xlab = "Age", ylab = "Age Nonlinearities", xlim = c(10, 95),
     yaxt = "n", xaxt = "n", main = "(a) Age Nonlinearities", cex.main = 1.5)
axis(side = 2, at = seq(from = -5, to = 5, by = 0.05), las = 1, cex.axis = 0.8)
axis(side = 1, at = seq(from = 10, to = 100, by = 15), las = 0, cex.axis = 0.8)

abline(a = 0, b = 0, lty = "dashed", col = "gray")

# categorical deviations in black (centered)
lines(x = merged_age_filled[, 1], y = merged_age_filled[, 2] - mean(merged_age_filled[, 2]), col = "black", lty = "solid")

AgeDev <- new_data

### PANEL (b): PERIOD NONLINEARITIES

new_data <- data.frame(p.index = p.index, pred = nonlin_dev(i.sp, H$X[, i.a] + H$X[, i.c]))

plot(x = new_data$p.index, y = new_data$pred - mean(new_data$pred), type = "l", ylim = c(-0.2, 0.155),
     xlab = "Period", ylab = "Period Nonlinearities", xlim = c(1972, 2024), col = "red",
     yaxt = "n", xaxt = "n", main = "(b) Period Nonlinearities", cex.main = 1.5)
axis(side = 2, at = seq(from = -5, to = 5, by = 0.05), las = 1, cex.axis = 0.8)
axis(side = 1, at = seq(from = 1975, to = 2025, by = 5), las = 0, cex.axis = 0.8)

abline(a = 0, b = 0, lty = "dashed", col = "gray")

# categorical deviations in black (centered)
lines(x = merged_period_filled[, 1], y = merged_period_filled[, 2] - mean(merged_period_filled[, 2]), col = "black", lty = "solid")

PeriodDev <- new_data

### PANEL (c): COHORT NONLINEARITIES

new_data <- data.frame(c.index = c.index, pred = nonlin_dev(i.sc, H$X[, i.c]))

plot(x = new_data$c.index, y = new_data$pred - mean(new_data$pred), type = "l", ylim = c(-0.1, 0.2), xlim = c(1880, 2007),
     xlab = "Cohort", ylab = "Cohort Nonlinearities", col = "red",
     yaxt = "n", xaxt = "n", main = "(c) Cohort Nonlinearities", cex.main = 1.5)
axis(side = 2, at = seq(from = -5, to = 5, by = 0.05), las = 1, cex.axis = 0.8)
axis(side = 1, at = seq(from = 1885, to = 2010, by = 15), las = 0, cex.axis = 0.8)

abline(a = 0, b = 0, lty = "dashed", col = "gray")

# categorical deviations in black (centered)
lines(x = merged_cohort_filled[, 1], y = merged_cohort_filled[, 2] - mean(merged_cohort_filled[, 2]), col = "black", lty = "solid")

CohortDev <- new_data

invisible(dev.off())

# ---- 21. Figure: Pnon_StratifyCohort.pdf ---------------------------------------------- #

# period nonlinearities
p_linear <- H$X[, i.a] + H$X[, i.c]
relevant_X <- cbind(p_linear, H$X[, i.sp])
relevant_X <- unique(relevant_X)
relevant_X <- relevant_X[order(relevant_X[, 1]), ]
relevant_X <- relevant_X[, -1]
intercept <- H$intercept
linear_predictor <- intercept + relevant_X %*% coef(m.lcsc)[i.sp]
predicted_probabilities <- exp(linear_predictor) / (1 + exp(linear_predictor))

p_new_data <- data.frame(p.index = p.index, pred = predicted_probabilities)

# expanding to the full age-period grid and building the age-cohort matrix
newdata <- expand.grid(a.index = a.index, p.index = p.index); newdata$c.index <- newdata$p.index - newdata$a.index
merged_data <- merge(p_new_data, newdata, by.x = "p.index", by.y = "p.index")
sorted_merged_data <- merged_data[, c("a.index", "p.index", "c.index", "pred")]

matAP.Pnon <- mean_by_ap(sorted_merged_data)
matAC.Pnon <- mean_by_ac(sorted_merged_data)

pdf(file.path("Figures", "Banks", "Pnon_StratifyCohort.pdf"), width = 9.75, height = 7.75)
par(mfrow = c(2, 2))

# matrix of period-nonlinearity careers (rows become cohorts after transposing)
d <- matAC.Pnon
d <- t(d)

# rows for the selected cohorts
row_numbers <- c(
  which(rownames(d) == "1920"),
  which(rownames(d) == "1940"),
  which(rownames(d) == "1960"),
  which(rownames(d) == "1980")
)

# one panel per selected cohort
for (i in row_numbers) {

  # blank plot
  plot(x = a.index, y = a.index, col = "black",
       type = "n", ylim = c(0.55, 0.9), xlim = c(5, 95),
       xlab = "Age", ylab = "Period Nonlinearities", yaxt = "n", xaxt = "n",
       main = rownames(d)[i],
       cex.main = 0.9, cex.lab = 0.8)

  # axis ticks
  axis(side = 2, at = seq(from = -2, to = 2, by = 0.1), las = 1, cex.axis = 0.8)
  axis(side = 1, at = seq(from = 15, to = 90, by = 15), las = 0, cex.axis = 0.8)

  # horizontal reference line at the intercept ON THE PROBABILITY SCALE, plogis(mu)
  abline(h = plogis(intercept), lty = "dashed", col = "darkgray")

  # the cohort's period-nonlinearity career
  lines(x = a.index, d[i, ], col = "black",
        lty = 1, lwd = 1)

  # 1990-1991 recession at this cohort's ages
  rect(xleft = 1989.5 - c.index[i], ybottom = par("usr")[3], xright = 1991.5 - c.index[i], ytop = par("usr")[4],
       border = NA, col = adjustcolor(set1_red, alpha.f = 0.1))

  # 2007-2009 recession at this cohort's ages
  rect(xleft = 2006.5 - c.index[i], ybottom = par("usr")[3], xright = 2009.5 - c.index[i], ytop = par("usr")[4],
       border = NA, col = adjustcolor(set1_red, alpha.f = 0.1))

}

invisible(dev.off())

# ---- 22. Done ------------------------------------------------------------------------- #

figs <- c("CohortCareers.raw", "2D_matAP.raw", "2D_matAP.apc_lcsc", "CohortCareersTuftean.apc_lcsc",
          "Logit_TwoPanel", "ThreePanel_APC_Gam", "2D_matAP.eta", "TwoPanel_Slopes_Gam",
          "ThreePanel_NonlinearitiesComparison", "Pnon_StratifyCohort")
fig.files <- file.path("Figures", "Banks", paste0(figs, ".pdf"))
out.files <- file.path("Output", c("banks_table1_fit_statistics.csv", "banks_table2_model_summary.csv",
                                   "banks_concurvity_blocks.csv", "tex/banks_table1.tex", "tex/banks_table2.tex"))
ok <- file.exists(c(fig.files, out.files))
cat("\n==== Output files ====\n")
cat(sprintf("  [%s] %s\n", ifelse(ok, "ok", "MISSING"), c(fig.files, out.files)), sep = "")
if (!all(ok)) warning("Some output files are missing.")
cat("Done: 01_confidence_banks_gss.R (figures model: ", m.lcsc$method, "; tables model: ", m.lcsc.ML$method,
    if (LEGACY_UBRE) "; LEGACY_UBRE = TRUE" else "", ")\n", sep = "")

## END OF R CODE
