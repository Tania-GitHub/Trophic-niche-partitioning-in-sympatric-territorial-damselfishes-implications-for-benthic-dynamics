# =========================================================
# Prey accumulation curves + model fitting
# Final figure in English, with combined legend
# and adjusted text sizes
# =========================================================

# Clear environment
rm(list = ls())

# Load packages
need <- c(
  "readr", "dplyr", "purrr", "stringr", "ggplot2", "vegan",
  "ragg", "dfoptim", "tibble", "tidyr", "patchwork"
)

to_install <- setdiff(need, rownames(installed.packages()))

if (length(to_install)) {
  install.packages(to_install, dependencies = TRUE)
}

invisible(lapply(need, library, character.only = TRUE))

# General settings
set.seed(123)

nperm <- 9999
infile <- "ITEM.csv"
out_dir <- "accumulation_results"

if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

# ---------- Helper functions ----------
make_safe_id <- function(x) {
  x <- iconv(x, to = "ASCII//TRANSLIT")
  x <- stringr::str_replace_all(x, "[^A-Za-z0-9]+", "_")
  x <- stringr::str_replace_all(x, "^_+|_+$", "")
  tolower(x)
}

write_csv_safe <- function(x, path) {
  ok <- TRUE
  msg <- NULL
  
  tryCatch(
    readr::write_csv(x, path),
    error = function(e) {
      ok <<- FALSE
      msg <<- conditionMessage(e)
    }
  )
  
  if (!ok) {
    ts <- format(Sys.time(), "%Y%m%d_%H%M%S")
    alt <- file.path(
      dirname(path),
      paste0(
        tools::file_path_sans_ext(basename(path)),
        "_", ts, ".csv"
      )
    )
    
    readr::write_csv(x, alt)
    
    message(
      "Could not write '", path,
      "'. Reason: ", msg,
      "\nAlternative file saved at: ", alt
    )
  }
}

species_color <- function(sp) {
  if (sp == "S. acapulcoensis") "#C16540"
  else if (sp == "S. flavilatus") "#00C1C8"
  else "#444444"
}

species_fullname <- function(sp) {
  if (sp == "S. acapulcoensis") return("Stegastes acapulcoensis")
  if (sp == "S. flavilatus") return("Stegastes flavilatus")
  return(sp)
}

# ---------- Read and clean data ----------
df <- readr::read_csv(
  infile,
  show_col_types = FALSE
)

stopifnot(all(c("Species", "Stomach") %in% names(df)))

item_cols <- setdiff(names(df), c("Species", "Stomach"))

df[item_cols] <- lapply(df[item_cols], function(x) {
  x <- suppressWarnings(as.numeric(x))
  x[is.na(x)] <- 0
  x
})

# ---------- Models ----------
f_clench <- function(n, par) {
  a <- exp(par[1])
  b <- exp(par[2])
  (a * n) / (1 + b * n)
}

f_negexp <- function(n, par) {
  a <- exp(par[1])
  b <- exp(par[2])
  (a * (1 - exp(-b * n))) / b
}

slope_clench <- function(nmax, par) {
  a <- exp(par[1])
  b <- exp(par[2])
  a / (1 + b * nmax)^2
}

slope_negexp <- function(nmax, par) {
  a <- exp(par[1])
  b <- exp(par[2])
  a * exp(-b * nmax)
}

asympt_common <- function(par) {
  a <- exp(par[1])
  b <- exp(par[2])
  a / b
}

# ---------- Metrics ----------
get_aic <- function(sse, n, k) {
  sse <- max(sse, .Machine$double.eps)
  n * log(sse / n) + 2 * k
}

get_aicc_from_aic <- function(aic, n, k) {
  if (n - k - 1 <= 0) return(Inf)
  aic + (2 * k * (k + 1)) / (n - k - 1)
}

# ---------- Fitting wrappers ----------
obj_sse <- function(f_model, x_n, y_sobs) {
  function(p) {
    pred <- f_model(x_n, p)
    sum((y_sobs - pred)^2)
  }
}

fit_with_optim_method <- function(x_n, y_sobs, f_model, par0, method) {
  obj <- obj_sse(f_model, x_n, y_sobs)
  
  res <- optim(
    par = par0,
    fn = obj,
    method = method,
    control = list(maxit = 1e5)
  )
  
  pred <- f_model(x_n, res$par)
  sse <- sum((y_sobs - pred)^2)
  sst <- sum((y_sobs - mean(y_sobs))^2)
  r2 <- 1 - sse / sst
  
  list(
    par = res$par,
    pred = pred,
    sse = sse,
    r2 = r2,
    converged = (res$convergence == 0),
    method_label = method
  )
}

fit_with_hj <- function(x_n, y_sobs, f_model, par0) {
  obj <- obj_sse(f_model, x_n, y_sobs)
  
  res <- dfoptim::hjkb(
    par = par0,
    fn = obj,
    control = list(maxfeval = 1e5)
  )
  
  pred <- f_model(x_n, res$par)
  sse <- sum((y_sobs - pred)^2)
  sst <- sum((y_sobs - mean(y_sobs))^2)
  r2 <- 1 - sse / sst
  
  conv <- if (!is.null(res$convergence)) {
    res$convergence == 0
  } else {
    TRUE
  }
  
  list(
    par = res$par,
    pred = pred,
    sse = sse,
    r2 = r2,
    converged = conv,
    method_label = "Hooke-Jeeves"
  )
}

fit_sequential_nm_bfgs <- function(x_n, y_sobs, f_model, par0) {
  fm1 <- fit_with_optim_method(
    x_n, y_sobs, f_model,
    par0,
    method = "Nelder-Mead"
  )
  
  fm2 <- fit_with_optim_method(
    x_n, y_sobs, f_model,
    fm1$par,
    method = "BFGS"
  )
  
  fm2$method_label <- "Nelder-Mead -> BFGS"
  fm2$converged <- fm1$converged && fm2$converged
  fm2
}

fit_model_strategy <- function(x_n, y_sobs, f_model, par0, strategy) {
  if (strategy == "quasi-newton") {
    return(fit_with_optim_method(x_n, y_sobs, f_model, par0, "BFGS"))
  }
  
  if (strategy == "simplex") {
    return(fit_with_optim_method(x_n, y_sobs, f_model, par0, "Nelder-Mead"))
  }
  
  if (strategy == "simplex+quasi-newton") {
    return(fit_sequential_nm_bfgs(x_n, y_sobs, f_model, par0))
  }
  
  if (strategy == "hooke-jeeves") {
    return(fit_with_hj(x_n, y_sobs, f_model, par0))
  }
  
  stop("Unknown strategy")
}

pretty_strategy <- function(s) {
  switch(
    s,
    "quasi-newton" = "Quasi-Newton",
    "simplex" = "Simplex",
    "simplex+quasi-newton" = "Simplex & Quasi-Newton",
    "hooke-jeeves" = "Hooke–Jeeves",
    s
  )
}

# ---------- Accumulation curve ----------
accum_for_species <- function(df_sp, nperm = 9999) {
  mat <- as.matrix(
    dplyr::select(df_sp, dplyr::all_of(item_cols))
  )
  
  acc <- vegan::specaccum(
    mat,
    method = "random",
    permutations = nperm
  )
  
  tibble::tibble(
    n = seq_along(acc$richness),
    Sobs = as.numeric(acc$richness),
    sd = as.numeric(acc$sd),
    lwr = pmax(Sobs - 1.96 * sd, 0),
    upr = Sobs + 1.96 * sd
  )
}

# ---------- Pipeline by species ----------
run_all_for_species <- function(sp_name) {
  message("Processing: ", sp_name)
  
  sp_id <- make_safe_id(sp_name)
  sp_dir <- file.path(out_dir, sp_id)
  
  if (!dir.exists(sp_dir)) {
    dir.create(sp_dir, recursive = TRUE)
  }
  
  df_sp <- df %>%
    dplyr::filter(Species == sp_name) %>%
    dplyr::arrange(Stomach)
  
  stopifnot(nrow(df_sp) > 1)
  
  acc_df <- accum_for_species(df_sp, nperm = nperm)
  
  x_n <- acc_df$n
  y <- acc_df$Sobs
  nmax <- max(x_n)
  S_last <- tail(y, 1)
  
  sp_col <- species_color(sp_name)
  
  obs_label <- paste0("Observed – ", species_fullname(sp_name))
  
  legend_limits <- c(
    "Clench model",
    "Observed – Stegastes acapulcoensis",
    "Observed – Stegastes flavilatus"
  )
  
  legend_labels <- c(
    expression("Clench model"),
    expression("Observed \u2013 " * italic("Stegastes acapulcoensis")),
    expression("Observed \u2013 " * italic("Stegastes flavilatus"))
  )
  
  # Initial values
  k <- min(5, length(y))
  slope0 <- max(sum((1:k) * y[1:k]) / sum((1:k)^2), 1e-3)
  A0 <- max(S_last * 1.1, slope0 * 3)
  b0 <- max(slope0 / A0, 1e-4)
  a0 <- max(slope0, 1e-4)
  
  par0_clench <- log(c(a0, b0))
  par0_negexp <- log(c(a0, b0))
  
  strategies <- c(
    "quasi-newton",
    "simplex",
    "simplex+quasi-newton",
    "hooke-jeeves"
  )
  
  all_rows <- list()
  
  for (s in strategies) {
    
    fitC <- fit_model_strategy(x_n, y, f_clench, par0_clench, s)
    fitE <- fit_model_strategy(x_n, y, f_negexp, par0_negexp, s)
    
    models <- list(
      Clench = list(
        fit = fitC,
        f = f_clench,
        slope = function(p) slope_clench(nmax, p),
        asympt = asympt_common
      ),
      NegExp = list(
        fit = fitE,
        f = f_negexp,
        slope = function(p) slope_negexp(nmax, p),
        asympt = asympt_common
      )
    )
    
    res_tbl_s <- purrr::imap_dfr(models, function(m, name) {
      par <- m$fit$par
      pred <- m$f(x_n, par)
      sse <- sum((y - pred)^2)
      sst <- sum((y - mean(y))^2)
      r2 <- 1 - sse / sst
      aic <- get_aic(sse, length(y), k = 2)
      aicc <- get_aicc_from_aic(aic, length(y), k = 2)
      slope_end <- m$slope(par)
      asympt <- m$asympt(par)
      completeness <- S_last / asympt
      
      tibble::tibble(
        Species = sp_name,
        Species_id = sp_id,
        Strategy = pretty_strategy(s),
        Strategy_key = s,
        Model = name,
        Fitting_method = m$fit$method_label,
        Converged = m$fit$converged,
        par1_log = par[1],
        par2_log = par[2],
        a_par = exp(par[1]),
        b_par = exp(par[2]),
        SSE = sse,
        R2 = r2,
        AIC = aic,
        AICc = aicc,
        S_asymptote = asympt,
        Final_observed_S = S_last,
        Completeness = completeness,
        Final_slope = slope_end
      )
    })
    
    all_rows[[s]] <- res_tbl_s
  }
  
  results_all <- dplyr::bind_rows(all_rows)
  
  best_overall <- results_all %>%
    dplyr::slice_min(SSE, n = 1, with_ties = FALSE)
  
  model_funs <- list(
    Clench = f_clench,
    NegExp = f_negexp
  )
  
  best_fun <- model_funs[[best_overall$Model[[1]]]]
  best_par <- c(
    best_overall$par1_log[[1]],
    best_overall$par2_log[[1]]
  )
  
  curve_best_df <- tibble::tibble(
    n = x_n,
    S_hat = best_fun(x_n, best_par)
  )
  
  obs_line_df <- acc_df %>%
    dplyr::transmute(
      n,
      y = Sobs,
      Series = obs_label
    )
  
  fit_line_df <- curve_best_df %>%
    dplyr::transmute(
      n,
      y = S_hat,
      Series = "Clench model"
    )
  
  obs_line_df$Series <- factor(obs_line_df$Series, levels = legend_limits)
  fit_line_df$Series <- factor(fit_line_df$Series, levels = legend_limits)
  
  label_x <- 31.5
  label_y <- if (sp_name == "S. acapulcoensis") 98 else 97
  
  p <- ggplot2::ggplot(acc_df, aes(n, Sobs)) +
    ggplot2::geom_ribbon(
      aes(ymin = lwr, ymax = upr),
      fill = sp_col,
      alpha = 0.15,
      show.legend = FALSE
    ) +
    ggplot2::geom_point(
      size = 2,
      color = sp_col,
      show.legend = FALSE
    ) +
    ggplot2::geom_line(
      data = obs_line_df,
      aes(n, y, color = Series, linetype = Series),
      linewidth = 0.7,
      show.legend = TRUE
    ) +
    ggplot2::geom_line(
      data = fit_line_df,
      aes(n, y, color = Series, linetype = Series),
      linewidth = 1.0,
      show.legend = TRUE
    ) +
    ggplot2::scale_color_manual(
      name = NULL,
      limits = legend_limits,
      labels = legend_labels,
      drop = FALSE,
      values = c(
        "Clench model" = "black",
        "Observed – Stegastes acapulcoensis" = "#C16540",
        "Observed – Stegastes flavilatus" = "#00C1C8"
      )
    ) +
    ggplot2::scale_linetype_manual(
      name = NULL,
      limits = legend_limits,
      labels = legend_labels,
      drop = FALSE,
      values = c(
        "Clench model" = "dashed",
        "Observed – Stegastes acapulcoensis" = "solid",
        "Observed – Stegastes flavilatus" = "solid"
      )
    ) +
    ggplot2::labs(
      x = "Number of stomachs examined",
      y = "Accumulated prey"
    ) +
    ggplot2::theme_classic(base_size = 11) +
    ggplot2::theme(
      axis.title = element_text(size = 11, color = "black"),
      axis.text = element_text(size = 9, color = "black"),
      legend.position = "top",
      legend.title = element_blank(),
      legend.text = element_text(size = 11, color = "black"),
      legend.key.size = grid::unit(0.9, "cm")
    ) +
    ggplot2::guides(
      color = ggplot2::guide_legend(
        override.aes = list(
          linetype = c("dashed", "solid", "solid"),
          linewidth = c(1.1, 1.1, 1.1),
          color = c("black", "#C16540", "#00C1C8")
        ),
        order = 1
      ),
      linetype = "none"
    ) +
    ggplot2::annotate(
      "text",
      x = label_x,
      y = label_y,
      label = if (sp_name == "S. acapulcoensis") "A" else "B",
      size = 4,
      fontface = "bold"
    )
  
  png_path <- file.path(sp_dir, paste0("curve_", sp_id, "_best.png"))
  pdf_path <- file.path(sp_dir, paste0("curve_", sp_id, "_best.pdf"))
  
  ggsave(
    png_path,
    p,
    width = 7,
    height = 5,
    dpi = 300,
    device = ragg::agg_png,
    bg = "white"
  )
  
  ggsave(
    pdf_path,
    p,
    width = 7,
    height = 5,
    device = cairo_pdf,
    bg = "white"
  )
  
  write_csv_safe(
    results_all,
    file.path(sp_dir, paste0("model_summary_", sp_id, ".csv"))
  )
  
  write_csv_safe(
    best_overall,
    file.path(sp_dir, paste0("best_global_", sp_id, ".csv"))
  )
  
  list(
    all = results_all,
    best_overall = best_overall,
    plot = p,
    dir = sp_dir
  )
}

# ---------- Run by species ----------
species <- unique(df$Species)

results <- purrr::map(species, run_all_for_species)

names(results) <- species

# ---------- Combine summaries ----------
summary_total <- dplyr::bind_rows(
  purrr::map(results, "all")
)

write_csv_safe(
  summary_total,
  file.path(out_dir, "model_summary_all_species.csv")
)

best_global <- dplyr::bind_rows(
  purrr::map(results, "best_overall")
)

write_csv_safe(
  best_global,
  file.path(out_dir, "best_global_all_species.csv")
)

# =========================================================
# Same scale in all figures
# Fixed y-axis from 15 to 100
# =========================================================

xmax_global <- max(
  purrr::map_dbl(
    results,
    ~ max(.x$plot$data$n, na.rm = TRUE)
  ),
  na.rm = TRUE
)

plots_no_title <- purrr::imap(
  results,
  ~ .x$plot +
    ggplot2::coord_cartesian(
      xlim = c(1, xmax_global),
      ylim = c(15, 100)
    )
)

plots_no_title <- plots_no_title[sort(names(plots_no_title))]

fig_2col <- patchwork::wrap_plots(
  plots_no_title,
  ncol = 2
) +
  patchwork::plot_layout(guides = "collect")

fig_2col <- fig_2col &
  ggplot2::theme(
    legend.position = "top",
    legend.text = element_text(size = 10, color = "black"),
    legend.key.size = grid::unit(0.9, "cm"),
    axis.text = element_text(size = 10, color = "black"),
    axis.title = element_text(size = 11, color = "black")
  )

print(fig_2col)

# ---------- Save combined figure ----------
ggplot2::ggsave(
  filename = "F:/Doctorado/TESIS version 3/Capitulo 2.- corregido/figuras/clench_model2.tiff",
  plot = fig_2col,
  device = ragg::agg_tiff,
  width = 18,
  height = 15,
  units = "cm",
  dpi = 500,
  compression = "lzw",
  bg = "white"
)

message("\nDone. Root folder: ", normalizePath(out_dir))