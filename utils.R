# utils.R
# Helper functions for MixCat model data construction and output extraction/plotting.
# All extract_* / plot_* functions accept an optional group_labels character vector
# (e.g. c("Rural", "Urban")) when the model was fitted with NGroup > 1.


#------------------------------------------------------------------------------
# Internal helpers
#------------------------------------------------------------------------------

.is_grouped <- function(model_data) model_data$NGroup > 1

.group_labels <- function(model_data, group_labels) {
  if (!is.null(group_labels)) group_labels else paste0("Group ", seq_len(model_data$NGroup))
}


#------------------------------------------------------------------------------
# Model data constructor
#------------------------------------------------------------------------------

# Build the Stan input list for any MixCat model.
# When group is NULL (default), a single stratum is assumed (NGroup = 1).
#
# Arguments:
#   titer        - numeric vector of antibody titers
#   age_group    - integer vector of age group indices (1:NAgeG) per individual
#   ageL         - numeric vector of age group lower bounds
#   ageU         - numeric vector of age group upper bounds
#   group    - integer vector of stratum indices (1:NGroup); NULL = unstratified (NGroup = 1)
#   prior_means  - length-4 vector: prior means for mu0, mu1, sd0, sd1
#   prior_sds    - length-4 vector: prior SDs for mu0, mu1, sd0, sd1
#   titer_pad    - padding below/above observed titer range for the fit grid
#   titer_step   - step size of the titer fit grid
make_model_data <- function(titer,
                            age_group,
                            ageL,
                            ageU,
                            group   = NULL,
                            prior_means = c(2, 3, 1, 1),
                            prior_sds   = c(1, 1, 1, 1),
                            titer_pad   = 1,
                            titer_step  = 0.2) {

  N     <- length(titer)
  NAgeG <- length(ageL)

  if (is.null(group)) group <- rep(1L, N)
  NGroup <- max(group)

  n_age <- matrix(0L, nrow = NGroup, ncol = NAgeG)
  for (g in seq_len(NGroup)) {
    n_age[g, ] <- as.integer(tabulate(age_group[group == g], nbins = NAgeG))
  }

  y_fit <- seq(min(titer) - titer_pad, max(titer) + titer_pad, by = titer_step)

  list(
    N           = N,
    y           = titer,
    ageG        = age_group,
    NAgeG       = NAgeG,
    ageL        = ageL,
    ageU        = ageU,
    NGroup      = NGroup,
    group        = group,
    n_age       = n_age,
    y_fit       = y_fit,
    predL       = length(y_fit),
    prior_means = prior_means,
    prior_sds   = prior_sds
  )
}


#------------------------------------------------------------------------------
# Seroprevalence by age — MixtureCatalytic
#------------------------------------------------------------------------------

extract_sero_age_cat <- function(draws, model_data, group_labels = NULL) {

  labels <- .group_labels(model_data, group_labels)
  ages   <- seq(min(model_data$ageL), max(model_data$ageU))

  curve_from_lambda <- function(lambda_draws) {
    foi  <- quantile(lambda_draws, c(0.5, 0.025, 0.975))
    prev <- t(1 - exp(-outer(foi, ages)))
    data.frame(age = ages, seroprev = prev[, 1], criL = prev[, 2], criU = prev[, 3])
  }

  result <- do.call(rbind, lapply(seq_len(model_data$NGroup), function(g) {
    df <- curve_from_lambda(draws$lambda[, g])
    df$group <- labels[g]
    df
  }))

  if (!.is_grouped(model_data)) {
    result$group <- NULL
    return(result)
  }

  overall_df       <- curve_from_lambda(draws$lambda_overall)
  overall_df$group <- "Overall"
  rbind(overall_df, result)
}

plot_sero_age_cat <- function(draws, model_data, group_labels = NULL) {

  predPrev <- extract_sero_age_cat(draws, model_data, group_labels)

  p <- ggplot(predPrev, aes(age, seroprev)) +
    theme_bw() +
    xlab("Age") + ylab("Estimated seroprevalence") +
    ylim(0, 1) +
    theme(text = element_text(size = 16))

  if (.is_grouped(model_data)) {
    p +
      geom_ribbon(aes(ymin = criL, ymax = criU, fill = group), alpha = 0.3, color = NA) +
      geom_line(aes(color = group), linewidth = 1) +
      labs(color = NULL, fill = NULL) +
      theme(legend.position = "top")
  } else {
    p +
      geom_ribbon(aes(ymin = criL, ymax = criU), fill = "indianred2", alpha = 0.4) +
      geom_line(color = "indianred2", linewidth = 1)
  }
}


#------------------------------------------------------------------------------
# Seroprevalence by age group — Mixture
#------------------------------------------------------------------------------

extract_sero_age_mix <- function(draws, model_data, group_labels = NULL) {

  labels     <- .group_labels(model_data, group_labels)
  age_labels <- paste0(model_data$ageL, "-", model_data$ageU)
  age_mids   <- (model_data$ageL + model_data$ageU) / 2

  result <- do.call(rbind, lapply(seq_len(model_data$NGroup), function(g) {
    sero_g <- draws$sero[, , g]
    data.frame(
      age_label = factor(age_labels, levels = age_labels),
      age_mid   = age_mids,
      seroprev  = apply(sero_g, 2, median),
      criL      = apply(sero_g, 2, quantile, 0.025),
      criU      = apply(sero_g, 2, quantile, 0.975),
      group     = labels[g]
    )
  }))

  if (!.is_grouped(model_data)) result$group <- NULL
  result
}

plot_sero_age_mix <- function(draws, model_data, group_labels = NULL) {

  predPrev <- extract_sero_age_mix(draws, model_data, group_labels)

  p <- ggplot(predPrev, aes(x = age_label, y = seroprev)) +
    theme_bw() +
    xlab("Age group") + ylab("Estimated seroprevalence") +
    ylim(0, 1) +
    theme(text = element_text(size = 16),
          axis.text.x = element_text(angle = 45, hjust = 1))

  if (.is_grouped(model_data)) {
    p +
      geom_linerange(aes(ymin = criL, ymax = criU, color = group),
                     position = position_dodge(width = 0.5), linewidth = 0.8) +
      geom_point(aes(color = group), size = 3,
                 position = position_dodge(width = 0.5)) +
      labs(color = NULL) +
      theme(legend.position = "top")
  } else {
    p +
      geom_linerange(aes(ymin = criL, ymax = criU), color = "indianred2",
                     linewidth = 0.8) +
      geom_point(color = "indianred2", size = 3)
  }
}


#------------------------------------------------------------------------------
# Overall and group-specific seroprevalence summary
#------------------------------------------------------------------------------

extract_sero <- function(draws, model_data, group_labels = NULL) {
  smry <- function(x) {
    q <- quantile(x, c(0.5, 0.025, 0.975))
    data.frame(median = q[[1]], criL = q[[2]], criU = q[[3]], row.names = NULL)
  }
  rows <- list(data.frame(label = "Overall", smry(draws$sero_all)))
  if (.is_grouped(model_data)) {
    labels <- .group_labels(model_data, group_labels)
    for (g in seq_len(model_data$NGroup))
      rows[[length(rows) + 1]] <- data.frame(label = labels[g], smry(draws$sero_grp[, g]))
  }
  do.call(rbind, rows)
}


#------------------------------------------------------------------------------
# Force of infection summary — MixtureCatalytic only
#------------------------------------------------------------------------------

extract_foi <- function(draws, model_data, group_labels = NULL) {
  labels <- .group_labels(model_data, group_labels)

  smry_lambda <- function(x) {
    q <- quantile(x, c(0.5, 0.025, 0.975))
    data.frame(median = q[[1]], criL = q[[2]], criU = q[[3]], row.names = NULL)
  }

  per_group <- do.call(rbind, lapply(seq_len(model_data$NGroup), function(g) {
    data.frame(group = labels[g], smry_lambda(draws$lambda[, g]))
  }))

  if (!.is_grouped(model_data)) {
    per_group$group <- "Overall"
    return(per_group)
  }

  overall_row <- data.frame(group = "Overall", smry_lambda(draws$lambda_overall))
  rbind(overall_row, per_group)
}


#------------------------------------------------------------------------------
# Overall distribution fit (seronegative / seropositive / overall)
#------------------------------------------------------------------------------

extract_dist_fit <- function(draws, model_data) {

  summarise_component <- function(mat, component) {
    data.frame(
      titer       = model_data$y_fit,
      density_med = apply(mat, 2, median),
      density_lo  = apply(mat, 2, quantile, 0.025),
      density_hi  = apply(mat, 2, quantile, 0.975),
      component   = component
    )
  }

  fit_df <- rbind(
    summarise_component(draws$fitNeg, "Seronegative"),
    summarise_component(draws$fitPos, "Seropositive"),
    summarise_component(draws$fitAll, "Overall")
  )
  fit_df$component <- factor(fit_df$component,
                             levels = c("Overall", "Seronegative", "Seropositive"))
  fit_df
}

plot_dist_fit <- function(draws, model_data) {

  fit_df <- extract_dist_fit(draws, model_data)
  pal    <- c("Overall" = "navy", "Seronegative" = "dodgerblue", "Seropositive" = "indianred2")

  ggplot() +
    geom_histogram(data = data.frame(titer = model_data$y),
                   aes(x = titer, y = after_stat(density)),
                   fill = "grey70", alpha = 0.6, bins = 25, color = "white") +
    geom_ribbon(data = fit_df,
                aes(x = titer, ymin = density_lo, ymax = density_hi, fill = component),
                alpha = 0.2) +
    geom_line(data = fit_df,
              aes(x = titer, y = density_med, color = component),
              linewidth = 1) +
    scale_color_manual(values = pal, name = "Component") +
    scale_fill_manual(values = pal, name = "Component") +
    theme_bw() +
    xlab("Titer") + ylab("Density") +
    theme(text = element_text(size = 14), legend.position = "top")
}


#------------------------------------------------------------------------------
# Observed vs predicted mean titer by age group
#------------------------------------------------------------------------------

extract_mean_titer <- function(draws, model_data, group_labels = NULL) {

  labels     <- .group_labels(model_data, group_labels)
  age_labels <- paste0(model_data$ageL, "-", model_data$ageU)
  y_rep_mat  <- draws$y_rep

  result <- do.call(rbind, lapply(seq_len(model_data$NGroup), function(g) {

    g_idx  <- which(model_data$group == g)
    g_y    <- model_data$y[g_idx]
    g_ageG <- model_data$ageG[g_idx]

    obs_df <- data.frame(
      age_label  = factor(age_labels, levels = age_labels),
      mean_titer = tapply(g_y, g_ageG, mean),
      se         = tapply(g_y, g_ageG, function(x) sd(x) / sqrt(length(x)))
    )
    obs_df$lo <- obs_df$mean_titer - 1.96 * obs_df$se
    obs_df$hi <- obs_df$mean_titer + 1.96 * obs_df$se

    pred_mean <- do.call(rbind, lapply(seq_along(model_data$ageL), function(a) {
      idx        <- which(model_data$group == g & model_data$ageG == a)
      draw_means <- rowMeans(y_rep_mat[, idx, drop = FALSE])
      data.frame(
        age_label = factor(age_labels[a], levels = age_labels),
        mean_med  = median(draw_means),
        mean_lo   = quantile(draw_means, 0.025),
        mean_hi   = quantile(draw_means, 0.975)
      )
    }))

    df <- rbind(
      data.frame(age_label = obs_df$age_label,
                 mean = obs_df$mean_titer, lo = obs_df$lo, hi = obs_df$hi,
                 type = "Observed"),
      data.frame(age_label = pred_mean$age_label,
                 mean = pred_mean$mean_med, lo = pred_mean$mean_lo, hi = pred_mean$mean_hi,
                 type = "Predicted")
    )
    if (.is_grouped(model_data)) df$group <- labels[g]
    df
  }))

  result
}

plot_mean_titer <- function(draws, model_data, group_labels = NULL) {

  plot_df <- extract_mean_titer(draws, model_data, group_labels)

  p <- ggplot(plot_df, aes(x = age_label, y = mean, color = type, group = type)) +
    geom_linerange(aes(ymin = lo, ymax = hi),
                   position = position_dodge(width = 0.5)) +
    geom_point(size = 2, position = position_dodge(width = 0.5)) +
    scale_color_manual(values = c("Observed" = "grey50", "Predicted" = "seagreen3"),
                       name = NULL) +
    theme_bw() +
    xlab("Age group") + ylab("Mean titer") +
    theme(text = element_text(size = 16), legend.position = "top",
          axis.text.x = element_text(angle = 45, hjust = 1))

  if (.is_grouped(model_data)) p + facet_wrap(~group) else p
}



#------------------------------------------------------------------------------
# Probability of being seropositive vs observed titer
#------------------------------------------------------------------------------

extract_prob_seropos <- function(draws, model_data, group_labels = NULL) {

  labels  <- .group_labels(model_data, group_labels)
  pcp     <- apply(draws$pC, c(1, 3), function(x) exp(x[2]) / sum(exp(x[1:2])))
  result  <- data.frame(
    obs      = model_data$y,
    ppos_med = apply(pcp, 2, median),
    ppos_lo  = apply(pcp, 2, quantile, 0.025),
    ppos_hi  = apply(pcp, 2, quantile, 0.975)
  )
  if (.is_grouped(model_data)) result$group <- labels[model_data$group]
  result
}

plot_prob_seropos <- function(draws, model_data, group_labels = NULL) {

  plot_df <- extract_prob_seropos(draws, model_data, group_labels)

  p <- ggplot(plot_df, aes(obs, ppos_med)) +
    theme_bw() +
    xlab("Observed titer") + ylab("Probability of being seropositive") +
    theme(text = element_text(size = 16))

  if (.is_grouped(model_data)) {
    p +
      geom_point(aes(color = group), alpha = 0.5) +
      labs(color = NULL) +
      theme(legend.position = "top")
  } else {
    p + geom_point(alpha = 0.5)
  }
}
