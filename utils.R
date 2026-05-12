# utils.R
# Helper functions for extracting and plotting MixCat model outputs.
# All functions take `draws` (output of rstan::extract(fit)) as their first argument
# and `model_data` (the list passed to stan()) as their second argument.
#
# Each section has an extract_* function returning a data frame and a
# plot_* function that calls it.


#------------------------------------------------------------------------------
# Model data constructors
#------------------------------------------------------------------------------

# Build the model data list for any MixCat stan model.
# When foi_group is NULL (default), a single stratum is assumed (NGroup = 1).
#
# Arguments:
#   titer        - numeric vector of antibody titers
#   age_group    - integer vector of age group indices (1:NAgeG) for each individual
#   ageL         - numeric vector of age group lower bounds
#   ageU         - numeric vector of age group upper bounds
#   foi_group    - integer vector of stratum indices (1:NGroup); NULL = unstratified
#   prior_means  - length-4 numeric vector: prior means for mu0, mu1, sd0, sd1
#   prior_sds    - length-4 numeric vector: prior SDs for mu0, mu1, sd0, sd1
#   titer_pad    - padding added below/above observed titer range for the fit grid
#   titer_step   - step size of the titer fit grid
make_model_data <- function(titer,
                            age_group,
                            ageL,
                            ageU,
                            foi_group   = NULL,
                            prior_means = c(2, 3, 1, 1),
                            prior_sds   = c(1, 1, 1, 1),
                            titer_pad   = 1,
                            titer_step  = 0.2) {

  N     <- length(titer)
  NAgeG <- length(ageL)

  if (is.null(foi_group)) foi_group <- rep(1L, N)
  NGroup <- max(foi_group)

  n_age <- matrix(0L, nrow = NGroup, ncol = NAgeG)
  for (g in seq_len(NGroup)) {
    n_age[g, ] <- as.integer(tabulate(age_group[foi_group == g], nbins = NAgeG))
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
    foiG        = foi_group,
    n_age       = n_age,
    y_fit       = y_fit,
    predL       = length(y_fit),
    prior_means = prior_means,
    prior_sds   = prior_sds
  )
}


#------------------------------------------------------------------------------
# Seroprevalence by age
#------------------------------------------------------------------------------

extract_seroprev <- function(draws, model_data) {

  predPrev <- data.frame(age      = seq(min(model_data$ageL), max(model_data$ageU)),
                         seroprev = NA,
                         criL     = NA,
                         criU     = NA)

  foi <- quantile(draws$lambda[, 1], c(0.5, 0.025, 0.975))
  predPrev[, 2:4] <- t(1 - exp(-outer(foi, predPrev$age)))

  predPrev
}

plot_seroprev <- function(draws, model_data) {

  predPrev <- extract_seroprev(draws, model_data)

  ggplot(predPrev, aes(age, seroprev)) +
    geom_ribbon(aes(ymin = criL, ymax = criU), fill = "indianred2", alpha = 0.4) +
    geom_line(col = "indianred2", linewidth = 1) +
    theme_bw() +
    xlab("Age") + ylab("Estimated seroprevalence") +
    ylim(0, 1) +
    theme(text = element_text(size = 16))
}


# Mixture model version: sero is estimated directly per age group
extract_seroprev_mix <- function(draws, model_data) {

  age_labels <- paste0(model_data$ageL, "-", model_data$ageU)
  age_mids   <- (model_data$ageL + model_data$ageU) / 2

  sero_g1 <- draws$sero[, , 1]
  data.frame(
    age_label = factor(age_labels, levels = age_labels),
    age_mid   = age_mids,
    seroprev  = apply(sero_g1, 2, median),
    criL      = apply(sero_g1, 2, quantile, 0.025),
    criU      = apply(sero_g1, 2, quantile, 0.975)
  )
}

plot_seroprev_mix <- function(draws, model_data) {

  predPrev <- extract_seroprev_mix(draws, model_data)

  ggplot(predPrev, aes(x = age_label, y = seroprev)) +
    geom_linerange(aes(ymin = criL, ymax = criU), color = "indianred2",
                   linewidth = 0.8) +
    geom_point(color = "indianred2", size = 3) +
    theme_bw() +
    xlab("Age group") + ylab("Estimated seroprevalence") +
    ylim(0, 1) +
    theme(text = element_text(size = 16),
          axis.text.x = element_text(angle = 45, hjust = 1))
}


# Mixture model grouped version: seroprevalence per age group, faceted by stratum
extract_seroprev_mix_grouped <- function(draws, model_data, group_labels = NULL) {

  labels     <- .group_labels(model_data, group_labels)
  age_labels <- paste0(model_data$ageL, "-", model_data$ageU)
  age_mids   <- (model_data$ageL + model_data$ageU) / 2

  do.call(rbind, lapply(seq_len(model_data$NGroup), function(g) {
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
}

plot_seroprev_mix_grouped <- function(draws, model_data, group_labels = NULL) {

  predPrev <- extract_seroprev_mix_grouped(draws, model_data, group_labels)

  ggplot(predPrev, aes(x = age_label, y = seroprev, color = group)) +
    geom_linerange(aes(ymin = criL, ymax = criU),
                   position = position_dodge(width = 0.5), linewidth = 0.8) +
    geom_point(size = 3, position = position_dodge(width = 0.5)) +
    theme_bw() +
    xlab("Age group") + ylab("Estimated seroprevalence") +
    ylim(0, 1) +
    labs(color = NULL) +
    theme(text = element_text(size = 16), legend.position = "top",
          axis.text.x = element_text(angle = 45, hjust = 1))
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

extract_mean_titer <- function(draws, model_data) {

  age_labels <- paste0(model_data$ageL, "-", model_data$ageU)
  y_rep_mat  <- draws$y_rep

  obs_df <- data.frame(
    age_label  = factor(age_labels, levels = age_labels),
    mean_titer = tapply(model_data$y, model_data$ageG, mean),
    se         = tapply(model_data$y, model_data$ageG, function(x) sd(x) / sqrt(length(x)))
  )
  obs_df$lo <- obs_df$mean_titer - 1.96 * obs_df$se
  obs_df$hi <- obs_df$mean_titer + 1.96 * obs_df$se

  pred_mean <- do.call(rbind, lapply(seq_along(model_data$ageL), function(a) {
    idx        <- which(model_data$ageG == a)
    draw_means <- rowMeans(y_rep_mat[, idx, drop = FALSE])
    data.frame(
      age_label = factor(age_labels[a], levels = age_labels),
      mean_med  = median(draw_means),
      mean_lo   = quantile(draw_means, 0.025),
      mean_hi   = quantile(draw_means, 0.975)
    )
  }))

  rbind(
    data.frame(age_label = obs_df$age_label,
               mean      = obs_df$mean_titer,
               lo        = obs_df$lo,
               hi        = obs_df$hi,
               type      = "Observed"),
    data.frame(age_label = pred_mean$age_label,
               mean      = pred_mean$mean_med,
               lo        = pred_mean$mean_lo,
               hi        = pred_mean$mean_hi,
               type      = "Predicted")
  )
}

plot_mean_titer <- function(draws, model_data) {

  plot_df <- extract_mean_titer(draws, model_data)

  ggplot(plot_df, aes(x = age_label, y = mean, color = type, group = type)) +
    geom_linerange(aes(ymin = lo, ymax = hi),
                   position = position_dodge(width = 0.5)) +
    geom_point(size = 2, position = position_dodge(width = 0.5)) +
    scale_color_manual(values = c("Observed" = "grey50", "Predicted" = "seagreen3"),
                       name = NULL) +
    theme_bw() +
    xlab("Age group") + ylab("Mean titer") +
    theme(text = element_text(size = 16), legend.position = "top",
          axis.text.x = element_text(angle = 45, hjust = 1))
}


#------------------------------------------------------------------------------
# Individual observed vs predicted titer
#------------------------------------------------------------------------------

extract_titer_pred <- function(draws, model_data) {
  data.frame(
    obs  = model_data$y,
    pred = apply(draws$y_rep, 2, median)
  )
}

plot_titer_pred <- function(draws, model_data) {

  plot_df <- extract_titer_pred(draws, model_data)

  ggplot(plot_df, aes(obs, pred)) +
    geom_point(alpha = 0.5) +
    geom_abline(slope = 1, intercept = 0, linewidth = 0.8) +
    theme_bw() +
    xlab("Observed titer") + ylab("Predicted titer") +
    theme(text = element_text(size = 16))
}


#------------------------------------------------------------------------------
# Probability of being seropositive vs observed titer
#------------------------------------------------------------------------------

extract_prob_seropos <- function(draws, model_data) {
  # pC is draws x 2 x N; compute P(seropos) = softmax of component 2
  pcp     <- apply(draws$pC, c(1, 3), function(x) exp(x[2]) / sum(exp(x[1:2])))
  data.frame(
    obs       = model_data$y,
    ppos_med  = apply(pcp, 2, median),
    ppos_lo   = apply(pcp, 2, quantile, 0.025),
    ppos_hi   = apply(pcp, 2, quantile, 0.975)
  )
}

plot_prob_seropos <- function(draws, model_data) {

  plot_df <- extract_prob_seropos(draws, model_data)

  ggplot(plot_df, aes(obs, ppos_med)) +
    geom_point(alpha = 0.5) +
    theme_bw() +
    xlab("Observed titer") + ylab("Probability of being seropositive") +
    theme(text = element_text(size = 16))
}


#------------------------------------------------------------------------------
# GROUPED MODEL (MixtureCatalyticGrouped)
# For models with foi estimated per group (e.g. urban/rural).
# All functions accept an optional group_labels character vector for readable
# group names (e.g. c("Rural", "Urban")); defaults to "Group 1", "Group 2" etc.
#
# Key array dimensions from rstan::extract:
#   draws$lambda:  n_draws x NGroup
#   draws$sero:    n_draws x NAgeG x NGroup
#   draws$fitNeg/fitPos/fitAll: n_draws x predL x NGroup
#   draws$y_rep, draws$pC: unchanged (per individual)
#------------------------------------------------------------------------------

.group_labels <- function(model_data, group_labels) {
  if (!is.null(group_labels)) group_labels else paste0("Group ", seq_len(model_data$NGroup))
}


# Seroprevalence curves by age, one per foi group
extract_seroprev_grouped <- function(draws, model_data, group_labels = NULL) {

  labels <- .group_labels(model_data, group_labels)

  do.call(rbind, lapply(seq_len(model_data$NGroup), function(g) {
    ages <- seq(min(model_data$ageL), max(model_data$ageU))
    foi  <- quantile(draws$lambda[, g], c(0.5, 0.025, 0.975))
    prev <- t(1 - exp(-outer(foi, ages)))
    data.frame(age      = ages,
               seroprev = prev[, 1],
               criL     = prev[, 2],
               criU     = prev[, 3],
               group    = labels[g])
  }))
}

plot_seroprev_grouped <- function(draws, model_data, group_labels = NULL) {

  predPrev <- extract_seroprev_grouped(draws, model_data, group_labels)

  ggplot(predPrev, aes(age, seroprev, color = group, fill = group)) +
    geom_ribbon(aes(ymin = criL, ymax = criU), alpha = 0.3, color = NA) +
    geom_line(linewidth = 1) +
    theme_bw() +
    xlab("Age") + ylab("Estimated seroprevalence") +
    ylim(0, 1) +
    labs(color = NULL, fill = NULL) +
    theme(text = element_text(size = 16), legend.position = "top")
}


# Observed vs predicted mean titer by age group, faceted by foi group
extract_mean_titer_grouped <- function(draws, model_data, group_labels = NULL) {

  labels     <- .group_labels(model_data, group_labels)
  age_labels <- paste0(model_data$ageL, "-", model_data$ageU)
  y_rep_mat  <- draws$y_rep

  do.call(rbind, lapply(seq_len(model_data$NGroup), function(g) {

    g_idx  <- which(model_data$foiG == g)
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
      idx        <- which(model_data$foiG == g & model_data$ageG == a)
      draw_means <- rowMeans(y_rep_mat[, idx, drop = FALSE])
      data.frame(
        age_label = factor(age_labels[a], levels = age_labels),
        mean_med  = median(draw_means),
        mean_lo   = quantile(draw_means, 0.025),
        mean_hi   = quantile(draw_means, 0.975)
      )
    }))

    rbind(
      data.frame(age_label = obs_df$age_label, mean = obs_df$mean_titer,
                 lo = obs_df$lo, hi = obs_df$hi, type = "Observed",  group = labels[g]),
      data.frame(age_label = pred_mean$age_label, mean = pred_mean$mean_med,
                 lo = pred_mean$mean_lo, hi = pred_mean$mean_hi, type = "Predicted", group = labels[g])
    )
  }))
}

plot_mean_titer_grouped <- function(draws, model_data, group_labels = NULL) {

  plot_df <- extract_mean_titer_grouped(draws, model_data, group_labels)

  ggplot(plot_df, aes(x = age_label, y = mean, color = type, group = type)) +
    geom_linerange(aes(ymin = lo, ymax = hi),
                   position = position_dodge(width = 0.5)) +
    geom_point(size = 2, position = position_dodge(width = 0.5)) +
    scale_color_manual(values = c("Observed" = "grey50", "Predicted" = "seagreen3"),
                       name = NULL) +
    facet_wrap(~group) +
    theme_bw() +
    xlab("Age group") + ylab("Mean titer") +
    theme(text = element_text(size = 16), legend.position = "top",
          axis.text.x = element_text(angle = 45, hjust = 1))
}


# Individual obs vs pred titer, colored by foi group
extract_titer_pred_grouped <- function(draws, model_data, group_labels = NULL) {
  labels <- .group_labels(model_data, group_labels)
  data.frame(
    obs   = model_data$y,
    pred  = apply(draws$y_rep, 2, median),
    group = labels[model_data$foiG]
  )
}

plot_titer_pred_grouped <- function(draws, model_data, group_labels = NULL) {

  plot_df <- extract_titer_pred_grouped(draws, model_data, group_labels)

  ggplot(plot_df, aes(obs, pred, color = group)) +
    geom_point(alpha = 0.5) +
    geom_abline(slope = 1, intercept = 0, linewidth = 0.8) +
    labs(color = NULL) +
    theme_bw() +
    xlab("Observed titer") + ylab("Predicted titer") +
    theme(text = element_text(size = 16), legend.position = "top")
}


# P(seropositive) vs observed titer, colored by foi group
extract_prob_seropos_grouped <- function(draws, model_data, group_labels = NULL) {
  labels  <- .group_labels(model_data, group_labels)
  pcp     <- apply(draws$pC, c(1, 3), function(x) exp(x[2]) / sum(exp(x[1:2])))
  data.frame(
    obs      = model_data$y,
    ppos_med = apply(pcp, 2, median),
    ppos_lo  = apply(pcp, 2, quantile, 0.025),
    ppos_hi  = apply(pcp, 2, quantile, 0.975),
    group    = labels[model_data$foiG]
  )
}

plot_prob_seropos_grouped <- function(draws, model_data, group_labels = NULL) {

  plot_df <- extract_prob_seropos_grouped(draws, model_data, group_labels)

  ggplot(plot_df, aes(obs, ppos_med, color = group)) +
    geom_point(alpha = 0.5) +
    labs(color = NULL) +
    theme_bw() +
    xlab("Observed titer") + ylab("Probability of being seropositive") +
    theme(text = element_text(size = 16), legend.position = "top")
}
