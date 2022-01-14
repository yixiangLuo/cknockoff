# plot the important quantities in the calibration inequality
draw_ineq <- function(ineq_invest){
  p_ineq <- ineq_invest$p_ineq %>%
    mutate(rej_type = as.character(2*kn_rej + nokn_rej)) %>%
    gather(side, value, ineq_L, ineq_R) %>%
    mutate(color = as.character(2 * (side == "ineq_R") + !kept))

  p_ineq %>% ggplot() +
    geom_point(aes(x = pval, y = value, color = color, shape = rej_type)) +
    # ylim(0, 0.01) +
    # xlim(0, 0.01) +
    labs(x = "right-sided p-value of j", y = "value of DP_j or b_j",
         shape = "Shape indicates: numerator in DP_j is") +
    scale_shape_manual(labels = c("0" = "0", "1" = "1 due to Tj",
                                  "2" = "1 due to kn", "3" = "1 due to both"),
                       values = c("0" = 20, "1" = 4, "2" = 1, "3" = 13)) +
    scale_color_manual(labels = c("0" = "DP_j (sample in our sampling region)", "1" = "DP_j (not in)",
                                  "2" = "b_j (sample in our sampling region)", "3" = "b_j (not in)"),
                       values = c("0" = "dodgerblue1", "1" = "dodgerblue4",
                                  "2" = "firebrick1", "3" = "firebrick4"))
}

# plot the knockoff feature statistics
draw_kn_stat <- function(ineq_invest, invest_j, SRL = F){
  eskn_rej_mc <- data.frame(ineq_invest$eskn_rej_mc) %>%
    gather(hypo, kn_rej, starts_with("X"))

  x <- ineq_invest$kn_stat_nodes$vjy_nodes
  y1 <- ineq_invest$kn_stat_nodes$kn_abs_stat_j
  y2 <- ineq_invest$kn_stat_nodes$kn_stat_thrs
  bandwidth <- max(abs(diff(x)))
  x_range <- c(min(x), max(x))
  F1 <- locpoly(x, y1, bandwidth=bandwidth, degree=1, gridsize = 100, range.x = x_range)
  F2 <- locpoly(x, y2, bandwidth=bandwidth*2, degree=1, gridsize = 100, range.x = x_range)

  data.frame(ineq_invest$kn_stats_mc) %>%
    mutate(pval_j = ineq_invest$p_ineq$pval,
           yvj = ineq_invest$vjy_mc) %>%
    gather(hypo, kn_stat, starts_with("X")) %>%
    mutate(this_j = (hypo == paste0("X", invest_j)),
           eskn_rej = eskn_rej_mc$kn_rej) %>%
    ggplot(aes(x = yvj , y = abs(kn_stat), group = hypo, color = this_j)) +
    geom_line(alpha = 0.5, aes(size = this_j), show.legend = FALSE) +
    scale_size_manual(values = c(0, 1), breaks = c(F, T)) +
    geom_point(aes(alpha = eskn_rej)) +
    scale_alpha_manual(values = c(0, 1), breaks = c(F, T), guide="none") +
    geom_line(data = data.frame(vjy = F1$x,
                                value = (F1$y),
                                hypo = invest_j,
                                this_j = T), aes(x = vjy, y = value), color = "black") +
    geom_point(data = data.frame(vjy = x,
                                 value = (y1),
                                 hypo = invest_j,
                                 this_j = T), aes(x = vjy, y = value), color = "black") +
    geom_line(data = data.frame(vjy = F2$x,
                                value = (F2$y),
                                hypo = invest_j,
                                this_j = T), aes(x = vjy, y = value), color = "blue") +
    geom_point(data = data.frame(vjy = x,
                                 value = (y2),
                                 hypo = invest_j,
                                 this_j = T), aes(x = vjy, y = value), color = "blue") +
    geom_line(data = data.frame(vjy = ineq_invest$vjy_mc,
                                kn_stat_thres = (ineq_invest$kn_stat_thres),
                                hypo = invest_j,
                                this_j = T), aes(x = vjy, y = kn_stat_thres), color = "blue") +
    labs(y = "|W_i|", color = "i = invest_j?")
}

# draw the CDF of the right-sided p-value from the t-statistics for feature j
draw_pval_cdf <- function(pval_mc, weights_mc, sample_region, df){
  x_seq <- seq(from = 0, to = 1, length.out = 500)

  emp_cdf <- sapply(x_seq, function(x){
    sum((pval_mc <= x) * weights_mc)
  }) / sum(weights_mc)

  region_left <- pt(sample_region$left, df = df, lower.tail = FALSE)
  region_right <- pt(sample_region$right, df = df, lower.tail = FALSE)
  ref_cdf <- sapply(x_seq, function(x){
    sum(x - region_left[region_left <= x]) - sum(x - region_right[region_right <= x])
  }) / sum(region_right - region_left)

  data <- data.frame(x = x_seq, CDF = emp_cdf, type = "empirical result") %>%
    rbind(data.frame(x = x_seq, CDF = ref_cdf, type = "theoretical prediction"))


  ggplot(data) +
    geom_line(aes(x = x, y = CDF, color = type)) +
    labs(title = "CDF of right-sided p-value of H_j")
}

# P(a + b^T z <= 0 | \|z\|) for z ~ n(0, I_{df + 1} \sigma^2). See 11.summary Lemma 5
prob_sphere_onto_dir <- function(df, a, b_norm2, z_norm2){
  return(pt(-a * sqrt(df) / sqrt(pmax(b_norm2 * z_norm2 - a^2, 0)), df = df))
}

# draw the CDF of the projection of the sample z onto a unit vector in the
# residue space of X
draw_yres_cdf <- function(y_Pi_Xk1_mc, weights_mc, z_norm2, df){
  bound <- min(sqrt(z_norm2), max(abs(y_Pi_Xk1_mc)) * 1.1)
  x_seq <- seq(from = -bound, to = bound, length.out = 500)

  emp_cdf <- sapply(x_seq, function(x){
    sum((y_Pi_Xk1_mc <= x) * weights_mc)
  }) / sum(weights_mc)

  ref_cdf <- sapply(x_seq, function(x){
    prob_sphere_onto_dir(df, -x, 1, z_norm2)
  })

  data <- data.frame(x = x_seq, CDF = emp_cdf, type = "empirical result") %>%
    rbind(data.frame(x = x_seq, CDF = ref_cdf, type = "theoretical prediction"))


  ggplot(data) +
    geom_line(aes(x = x, y = CDF, color = type)) +
    labs(title = "CDF of z projected on to a unit vector in the residue space of X")
}
