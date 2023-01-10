set.seed(1)
options(digits = 2)
if (!require("pacman")) install.packages("pacman")
p_load("tidyverse", "writexl", "cowplot", "ggdist")

# Compute OR given RR and risk in control group
or_from_rr <- function(risk_ratio, untreated_risk) {
    risk_ratio * (1 - untreated_risk) / (1 - risk_ratio * untreated_risk)
}
# Compute OR given number of events in each group
or_from_events <- function(
    events_treated,
    events_control,
    n_treated,
    n_control
) {
    odds_treated <- events_treated / (n_treated - events_treated)
    odds_control  <- events_control / (n_control - events_control)
    odds_treated / odds_control
}

# Compute SE of logOR given 95% CI
se_from_ci <- function(lower_limit, upper_limit) {
    (upper_limit - lower_limit) / 3.92
}
# Compute SE of logOR given CDF (ie, Pr(logOR < threshold))
se_from_cdf <- function(mean, expected_cdf, threshold = log(1)) {
    data.frame(mean, se = runif(10^6)) |>
    dplyr::mutate(cdf = round(pnorm(threshold, mean, se), 4)) |>
    dplyr::filter(cdf == expected_cdf) |>
    dplyr::group_by(cdf) |>
    dplyr::summarise(se = mean(se)) |>
    dplyr::transmute(se) |> as.numeric()
}
# Compute SE of logOR given number of events in each group
se_from_events <- function(
    events_treated,
    events_control,
    n_treated,
    n_control
) {
    sqrt(
        1 / events_treated +
        1 / (n_treated - events_treated) +
        1 / events_control +
        1 / (n_control - events_control)
    )
}

# Formulate empirical priors
empirical_priors <- data.frame(
    study = c("CASTLE-AF", "HF Meta-Analysis", "Non-HF Meta-Analysis"),
    mean = c(
        log(or_from_events(51, 82, 179, 184)),
        log(or_from_rr(0.52, 65 / 335)),
        log(or_from_rr(0.67, 8 / 337))
    ),
    se = c(
        se_from_events(51, 82, 179, 184),
        se_from_ci(
            log(or_from_rr(0.35, 65 / 335)),
            log(or_from_rr(0.76, 65 / 335))
        ),
        se_from_ci(
            log(or_from_rr(0.23, 8 / 337)),
            log(or_from_rr(1.99, 8 / 337))
        )
    )
)
write_xlsx(empirical_priors, "data-01_empirical-priors.xlsx")
# Formulate theoretical standard priors
standard_priors <- data.frame(
    belief = c(rep("Neutral", 3), rep("Optimistic", 3), rep("Pessimistic", 3)),
    strength = rep(c("Strong", "Moderate", "Weak"), 3),
    mean = c(
        rep(log(1), 3),
        rep(log(or_from_events(89, 101, 1108, 1096)), 3),
        rep(-log(or_from_events(89, 101, 1108, 1096)), 3)
    ),
    se = c(
        se_from_ci(log(1 / 1.5), log(1.5)),
        se_from_ci(log(1 / 2), log(2)),
        se_from_ci(log(1 / 3), log(3)),
        rep(c(
            se_from_cdf(log(or_from_events(89, 101, 1108, 1096)), 0.95, log(1)),
            se_from_cdf(log(or_from_events(89, 101, 1108, 1096)), 0.85, log(1)),
            se_from_cdf(log(or_from_events(89, 101, 1108, 1096)), 0.70, log(1))
        ), 2)
    )
)
write_xlsx(standard_priors, "data-02_standard-priors.xlsx")

# Set CABANA data
outcomes <- data.frame(
    ablation = c(
        rep(1, 1108),
        rep(0, 1096)
    ),
    primary = c(
        rep(1, 89), rep(0, 1108 - 89),
        rep(1, 101), rep(0, 1096 - 101)
    ),
    death = c(
        rep(1, 58), rep(0, 1108 - 58),
        rep(1, 67), rep(0, 1096 - 67)
    )
)
# Set priors
neutral_strong <- prior(normal(0, 0.207), class = "b")
neutral_moderate <- prior(normal(0, 0.354), class = "b")
neutral_weak <- prior(normal(0, 0.561), class = "b")
optimistic_strong <- prior(normal(-0.15, 0.091), class = "b")
optimistic_moderate <- prior(normal(-0.15, 0.145), class = "b")
optimistic_weak <- prior(normal(-0.15, 0.287), class = "b")
pessimistic_strong <- prior(normal(0.15, 0.091), class = "b")
pessimistic_moderate <- prior(normal(0.15, 0.145), class = "b")
pessimistic_weak <- prior(normal(0.15, 0.287), class = "b")
# Get effective sample sizes for theoretical priors
brm(
    formula = primary ~ ablation,
    prior = neutral_strong,
    data = outcomes,
    family = "bernoulli",
    seed = 123
) |> summary() # ESS = 3454
brm(
    formula = primary ~ ablation,
    prior = neutral_moderate,
    data = outcomes,
    family = "bernoulli",
    seed = 123
) |> summary() # ESS = 3331
brm(
    formula = primary ~ ablation,
    prior = neutral_weak,
    data = outcomes,
    family = "bernoulli",
    seed = 123
) |> summary() # ESS = 3615
brm(
    formula = primary ~ ablation,
    prior = optimistic_strong,
    data = outcomes,
    family = "bernoulli",
    seed = 123
) |> summary() # ESS = 3646
brm(
    formula = primary ~ ablation,
    prior = optimistic_moderate,
    data = outcomes,
    family = "bernoulli",
    seed = 123
) |> summary() # ESS = 3750
brm(
    formula = primary ~ ablation,
    prior = optimistic_weak,
    data = outcomes,
    family = "bernoulli",
    seed = 123
) |> summary() # ESS = 3278
brm(
    formula = primary ~ ablation,
    prior = pessimistic_strong,
    data = outcomes,
    family = "bernoulli",
    seed = 123
) |> summary() # ESS = 3232
brm(
    formula = primary ~ ablation,
    prior = pessimistic_moderate,
    data = outcomes,
    family = "bernoulli",
    seed = 123
) |> summary() # ESS = 3703
brm(
    formula = primary ~ ablation,
    prior = pessimistic_weak,
    data = outcomes,
    family = "bernoulli",
    seed = 123
) |> summary() # ESS = 3114
# Set data for effective sample sizes
data_text <- data.frame(
    belief = c("Neutral", "Optimistic", "Pessimistic"),
    x = c(-0.8, 0.2, -0.8), y = rep(3.2, 3),
    label = c(
        "ESS:\nStrong, 3454\nModerate, 3331\nWeak, 3615",
        "ESS:\nStrong, 3646\nModerate, 3750\nWeak, 3278",
        "ESS:\nStrong, 3232\nModerate, 3703\nWeak, 3114"
    )
)

# Plot empirical priors
plot_empirical <- data.frame(
    study = as.vector(apply(empirical_priors[1], 1, function(x) rep(x, 10^6))),
    x = apply(empirical_priors[2:3], 1, function(x) rnorm(10^6, x[1], x[2])) |>
    as.vector()
) |>
ggplot(aes(x, linetype = study)) +
stat_density(geom = "line", position = "identity") +
geom_vline(xintercept = 0) +
labs(x = "\nln(OR)",  y = "Density\n", linetype = "Study") +
scale_linetype_manual(
    values = c("solid", "dashed", "dotted"),
    limits = c("CASTLE-AF", "HF Meta-Analysis", "Non-HF Meta-Analysis")
) +
scale_x_continuous(expand = c(0, 0)) +
scale_y_continuous(limits = c(0, 2), breaks = c(0, 1, 2), expand = c(0, 0)) +
theme_bw() +
theme(
    legend.position = c(0.8, 0.5),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.background = element_rect(color = "black"),
    plot.margin = margin(0.3, 0.3, 0.3, 0.3, "in"),
    axis.title = element_text(face = "bold", size = 14),
    strip.text = element_text(size = 12),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10)
)
# Plot theoretical priors
plot_standard <- data.frame(
    belief = as.vector(apply(standard_priors[1], 1, function(x) rep(x, 10^6))),
    strength = apply(standard_priors[2], 1, function(x) rep(x, 10^6)) |>
    as.vector(),
    x = apply(standard_priors[3:4], 1, function(x) rnorm(10^6, x[1], x[2])) |>
    as.vector()
) |>
ggplot(aes(x, linetype = strength)) +
stat_density(geom = "line", position = "identity") +
geom_vline(xintercept = 0) +
facet_wrap(belief ~ ., dir = "v", scales = "free") +
labs(x = "\nln(OR)",  y = "Density\n", linetype = "Strength") +
scale_linetype_manual(
    values = c("solid", "dashed", "dotted"),
    limits = c("Strong", "Moderate", "Weak")
) +
geom_text(
  data = data_text,
  mapping = aes(x, y, label = label),
  inherit.aes = FALSE,
  hjust = 0,
  size = 3.5
) +
scale_x_continuous(limits = c(-1, 1), expand = c(0, 0)) +
scale_y_continuous(limits = c(0, 5), breaks = seq(0, 5), expand = c(0, 0)) +
theme_bw() +
theme(
    legend.position = "top",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.background = element_rect(color = "black"),
    plot.margin = margin(0.3, 0.3, 0.3, 0.3, "in"),
    axis.title = element_text(face = "bold", size = 14),
    strip.text = element_text(size = 12),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10),
    panel.spacing = unit(1.2, "lines")
)
# Join panels
pdf("figure-01_priors-used.pdf", width = 12, height = 12)
plot_grid(
    plot_empirical, plot_standard,
    labels = "AUTO",
    rel_widths = c(1.5, 1)
)
dev.off()