set.seed(1)
options(digits = 2)
if (!require("pacman")) install.packages("pacman")
p_load("tidyverse", "readxl", "cowplot", "brms", "insight")

# Set data
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
read_xlsx("data-02_standard-priors.xlsx")
neutral_strong <- prior(normal(0, 0.207), class = "b")
neutral_moderate <- prior(normal(0, 0.354), class = "b")
neutral_weak <- prior(normal(0, 0.561), class = "b")
optimistic_strong <- prior(normal(-0.15, 0.091), class = "b")
optimistic_moderate <- prior(normal(-0.15, 0.145), class = "b")
optimistic_weak <- prior(normal(-0.15, 0.287), class = "b")
pessimistic_strong <- prior(normal(0.15, 0.091), class = "b")
pessimistic_moderate <- prior(normal(0.15, 0.145), class = "b")
pessimistic_weak <- prior(normal(0.15, 0.287), class = "b")

# Fit Bayesian models for primary endpoint
fit_neutral_strong <- brm(
    formula = primary ~ ablation,
    prior = neutral_strong,
    data = outcomes,
    family = "bernoulli",
    seed = 123
)
fit_neutral_moderate <- brm(
    formula = primary ~ ablation,
    prior = neutral_moderate,
    data = outcomes,
    family = "bernoulli",
    seed = 123
)
fit_neutral_weak <- brm(
    formula = primary ~ ablation,
    prior = neutral_weak,
    data = outcomes,
    family = "bernoulli",
    seed = 123
)
fit_optimistic_strong <- brm(
    formula = primary ~ ablation,
    prior = optimistic_strong,
    data = outcomes,
    family = "bernoulli",
    seed = 123
)
fit_optimistic_moderate <- brm(
    formula = primary ~ ablation,
    prior = optimistic_moderate,
    data = outcomes,
    family = "bernoulli",
    seed = 123
)
fit_optimistic_weak <- brm(
    formula = primary ~ ablation,
    prior = optimistic_weak,
    data = outcomes,
    family = "bernoulli",
    seed = 123
)
fit_pessimistic_strong <- brm(
    formula = primary ~ ablation,
    prior = pessimistic_strong,
    data = outcomes,
    family = "bernoulli",
    seed = 123
)
fit_pessimistic_moderate <- brm(
    formula = primary ~ ablation,
    prior = pessimistic_moderate,
    data = outcomes,
    family = "bernoulli",
    seed = 123
)
fit_pessimistic_weak <- brm(
    formula = primary ~ ablation,
    prior = pessimistic_weak,
    data = outcomes,
    family = "bernoulli",
    seed = 123
)

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
# Compute cdf (ie, Pr(logOR < threshold)) from Bayesian model
pr_benefit <- function(model, threshold = log(1)) {
    round(abs(pnorm(
        threshold,
        summary(model)$fixed[2, 1],
        summary(model)$fixed[2, 2]
    )) * 100, 0)
}
# Compute probabilities of benefit and harm
pr_effects <- function(model, threshold) {c(
    # Compute probabilities of benefit
    any_benefit = pr_benefit(model),
    major_benefit = pr_benefit(model, threshold),
    minor_benefit = pr_benefit(model) - pr_benefit(model, threshold),
    # Compute probabilities of harm
    any_harm = 100 - pr_benefit(model),
    major_harm = 100 - pr_benefit(model, -threshold),
    minor_harm = pr_benefit(model, -threshold) - pr_benefit(model)
)}
# Compute RR from OR
rr_from_or <- function(odds_ratio, untreated_risk) {round(
    odds_ratio / ((1 - untreated_risk) + untreated_risk * odds_ratio), 2
)}

# Set the threshold for clinically relevant events
threshold <- log(or_from_events(89, 101, 1108, 1096))
# Get data to plot standard priors
standard_priors <- read_xlsx("data-02_standard-priors.xlsx")
data_standard <- data.frame(
    belief = as.vector(apply(standard_priors[1], 1, function(x) rep(x, 10^6))),
    strength = apply(standard_priors[2], 1, function(x) rep(x, 10^6)) |>
    as.vector(),
    x = apply(standard_priors[3:4], 1, function(x) rnorm(10^6, x[1], x[2])) |>
    as.vector()
)
# Get data for posteriors
data_posteriors <- rbind(
data.frame(density(get_parameters(fit_neutral_strong)[, 2])[c("x", "y")]),
data.frame(density(get_parameters(fit_neutral_moderate)[, 2])[c("x", "y")]),
data.frame(density(get_parameters(fit_neutral_weak)[, 2])[c("x", "y")]),
data.frame(density(get_parameters(fit_optimistic_strong)[, 2])[c("x", "y")]),
data.frame(density(get_parameters(fit_optimistic_moderate)[, 2])[c("x", "y")]),
data.frame(density(get_parameters(fit_optimistic_weak)[, 2])[c("x", "y")]),
data.frame(density(get_parameters(fit_pessimistic_strong)[, 2])[c("x", "y")]),
data.frame(density(get_parameters(fit_pessimistic_moderate)[, 2])[c("x", "y")]),
data.frame(density(get_parameters(fit_pessimistic_weak)[, 2])[c("x", "y")])
) |>
mutate(region = case_when(
        (x < threshold) ~ "Major Benefit",
        (threshold <= x & x < log(1)) ~ "Minor Benefit",
        (log(1) <= x & x < -threshold) ~ "Minor Harm",
        (-threshold <= x) ~ "Major Harm",
        TRUE ~ "Other"
))
data_posteriors$belief <- c(
    rep("Neutral", nrow(data_posteriors) / 3),
    rep("Optimistic", nrow(data_posteriors) / 3),
    rep("Pessimistic", nrow(data_posteriors) / 3)
)
data_posteriors$strength <- c(rep(c(
    rep("Strong", nrow(data_posteriors) / 9),
    rep("Moderate", nrow(data_posteriors) / 9),
    rep("Weak", nrow(data_posteriors) / 9))
, 3))
# Get # Get data for annotations
summary(fit_neutral_strong)
pr_effects(fit_neutral_strong, threshold)
summary(fit_neutral_moderate)
pr_effects(fit_neutral_moderate, threshold)
summary(fit_neutral_weak)
pr_effects(fit_neutral_weak, threshold)
summary(fit_optimistic_strong)
pr_effects(fit_optimistic_strong, threshold)
summary(fit_optimistic_moderate)
pr_effects(fit_optimistic_moderate, threshold)
summary(fit_optimistic_weak)
pr_effects(fit_optimistic_weak, threshold)
summary(fit_pessimistic_strong)
pr_effects(fit_pessimistic_strong, threshold)
summary(fit_pessimistic_moderate)
pr_effects(fit_pessimistic_moderate, threshold)
summary(fit_pessimistic_weak)
pr_effects(fit_pessimistic_weak, threshold)
data_text <- data.frame(
    belief = c(rep("Neutral", 3), rep("Optimistic", 3), rep("Pessimistic", 3)),
    strength = rep(c("Strong", "Moderate", "Weak"), 3),
    x = rep(-1.8, 9), y = rep(3.2, 9),
    label = c(
        paste0( # Neutral strong
            "Median RR = ",
            rr_from_or(exp(-0.10), 101 / 1096),
            "\n95% CrI: (", rr_from_or(exp(-0.34), 101 / 1096),
            "; ", rr_from_or(exp(0.15), 101 / 1096), ")\n\n",
            "Pr(Major Benefit) = 33%\n", "Pr(Minor Benefit) = 45%\n",
            "Pr(Minor Harm) = 20%\n", "Pr(Major Harm) = 2%"
        ),
        paste0( # Neutral moderate
            "Median RR = ",
            rr_from_or(exp(-0.13), 101 / 1096),
            "\n95% CrI: (", rr_from_or(exp(-0.39), 101 / 1096),
            "; ", rr_from_or(exp(0.14), 101 / 1096), ")\n\n",
            "Pr(Major Benefit) = 43%\n", "Pr(Minor Benefit) = 39%\n",
            "Pr(Minor Harm) = 16%\n", "Pr(Major Harm) = 2%"
        ),
        paste0( # Neutral weak
            "Median RR = ",
            rr_from_or(exp(-0.14), 101 / 1096),
            "\n95% CrI: (", rr_from_or(exp(-0.43), 101 / 1096),
            "; ", rr_from_or(exp(0.15), 101 / 1096), ")\n\n",
            "Pr(Major Benefit) = 47%\n", "Pr(Minor Benefit) = 36%\n",
            "Pr(Minor Harm) = 15%\n", "Pr(Major Harm) = 2%"
        ),
        paste0( # Optimistic strong
            "Median RR = ",
            rr_from_or(exp(-0.15), 101 / 1096),
            "\n95% CrI: (", rr_from_or(exp(-0.30), 101 / 1096),
            "; ", rr_from_or(exp(-0.01), 101 / 1096), ")\n\n",
            "Pr(Major Benefit) = 50%\n", "Pr(Minor Benefit) = 47%\n",
            "Pr(Minor Harm) = 3%\n", "Pr(Major Harm) = 0%"
        ),
        paste0( # Optimistic moderate
            "Median RR = ",
            rr_from_or(exp(-0.15), 101 / 1096),
            "\n95% CrI: (", rr_from_or(exp(-0.36), 101 / 1096),
            "; ", rr_from_or(exp(0.05), 101 / 1096), ")\n\n",
            "Pr(Major Benefit) = 51%\n", "Pr(Minor Benefit) = 42%\n",
            "Pr(Minor Harm) = 7%\n", "Pr(Major Harm) = 0%"
        ),
        paste0( # Optimistic weak
            "Median RR = ",
            rr_from_or(exp(-0.15), 101 / 1096),
            "\n95% CrI: (", rr_from_or(exp(-0.42), 101 / 1096),
            "; ", rr_from_or(exp(0.12), 101 / 1096), ")\n\n",
            "Pr(Major Benefit) = 49%\n", "Pr(Minor Benefit) = 37%\n",
            "Pr(Minor Harm) = 12%\n", "Pr(Major Harm) = 2%"
        ),
        paste0( # Pessimistic strong
            "Median RR = ",
            rr_from_or(exp(0.07), 101 / 1096),
            "\n95% CrI: (", rr_from_or(exp(-0.08), 101 / 1096),
            "; ", rr_from_or(exp(0.23), 101 / 1096), ")\n\n",
            "Pr(Major Benefit) = 0%\n", "Pr(Minor Benefit) = 19%\n",
            "Pr(Minor Harm) = 65%\n", "Pr(Major Harm) = 16%"
        ),
        paste0( # Pessimistic moderate
            "Median RR = ",
            rr_from_or(exp(0), 101 / 1096),
            "\n95% CrI: (", rr_from_or(exp(-0.19), 101 / 1096),
            "; ", rr_from_or(exp(0.21), 101 / 1096), ")\n\n",
            "Pr(Major Benefit) = 6%\n", "Pr(Minor Benefit) = 42%\n",
            "Pr(Minor Harm) = 44%\n", "Pr(Major Harm) = 8%"
        ),
        paste0( # Pessimistic weak
            "Median RR = ",
            rr_from_or(exp(-0.08), 101 / 1096),
            "\n95% CrI: (", rr_from_or(exp(-0.36), 101 / 1096),
            "; ", rr_from_or(exp(0.19), 101 / 1096), ")\n\n",
            "Pr(Major Benefit) = 32%\n", "Pr(Minor Benefit) = 40%\n",
            "Pr(Minor Harm) = 23%\n", "Pr(Major Harm) = 5%"
        )
    )
)
# Plot posterior distributions for primary endpoint
pdf("figure-04_posteriors-standard-primary.pdf", width = 12, height = 12)
data_posteriors |>
ggplot(aes(x, y, fill = region)) +
geom_line(linetype = "solid") + geom_area() + geom_vline(xintercept = 0) +
facet_grid(
    factor(belief, levels = c("Optimistic", "Neutral", "Pessimistic")) ~
    factor(strength, levels = c("Strong", "Moderate", "Weak")),
    scales = "free", switch = "y",
) +
scale_fill_manual(
    values = c("#638CE3", "#779BE7", "#A480CF", "#946BC7"),
    limits = c("Major Benefit", "Minor Benefit", "Minor Harm", "Major Harm")
) +
stat_density(
    data = data_standard,
    mapping = aes(x),
    inherit.aes = FALSE,
    geom = "line",
    position = "identity",
    linetype = "dashed"
) +
labs(x = "\nln(OR)",  y = "Density\n", fill = "") +
geom_text(
  data = data_text,
  mapping = aes(x, y, label = label),
  inherit.aes = FALSE,
  hjust = 0
) +
scale_x_continuous(limits = c(-2, 1), breaks = seq(-2, 1), expand = c(0, 0)) +
scale_y_continuous(limits = c(0, 5), expand = c(0, 0), position = "right") +
theme_bw() +
theme(
    legend.position = "top",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.background = element_rect(color = "black"),
    plot.margin = margin(0.3, 0.3, 0.3, 0.3, "in"),
    axis.title = element_text(face = "bold", size = 14),
    strip.text = element_text(size = 12),
    legend.text = element_text(size = 12),
    panel.spacing = unit(1.2, "lines")
)
dev.off()