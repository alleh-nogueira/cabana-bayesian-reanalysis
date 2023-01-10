set.seed(1)
options(digits = 2)
if (!require("pacman")) install.packages("pacman")
p_load("tidyverse", "cowplot", "brms", "insight")

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

# Set prior
flat <- prior(normal(0, 100), class = "b")

# Fit Bayesian models
primary_analysis <- brm(
    formula = primary ~ ablation,
    prior = flat,
    data = outcomes,
    family = "bernoulli",
    seed = 123
)
death_analysis <- brm(
    formula = death ~ ablation,
    prior = flat,
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

# Set the threshold for clinically relevant events
threshold <- log(or_from_events(89, 101, 1108, 1096))

# Compute probabilities of effects for primary endpoint
pr_effects(primary_analysis, threshold)
# Plot posterior distribution for primary endpoint
plot_primary <-
data.frame(density(get_parameters(primary_analysis)[, 2])[c("x", "y")]) |>
mutate(
    region = case_when(
        (x < threshold) ~ "Pr(Major Benefit) = 51%",
        (threshold <= x & x < log(1)) ~ "Pr(Minor Benefit) = 33%",
        (log(1) <= x & x < -threshold) ~ "Pr(Minor Harm) = 14%",
        (-threshold <= x) ~ "Pr(Major Harm) = 2%",
        TRUE ~ "Other"
    )
) |>
ggplot(aes(x, y, fill = region)) +
geom_line() + geom_area() +
scale_fill_manual(
    values = c("#638CE3", "#779BE7", "#A480CF", "#946BC7"),
    limits = c(
        "Pr(Major Benefit) = 51%",
        "Pr(Minor Benefit) = 33%",
        "Pr(Minor Harm) = 14%",
        "Pr(Major Harm) = 2%"
    )
) +
geom_vline(xintercept = 0) +
labs(x = "\nln(OR)",  y = "Density\n", fill = "") +
scale_x_continuous(
    limits = c(-1, 0.6),
    breaks = seq(-1, 0.6, 0.2),
    expand = c(0, 0)
) +
scale_y_continuous(limits = c(0, 4), expand = c(0, 0)) +
theme_cowplot() +
theme(
    legend.position = c(0.1, 0.8),
    plot.background = element_rect(color = "black"),
    plot.margin = margin(0.3, 0.3, 0.3, 0.3, "in"),
    axis.title = element_text(face = "bold", size = 14),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10)
)

# Compute probabilities of effects for mortality prevention
pr_effects(death_analysis, threshold)
# Plot posterior distribution for primary endpoint
plot_death <-
data.frame(density(get_parameters(death_analysis)[, 2])[c("x", "y")]) |>
mutate(
    region = case_when(
        (x < threshold) ~ "Pr(Major Benefit) = 54%",
        (threshold <= x & x < log(1)) ~ "Pr(Minor Benefit) = 28%",
        (log(1) <= x & x < -threshold) ~ "Pr(Minor Harm) = 14%",
        (-threshold <= x) ~ "Pr(Major Harm) = 4%",
        TRUE ~ "Other"
    )
) |>
ggplot(aes(x, y, fill = region)) +
geom_line() + geom_area() +
scale_fill_manual(
    values = c("#638CE3", "#779BE7", "#A480CF", "#946BC7"),
    limits = c(
        "Pr(Major Benefit) = 54%",
        "Pr(Minor Benefit) = 28%",
        "Pr(Minor Harm) = 14%",
        "Pr(Major Harm) = 4%"
    )
) +
geom_vline(xintercept = 0) +
labs(x = "\nln(OR)",  y = "Density\n", fill = "") +
scale_x_continuous(
    limits = c(-1, 0.6),
    breaks = seq(-1, 0.6, 0.2),
    expand = c(0, 0)
) +
scale_y_continuous(limits = c(0, 4), expand = c(0, 0)) +
theme_cowplot() +
theme(
    legend.position = c(0.1, 0.8),
    plot.background = element_rect(color = "black"),
    plot.margin = margin(0.3, 0.3, 0.3, 0.3, "in"),
    axis.title = element_text(face = "bold", size = 14),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10)
)

# Join panels
pdf("figure-02_posteriors-non-informative.pdf", width = 10, height = 8)
plot_grid(
plot_primary, plot_death,
labels = "AUTO"
)
dev.off()