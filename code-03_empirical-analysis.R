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
read_xlsx("data-01_empirical-priors.xlsx")
castle <- prior(normal(-0.702, 0.314), class = "b")
hf <- prior(normal(-0.763, 0.312), class = "b")
nonhf <- prior(normal(-0.408, 0.561), class = "b")

# Fit Bayesian models for primary endpoint
primary_castle <- brm(
    formula = primary ~ ablation,
    prior = castle,
    data = outcomes,
    family = "bernoulli",
    seed = 123
)
primary_hf <- brm(
    formula = primary ~ ablation,
    prior = hf,
    data = outcomes,
    family = "bernoulli",
    seed = 123
)
primary_nonhf <- brm(
    formula = primary ~ ablation,
    prior = nonhf,
    data = outcomes,
    family = "bernoulli",
    seed = 123
)
# Fit Bayesian models for mortality
death_castle <- brm(
    formula = death ~ ablation,
    prior = castle,
    data = outcomes,
    family = "bernoulli",
    seed = 123
)
death_hf <- brm(
    formula = death ~ ablation,
    prior = hf,
    data = outcomes,
    family = "bernoulli",
    seed = 123
)
death_nonhf <- brm(
    formula = death ~ ablation,
    prior = nonhf,
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
# Get data for plotting priors
empirical_priors <- read_xlsx("data-01_empirical-priors.xlsx")
data_empirical <- data.frame(
    study = as.vector(apply(empirical_priors[1], 1, function(x) rep(x, 10^6))),
    x = apply(empirical_priors[2:3], 1, function(x) rnorm(10^6, x[1], x[2])) |>
    as.vector()
)

# Get data for primary endpoint
data_primary <- rbind(
data.frame(density(get_parameters(primary_castle)[, 2])[c("x", "y")]),
data.frame(density(get_parameters(primary_hf)[, 2])[c("x", "y")]),
data.frame(density(get_parameters(primary_nonhf)[, 2])[c("x", "y")])
) |>
mutate(region = case_when(
        (x < threshold) ~ "Major Benefit",
        (threshold <= x & x < log(1)) ~ "Minor Benefit",
        (log(1) <= x & x < -threshold) ~ "Minor Harm",
        (-threshold <= x) ~ "Major Harm",
        TRUE ~ "Other"
))
data_primary$study <- c(
    rep("CASTLE-AF", nrow(data_primary) / 3),
    rep("HF Meta-Analysis", nrow(data_primary) / 3),
    rep("Non-HF Meta-Analysis", nrow(data_primary) / 3)
)
# Get data for annotations
summary(primary_castle)
summary(primary_hf)
summary(primary_nonhf)
pr_effects(primary_castle, threshold)
pr_effects(primary_hf, threshold)
pr_effects(primary_nonhf, threshold)
data_text <- data.frame(
    study = c("CASTLE-AF", "HF Meta-Analysis", "Non-HF Meta-Analysis"),
    x = rep(-1.7, 3), y = rep(3.2, 3),
    label = c(
        paste0(
            "Median RR = ",
            rr_from_or(exp(-0.25), 101 / 1096),
            "\n95% CrI: (", rr_from_or(exp(-0.52), 101 / 1096),
            "; ", rr_from_or(exp(0.02), 101 / 1096), ")\n\n",
            "Pr(Major Benefit) = 77%\n", "Pr(Minor Benefit) = 20%\n",
            "Pr(Minor Harm) = 3%\n", "Pr(Major Harm) = 0%"
        ),
        paste0(
            "Median RR = ",
            rr_from_or(exp(-0.27), 101 / 1096),
            "\n95% CrI: (", rr_from_or(exp(-0.54), 101 / 1096),
            "; ", rr_from_or(exp(0), 101 / 1096), ")\n\n",
            "Pr(Major Benefit) = 80%\n", "Pr(Minor Benefit) = 17%\n",
            "Pr(Minor Harm) = 3%\n", "Pr(Major Harm) = 0%"
        ),
        paste0(
            "Median RR = ", 
            rr_from_or(exp(-0.18), 101 / 1096),
            "\n95% CrI: (", rr_from_or(exp(-0.46), 101 / 1096),
            "; ", rr_from_or(exp(0.11), 101 / 1096), ")\n\n",
            "Pr(Major Benefit) = 57%\n", "Pr(Minor Benefit) = 31%\n",
            "Pr(Minor Harm) = 12%\n", "Pr(Major Harm) = 1%"
        )
    )
)
# Plot posterior distributions for primary endpoint
posteriors_primary <- data_primary |>
ggplot(aes(x, y, fill = region)) +
geom_line(linetype = "solid") + geom_area() + geom_vline(xintercept = 0) +
facet_wrap(. ~ study, dir = "h", scales = "free") +
scale_fill_manual(
    values = c("#638CE3", "#779BE7", "#A480CF", "#946BC7"),
    limits = c("Major Benefit", "Minor Benefit", "Minor Harm", "Major Harm")
) +
stat_density(
    data = data_empirical,
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
scale_x_continuous(limits = c(-2, 1), expand = c(0, 0)) +
scale_y_continuous(limits = c(0, 4), expand = c(0, 0)) +
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

# Get data for mortality
data_death <- rbind(
data.frame(density(get_parameters(death_castle)[, 2])[c("x", "y")]),
data.frame(density(get_parameters(death_hf)[, 2])[c("x", "y")]),
data.frame(density(get_parameters(death_nonhf)[, 2])[c("x", "y")])
) |>
mutate(region = case_when(
        (x < threshold) ~ "Major Benefit",
        (threshold <= x & x < log(1)) ~ "Minor Benefit",
        (log(1) <= x & x < -threshold) ~ "Minor Harm",
        (-threshold <= x) ~ "Major Harm",
        TRUE ~ "Other"
))
data_death$study <- c(
    rep("CASTLE-AF", nrow(data_death) / 3),
    rep("HF Meta-Analysis", nrow(data_death) / 3),
    rep("Non-HF Meta-Analysis", nrow(data_death) / 3)
)
# Get data for annotations
summary(death_castle)
summary(death_hf)
summary(death_nonhf)
pr_effects(death_castle, threshold)
pr_effects(death_hf, threshold)
pr_effects(death_nonhf, threshold)
data_text <- data.frame(
    study = c("CASTLE-AF", "HF Meta-Analysis", "Non-HF Meta-Analysis"),
    x = rep(-1.7, 3), y = rep(3.2, 3),
    label = c(
        paste0(
            "Median RR = ",
            rr_from_or(exp(-0.30), 101 / 1096),
            "\n95% CrI: (", rr_from_or(exp(-0.63), 101 / 1096),
            "; ", rr_from_or(exp(0.02), 101 / 1096), ")\n\n",
            "Pr(Major Benefit) = 83%\n", "Pr(Minor Benefit) = 14%\n",
            "Pr(Minor Harm) = 3%\n", "Pr(Major Harm) = 0%"
        ),
        paste0(
            "Median RR = ",
            rr_from_or(exp(-0.32), 101 / 1096),
            "\n95% CrI: (", rr_from_or(exp(-0.65), 101 / 1096),
            "; ", rr_from_or(exp(-0.01), 101 / 1096), ")\n\n",
            "Pr(Major Benefit) = 85%\n", "Pr(Minor Benefit) = 12%\n",
            "Pr(Minor Harm) = 3%\n", "Pr(Major Harm) = 0%"
        ),
        paste0(
            "Median RR = ",
            rr_from_or(exp(-0.19), 101 / 1096),
            "\n95% CrI: (", rr_from_or(exp(-0.54), 101 / 1096),
            "; ", rr_from_or(exp(0.15), 101 / 1096), ")\n\n",
            "Pr(Major Benefit) = 59%\n", "Pr(Minor Benefit) = 27%\n",
            "Pr(Minor Harm) = 11%\n", "Pr(Major Harm) = 3%"
        )
    )
)
# Plot posterior distributions for primary endpoint
posteriors_death <- data_primary |>
ggplot(aes(x, y, fill = region)) +
geom_line(linetype = "solid") + geom_area() + geom_vline(xintercept = 0) +
facet_wrap(. ~ study, dir = "h", scales = "free") +
scale_fill_manual(
    values = c("#638CE3", "#779BE7", "#A480CF", "#946BC7"),
    limits = c("Major Benefit", "Minor Benefit", "Minor Harm", "Major Harm")
) +
stat_density(
    data = data_empirical,
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
scale_x_continuous(limits = c(-2, 1), expand = c(0, 0)) +
scale_y_continuous(limits = c(0, 4), expand = c(0, 0)) +
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

pdf("figure-03_posteriors-empirical.pdf", width = 14, height = 14)
plot_grid(
    posteriors_primary, posteriors_death,
    labels = "AUTO", ncol = 1
)
dev.off()