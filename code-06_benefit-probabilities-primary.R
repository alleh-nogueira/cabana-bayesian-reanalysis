set.seed(1)
options(digits = 2)
if (!require("pacman")) install.packages("pacman")
p_load("tidyverse", "readxl", "writexl", "brms", "insight")

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

# Set flat prior
flat <- prior(normal(0, 100), class = "b")
# Set evidence-based priors
castle <- prior(normal(-0.702, 0.222), class = "b")
hf <- prior(normal(-0.763, 0.221), class = "b")
nonhf <- prior(normal(-0.408, 0.561), class = "b")
# Set theoretical priors
neutral_strong <- prior(normal(0, 0.207), class = "b")
neutral_moderate <- prior(normal(0, 0.354), class = "b")
neutral_weak <- prior(normal(0, 0.561), class = "b")
optimistic_strong <- prior(normal(-0.15, 0.091), class = "b")
optimistic_moderate <- prior(normal(-0.15, 0.145), class = "b")
optimistic_weak <- prior(normal(-0.15, 0.287), class = "b")
pessimistic_strong <- prior(normal(0.15, 0.091), class = "b")
pessimistic_moderate <- prior(normal(0.15, 0.145), class = "b")
pessimistic_weak <- prior(normal(0.15, 0.287), class = "b")

# Fit Bayesian regression of non-informative prior
fit_flat <- brm(
    formula = primary ~ ablation,
    prior = flat,
    data = outcomes,
    family = "bernoulli",
    seed = 123
)
# Fit Bayesian regression of evidence-based priors
fit_castle <- brm(
    formula = primary ~ ablation,
    prior = castle,
    data = outcomes,
    family = "bernoulli",
    seed = 123
)
fit_hf <- brm(
    formula = primary ~ ablation,
    prior = hf,
    data = outcomes,
    family = "bernoulli",
    seed = 123
)
fit_nonhf <- brm(
    formula = primary ~ ablation,
    prior = nonhf,
    data = outcomes,
    family = "bernoulli",
    seed = 123
)
# Fit Bayesian regression of theoretical priors
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

# Compute OR given RR and risk in control group
or_from_rr <- function(risk_ratio, untreated_risk = 101 / 1096) {
    risk_ratio * (1 - untreated_risk) / (1 - risk_ratio * untreated_risk)
}
# Compute cdf (ie, Pr(logOR < threshold)) from Bayesian model
pr_benefit <- function(model, threshold = log(1)) {
    round(abs(pnorm(
        threshold,
        summary(model)$fixed[2, 1],
        summary(model)$fixed[2, 2]
    )) * 100, 0)
}
# Compute RR from OR
rr_from_or <- function(odds_ratio, untreated_risk = 101 / 1096) {round(
    odds_ratio / ((1 - untreated_risk) + untreated_risk * odds_ratio), 2
)}
# Compute median RR (95% CrI) from Bayesian regression model
rr_from_brm <- function(model, untreated_risk = 101 / 1096) {
    paste0(
            rr_from_or(exp(summary(model)$fixed[2, 1]), untreated_risk), " (",
            rr_from_or(exp(summary(model)$fixed[2, 3]), untreated_risk), "; ",
            rr_from_or(exp(summary(model)$fixed[2, 4]), untreated_risk), ")"
        )
}

# Compute Pr(Benefit) of each model for various thresholds
data.frame(
    "Prior Belief and Strength" = c(
        "Noninformative",
        "CASTLE-AF",
        "HF Meta-Analysis",
        "Non-HF Meta-Analysis",
        "Optimistic Strong",
        "Optimistic Moderate",
        "Optimistic Weak",
        "Neutral Strong",
        "Neutral Moderate",
        "Neutral Weak",
        "Pessimistic Strong",
        "Pessimistic Moderate",
        "Pessimistic Weak"
    ),
    "Posterior Median RR and CrI" = c(
        rr_from_brm(fit_flat),
        rr_from_brm(fit_castle),
        rr_from_brm(fit_hf),
        rr_from_brm(fit_nonhf),
        rr_from_brm(fit_optimistic_strong),
        rr_from_brm(fit_optimistic_moderate),
        rr_from_brm(fit_optimistic_weak),
        rr_from_brm(fit_neutral_strong),
        rr_from_brm(fit_neutral_moderate),
        rr_from_brm(fit_neutral_weak),
        rr_from_brm(fit_pessimistic_strong),
        rr_from_brm(fit_pessimistic_moderate),
        rr_from_brm(fit_pessimistic_weak)
    ),
    "Pr RRR above 0" = c(
        pr_benefit(fit_flat),
        pr_benefit(fit_castle),
        pr_benefit(fit_hf),
        pr_benefit(fit_nonhf),
        pr_benefit(fit_optimistic_strong),
        pr_benefit(fit_optimistic_moderate),
        pr_benefit(fit_optimistic_weak),
        pr_benefit(fit_neutral_strong),
        pr_benefit(fit_neutral_moderate),
        pr_benefit(fit_neutral_weak),
        pr_benefit(fit_pessimistic_strong),
        pr_benefit(fit_pessimistic_moderate),
        pr_benefit(fit_pessimistic_weak)
    ),
    "Pr RRR above 0.05" = c(
        pr_benefit(fit_flat, log(or_from_rr(0.95))),
        pr_benefit(fit_castle, log(or_from_rr(0.95))),
        pr_benefit(fit_hf, log(or_from_rr(0.95))),
        pr_benefit(fit_nonhf, log(or_from_rr(0.95))),
        pr_benefit(fit_optimistic_strong, log(or_from_rr(0.95))),
        pr_benefit(fit_optimistic_moderate, log(or_from_rr(0.95))),
        pr_benefit(fit_optimistic_weak, log(or_from_rr(0.95))),
        pr_benefit(fit_neutral_strong, log(or_from_rr(0.95))),
        pr_benefit(fit_neutral_moderate, log(or_from_rr(0.95))),
        pr_benefit(fit_neutral_weak, log(or_from_rr(0.95))),
        pr_benefit(fit_pessimistic_strong, log(or_from_rr(0.95))),
        pr_benefit(fit_pessimistic_moderate, log(or_from_rr(0.95))),
        pr_benefit(fit_pessimistic_weak, log(or_from_rr(0.95)))
    ),
    "Pr RRR above 0.10" = c(
        pr_benefit(fit_flat, log(or_from_rr(0.9))),
        pr_benefit(fit_castle, log(or_from_rr(0.9))),
        pr_benefit(fit_hf, log(or_from_rr(0.9))),
        pr_benefit(fit_nonhf, log(or_from_rr(0.9))),
        pr_benefit(fit_optimistic_strong, log(or_from_rr(0.9))),
        pr_benefit(fit_optimistic_moderate, log(or_from_rr(0.9))),
        pr_benefit(fit_optimistic_weak, log(or_from_rr(0.9))),
        pr_benefit(fit_neutral_strong, log(or_from_rr(0.9))),
        pr_benefit(fit_neutral_moderate, log(or_from_rr(0.9))),
        pr_benefit(fit_neutral_weak, log(or_from_rr(0.9))),
        pr_benefit(fit_pessimistic_strong, log(or_from_rr(0.9))),
        pr_benefit(fit_pessimistic_moderate, log(or_from_rr(0.9))),
        pr_benefit(fit_pessimistic_weak, log(or_from_rr(0.9)))
    ),
    "Pr RRR above 0.15" = c(
        pr_benefit(fit_flat, log(or_from_rr(0.87))),
        pr_benefit(fit_castle, log(or_from_rr(0.85))),
        pr_benefit(fit_hf, log(or_from_rr(0.85))),
        pr_benefit(fit_nonhf, log(or_from_rr(0.85))),
        pr_benefit(fit_optimistic_strong, log(or_from_rr(0.85))),
        pr_benefit(fit_optimistic_moderate, log(or_from_rr(0.85))),
        pr_benefit(fit_optimistic_weak, log(or_from_rr(0.85))),
        pr_benefit(fit_neutral_strong, log(or_from_rr(0.85))),
        pr_benefit(fit_neutral_moderate, log(or_from_rr(0.85))),
        pr_benefit(fit_neutral_weak, log(or_from_rr(0.85))),
        pr_benefit(fit_pessimistic_strong, log(or_from_rr(0.85))),
        pr_benefit(fit_pessimistic_moderate, log(or_from_rr(0.85))),
        pr_benefit(fit_pessimistic_weak, log(or_from_rr(0.85)))
    ),
    "Pr RRR above 0.20" = c(
        pr_benefit(fit_flat, log(or_from_rr(0.8))),
        pr_benefit(fit_castle, log(or_from_rr(0.8))),
        pr_benefit(fit_hf, log(or_from_rr(0.8))),
        pr_benefit(fit_nonhf, log(or_from_rr(0.8))),
        pr_benefit(fit_optimistic_strong, log(or_from_rr(0.8))),
        pr_benefit(fit_optimistic_moderate, log(or_from_rr(0.8))),
        pr_benefit(fit_optimistic_weak, log(or_from_rr(0.8))),
        pr_benefit(fit_neutral_strong, log(or_from_rr(0.8))),
        pr_benefit(fit_neutral_moderate, log(or_from_rr(0.8))),
        pr_benefit(fit_neutral_weak, log(or_from_rr(0.8))),
        pr_benefit(fit_pessimistic_strong, log(or_from_rr(0.8))),
        pr_benefit(fit_pessimistic_moderate, log(or_from_rr(0.8))),
        pr_benefit(fit_pessimistic_weak, log(or_from_rr(0.8)))
    )
) |> write_xlsx("table-01_benefit-probabilities-primary.xlsx")