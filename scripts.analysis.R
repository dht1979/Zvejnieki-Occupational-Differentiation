## ============================================================
## Femora Pendant vs None × Sex (Female vs Male)
## Bayesian 95% credible intervals (mean dot + CI bars)
## Layout: Ix | Iy (row1), J | (IxIy + legend) (row2)
## Vertical PNG @ 600 dpi
## ============================================================

library(dplyr)
library(tidyr)
library(ggplot2)
library(readr)

if (!requireNamespace("patchwork", quietly = TRUE)) install.packages("patchwork")
if (!requireNamespace("forcats", quietly = TRUE)) install.packages("forcats")
library(patchwork)
library(forcats)

## ---- BASE DIRECTORY ----
base_dir <- "C:whatever your path is"

## File in base_dir (preferred)
in_path <- file.path(base_dir, "Femora Pendant All.csv")

## If you are running from the uploaded sandbox copy instead, use this:
# in_path <- "/mnt/data/Femora Pendant All.csv"

out_png <- file.path(base_dir, "Femora_PendantVsNone_BySex_CrI_Vertical.png")

## ============================================================
## Read + standardize
## Expected columns: Group, Sex, plus numeric vars (Ix, Iy, J, IxIy)
## ============================================================

dat_raw <- read_csv(in_path, show_col_types = FALSE)

names(dat_raw) <- trimws(names(dat_raw))
names(dat_raw) <- gsub("^Ix\\s*/\\s*Iy$", "IxIy", names(dat_raw))

stopifnot(all(c("Group", "Sex") %in% names(dat_raw)))

dat <- dat_raw %>%
  mutate(
    Group = as.factor(Group),
    Sex   = as.factor(Sex)
  )

## Enforce Sex order Female then Male (handles common variants)
dat$Sex <- fct_recode(dat$Sex,
                      Female = "F", Female = "f", Female = "female", Female = "FEMALE",
                      Male   = "M", Male   = "m", Male   = "male",   Male   = "MALE")
dat$Sex <- factor(dat$Sex, levels = c("Female", "Male"))

## Enforce Group order None then Pendant if present
dat$Group <- fct_recode(dat$Group,
                        None    = "none", None = "NONE",
                        Pendant = "pendant", Pendant = "PENDANT")
if (all(c("None","Pendant") %in% levels(dat$Group))) {
  dat$Group <- fct_relevel(dat$Group, "None", "Pendant")
}

## Numeric variables
vars <- setdiff(names(dat), c("Group", "Sex"))
dat[vars] <- lapply(dat[vars], function(x) suppressWarnings(as.numeric(x)))

## ============================================================
## Bayesian posterior mean sampler (Normal + Inv-Gamma; weak priors)
## ============================================================

rinvgamma <- function(n, shape, rate) 1 / rgamma(n, shape = shape, rate = rate)

posterior_mu_draws <- function(y, n_draws = 4000) {
  y <- y[is.finite(y)]
  n <- length(y)
  if (n < 2) return(rep(NA_real_, n_draws))

  ybar <- mean(y)
  s2   <- var(y)

  kappa0 <- 0.001
  alpha0 <- 0.001
  beta0  <- max(s2, 1e-8) * 0.001

  kappa_n <- kappa0 + n
  mu_n    <- ybar
  alpha_n <- alpha0 + n/2
  beta_n  <- beta0 + 0.5 * (n - 1) * s2

  sig2 <- rinvgamma(n_draws, alpha_n, beta_n)
  rnorm(n_draws, mean = mu_n, sd = sqrt(sig2 / kappa_n))
}

set.seed(1)
N_DRAWS <- 4000

## ============================================================
## Long + summarize
## ============================================================

long <- dat %>%
  pivot_longer(cols = all_of(vars),
               names_to = "Variable",
               values_to = "Value") %>%
  filter(is.finite(Value))

summ <- long %>%
  group_by(Sex, Group, Variable) %>%
  summarize(
    mu_draw = list(posterior_mu_draws(Value, N_DRAWS)),
    .groups = "drop"
  ) %>%
  mutate(
    mean  = sapply(mu_draw, mean, na.rm = TRUE),
    lower = sapply(mu_draw, quantile, probs = 0.025, na.rm = TRUE),
    upper = sapply(mu_draw, quantile, probs = 0.975, na.rm = TRUE)
  ) %>%
  select(Sex, Group, Variable, mean, lower, upper)

## Keep variable order Ix, Iy, J, IxIy if present
requested <- c("Ix","Iy","J","IxIy")
present <- unique(as.character(summ$Variable))
var_order <- requested[requested %in% present]
summ$Variable <- factor(as.character(summ$Variable), levels = var_order)

## ============================================================
## Panel function; I made the legends in Photoshop, too annoying to do here
## ============================================================

make_panel <- function(v) {
  df <- summ %>% filter(Variable == v)
  dodge <- position_dodge(width = 0.30)

  ggplot(df, aes(x = Sex, y = mean)) +
    geom_errorbar(
      aes(ymin = lower, ymax = upper, color = Group),
      position = dodge, width = 0.10, linewidth = 0.9
    ) +
    geom_point(
      aes(group = Group),
      position = dodge, size = 2.5, color = "black"
    ) +
    labs(x = NULL, y = v, color = "Group") +
    theme_classic(base_size = 11) +
    theme(
      axis.line  = element_line(linewidth = 0.6),
      axis.ticks = element_line(linewidth = 0.6),
      axis.text  = element_text(size = 10),
      axis.title = element_text(size = 11, face = "bold"),
      plot.margin = margin(6, 10, 6, 6),
      legend.position = "none"
    )
}

p_Ix   <- make_panel("Ix")
p_Iy   <- make_panel("Iy")
p_J    <- make_panel("J")
p_IxIy <- make_panel("IxIy")

## ============================================================
## Legend panel (blank space right of IxIy for the legend)
## ============================================================

legend_plot <- ggplot(
  summ %>% filter(Variable == "Ix"),
  aes(x = Sex, y = mean, color = Group)
) +
  geom_point(size = 2.6) +
  labs(color = "Mortuary group") +
  theme_void(base_size = 12) +
  theme(
    legend.position = "left",
    legend.title = element_text(face = "bold"),
    legend.text = element_text(size = 11),
    legend.key = element_blank()
  )

legend_g <- patchwork::wrap_elements(full = patchwork::get_legend(legend_plot))

legend_panel <- ggplot() + theme_void()
legend_panel <- legend_panel +
  inset_element(
    legend_g,
    left = 0.15, right = 0.95,
    bottom = 0.30, top = 0.70
  )

## ============================================================
## Layout: Ix|Iy; J | (IxIy + legend)
## ============================================================

row1 <- p_Ix | p_Iy
row2 <- p_J  | (p_IxIy | legend_panel) + plot_layout(widths = c(2, 1))

fig <- row1 / row2

## ============================================================
## Save PNG (vertical)
## ============================================================

ggsave(out_png, fig, width = 8.5, height = 11, units = "in", dpi = 600)

print(fig)
message("Saved: ", out_png)

