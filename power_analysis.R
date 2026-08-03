library(pwr)
library(dplyr)


data <- read.csv("data.csv", sep = ";")

# group size
n_PD <- 69
n_HC <- 118

# power analysis for each odorant
data <- data %>%
  rowwise() %>%
  mutate(
    PD_rate = PD_corr / n_PD,
    HC_rate = HC_corr / n_HC,
    effect_h = ES.h(p1 = PD_rate, p2 = HC_rate),
    power = pwr.2p2n.test(h = abs(effect_h), n1 = n_PD, n2 = n_HC,
                          sig.level = 0.05)$power
  ) %>%
  ungroup()

# output
print(data %>% select(Odor, PD_corr, HC_corr, PD_rate, HC_rate, effect_h, power))

write.csv(data, "power_analysis_all_odors.csv", row.names = FALSE)