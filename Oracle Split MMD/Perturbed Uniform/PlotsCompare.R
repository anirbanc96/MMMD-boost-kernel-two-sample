################################################################################
##################### Code for generating all the plots ########################


library(tidyverse)

# read the power values
multi1 <- read_csv("MultiPower-LAP.csv")[,-c(1,2)]
multi2 <- read_csv("MultiPower-GEXP.csv")[,-c(1,2)]
multi3 <- read_csv("MultiPower-MIXED.csv")[,-c(1,2)]
single1 <- read_csv("SinglePower-LAP.csv")[,-c(1,2)]
single2 <- read_csv("SinglePower-GAUSS.csv")[,-c(1,2)]
FR <- read_csv("Power-FR.csv")[,-c(1,2)]

# read MMDAgg-Paper results

MMDAgg.results <- read_csv("CopyOfresults.csv")
MMDAgg.gaussian.power <- MMDAgg.results %>%
  filter(experiment == "1: uniform alternative") %>%
  filter(kernel_type == "gaussian") %>%
  filter(l_minus_l_plus == "(-6, -2)") %>%
  filter(function_type == "increasing") %>%
  filter(d == 1) %>%
  select(c(perturbation_or_Qi,power))

MMDAgg.laplace.power <- MMDAgg.results %>%
  filter(experiment == "1: uniform alternative") %>%
  filter(kernel_type == "laplace") %>%
  filter(l_minus_l_plus == "(-6, -2)") %>%
  filter(function_type == "increasing") %>%
  filter(d == 1) %>%
  select(c(perturbation_or_Qi,power))

# mean of power values

multi1.mean <- apply(multi1, 1, mean)
multi2.mean <- apply(multi2, 1, mean)
multi3.mean <- apply(multi3, 1, mean)
single1.mean <- apply(single1, 1, mean)
single2.mean <- apply(single2, 1, mean)
FR.mean <- apply(FR, 1, mean)

# sd of power values
#multi1.sd <- apply(multi1, 1, sd)
#multi2.sd <- apply(multi2, 1, sd)
#multi3.sd <- apply(multi3, 1, sd)
#single1.sd <- apply(single1, 1, sd)
#single2.sd <- apply(single2, 1, sd)
#FR.sd <- apply(FR, 1, sd)


# upper and lower ranges of multiple kernel power
#multi1.up <- multi1.mean + multi1.sd
#multi1.low <- multi1.mean - multi1.sd

#multi2.up <- multi2.mean + multi2.sd
#multi2.low <- multi2.mean - multi2.sd

#multi3.up <- multi3.mean + multi3.sd
#multi3.low <- multi3.mean - multi3.sd

# upper and lower ranges of single kernel power
#single1.low <- single1.mean - single1.sd
#single1.up <- single1.mean + single1.sd

#single2.low <- single2.mean - single2.sd
#single2.up <- single2.mean + single2.sd

# upper and lower ranges of graph tests
#FR.up <- FR.mean + FR.sd
#FR.low <- FR.mean - FR.sd

d <- (read_csv("MultiPower-GEXP.csv")[2])%>%pull(1)

power.tibble <- tibble(power = c(single1.mean,single2.mean, multi1.mean,
                                 multi2.mean, multi3.mean, FR.mean,
                                 MMDAgg.gaussian.power$power, MMDAgg.laplace.power$power),
                       #up = c(single1.up,single2.up, multi1.up, multi2.up,
                      #        multi3.up, FR.up),
                       #low = c(single1.low,single2.low, multi1.low, multi2.low,
                      #         multi3.low, FR.low),
                       dim = c(d,d,d,d,d,d,d,d),
                       group = c(rep("Single Laplace Kernel", length(d)),
                                 rep("Single Gaussian Kernel", length(d)),
                                 rep("Multiple Laplace Kernels", length(d)),
                                 rep("Multiple Gaussian Kernels", length(d)),
                                 rep("Multiple Mixed Kernels", length(d)),
                                 rep("Friedman Rafsky", length(d)),
                                 rep("MMDAgg Increasing Gaussian", length(d)),
                                 rep("MMDAgg Increasing Laplace", length(d))))

# plotting the power along with sd
(power.plot <- power.tibble %>%
    ggplot(aes(x = dim, y = power, group = group, col = group, fill = group)) +
    geom_point(size = 2) +
    geom_line(size = 1) +
    #geom_ribbon(aes(ymin = low, ymax = up), alpha = 0.3, linetype = 0)+
    labs(x = "Number of Perturbations", y = "Power", color = "Test", fill = "Test") +
    theme_bw()
)

# saving the plot
ggsave(plot = power.plot, filename = "MMD-Agg-Plot.pdf", device = "pdf",
       width = 5, height = 3)

# Power table
power.table <- tibble("Number of Perturbations" = d,
                      "Friedman Rafsky" = FR.mean,
                      "Single Laplace Kernel" = single1.mean,
                      "Single Gaussian Kernel" = single2.mean,
                      "Multiple Laplace Kernel" = multi1.mean,
                      "Multiple Gaussian Kernel" = multi2.mean,
                      "Multiple Mixed Kernel" = multi3.mean,
                      "MMDAgg Increasing Gaussian" = MMDAgg.gaussian.power$power,
                      "MMDAgg Increasing Laplace" = MMDAgg.laplace.power$power)
library(kableExtra)
power.table %>%
  kableExtra::kable(format = "latex",
                    row.names = NA, booktabs = T, linesep = "", escape = F,
                    digits = 2, align = rep("c",4)) %>%
    column_spec(column = 1:9, width = "0.7in") %>%
  kableExtra::save_kable("MMDAgg-Table.pdf")
