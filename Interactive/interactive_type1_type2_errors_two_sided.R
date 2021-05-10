##########################################################
# TWO-SIDED TEST - INTERACTIVE TYPE I AND II ERROR PLOT
##########################################################


library(manipulate)
library(ggplot2)


# Dynamic ggplot
two_sided <- function(sigma, mean_H1, n, alpha) {
  
  mean_H0 <- 0
  std_err <- sigma / sqrt(n)
  
  # Create the ggplot
  xmin <- min(mean_H0, mean_H1) - 4 * std_err
  xmax <- max(mean_H0, mean_H1) + 4 * std_err
  g = ggplot(data.frame(x = c(xmin, xmax)), aes(x = x)) + ylab("density")
  
  # Plot the null hypothesis distribution
  g = g + stat_function(fun = dnorm, geom = "line", size = 1, 
                        args = list(mean = mean_H0, sd = std_err), 
                        aes(colour = "H0"))
  
  # Plot the alternative hypothesis distribution
  g = g + stat_function(fun = dnorm, geom = "line", size = 1, 
                        args = list(mean = mean_H1, sd = std_err), 
                        aes(colour = "H1"))
  
  # Plot the critical values at level alpha
  alpha_quantile_lower = mean_H0 - qnorm(1 - alpha/2) * std_err
  alpha_quantile_upper = mean_H0 + qnorm(1 - alpha/2) * std_err
  g = g + geom_vline(xintercept=alpha_quantile_lower, size = 1)
  g = g + geom_vline(xintercept=alpha_quantile_upper, size = 1)
  
  # Plot the type I error
  x_alpha_lower <- seq(xmin, alpha_quantile_lower, 0.01)
  y_alpha_lower <- dnorm(x_alpha_lower, mean = mean_H0, sd = std_err)
  g = g + geom_area(aes(x = x, y = y), data = data.frame(x = x_alpha_lower, y = y_alpha_lower), 
                    fill = "orange1", alpha = 0.5)
  
  x_alpha_upper <- seq(alpha_quantile_upper, xmax, 0.01)
  y_alpha_upper <- dnorm(x_alpha_upper, mean = mean_H0, sd = std_err)
  g = g + geom_area(aes(x = x, y = y), data = data.frame(x = x_alpha_upper, y = y_alpha_upper), 
                    fill = "orange1", alpha = 0.5)
  
  # Plot the type II error
  # if (mean_H1 >= 0)
  #   x_beta <- seq(xmin, alpha_quantile_upper, 0.01)
  # else
  #   x_beta <- seq(alpha_quantile_lower, xmax, 0.01)
  
  x_beta <- seq(alpha_quantile_lower, alpha_quantile_upper, 0.01)
  
  y_beta <- dnorm(x_beta, mean = mean_H1, sd = std_err)
  g = g + geom_area(aes(x = x, y = y), data = data.frame(x = x_beta, y = y_beta), 
                    fill = "skyblue1", alpha = 0.5)
  
  
  # Set background colour to white
  g = g + theme_bw()
  
  # Set title
  g = g + ggtitle(label = "Type I and II Errors")
  
  # Set legend
  g = g + scale_colour_manual("Hypotheses", values = c("orange3", "skyblue3"))
  
  # Show plot
  g
}

# Manipulates the dynamic plot values via sliders
manipulate(
  two_sided(sigma, mean_H1, n, alpha),
  sigma = slider(1, 5, step = 1, initial = 3),
  mean_H1 = slider(-5, 5, step = 1, initial = 2),
  n = slider(1, 50, step = 1, initial = 12),
  alpha = slider(0.01, 0.2, step = 0.01, initial = 0.05)
)
