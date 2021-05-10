##########################################################
# ONE-SIDED TEST - INTERACTIVE TYPE I AND II ERROR PLOT
##########################################################


library(manipulate)
library(ggplot2)


# Dynamic ggplot
one_sided <- function(sigma, mean_H1, n, alpha) {
  
  mean.H0 <- 0
  std_err <- sigma / sqrt(n)
  
  # Create the ggplot
  xmax <- mean_H1 + 4 * std_err
  xmin <- -3
  g = ggplot(data.frame(x = c(xmin, xmax)), aes(x = x)) + ylab("density")
  
  # Plot the null hypothesis distribution
  g = g + stat_function(fun = dnorm, geom = "line", size = 1, 
                        args = list(mean = mean.H0, sd = std_err), 
                        aes(colour = "H0"))
  
  # Plot the alternative hypothesis distribution
  g = g + stat_function(fun = dnorm, geom = "line", size = 1, 
                        args = list(mean = mean_H1, sd = std_err), 
                        aes(colour = "H1"))
  
  # Plot the critical value at level alpha
  alpha_quantile = mean.H0 + qnorm(1 - alpha) * std_err
  g = g + geom_vline(xintercept=alpha_quantile, size = 1)
  
  # Plot the type I error
  x_alpha <- seq(alpha_quantile, xmax, 0.01)
  y_alpha <- dnorm(x_alpha, mean = mean.H0, sd = std_err)
  g = g + geom_area(aes(x = x, y = y), data = data.frame(x = x_alpha, y = y_alpha), 
                    fill = "orange1", alpha = 0.5)
  
  # Plot the type II error
  x_beta <- seq(xmin, alpha_quantile, 0.01)
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
  one_sided(sigma, mean_H1, n, alpha),
  sigma = slider(1, 5, step = 1, initial = 3),
  mean_H1 = slider(1, 5, step = 1, initial = 2),
  n = slider(1, 50, step = 1, initial = 16),
  alpha = slider(0.01, 0.2, step = 0.01, initial = 0.05)
)
