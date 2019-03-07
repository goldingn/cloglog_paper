# example illustrating spatial invariance of the cloglog model, but not the
# logit or probit

# clear workspace and set RNG seed
rm(list = ls())
set.seed(1)

library(RColorBrewer)

# function to simulate presence-absence at a given scale
sim_pres <- function (x,
                      model = c('logit', 'probit', 'cloglog'),
                      param = c(alpha = -3, beta = 2),
                      A = 1) {
  # grab the right link function
  model <- match.arg(model)
  link <- switch(model,
                 logit = plogis,
                 probit = pnorm,
                 cloglog = function(x) 1 - exp(-exp(x)))
  
  # make sure A is an positive integer (number of sub-quadrats to include in larger
  # quadrats)
  A <- round(A)
  stopifnot(A > 0)
  
  # calculate eta common to all models
  eta <- param[1] + param[2] * x
  
  # probability of presence under that model
  p <- link(eta)
  
  # simulate presence/absence un sub-quadrats
  y <- replicate(A, rbinom(length(x), 1, p))
  
  # return presence in any quadrat
  apply(y, 1, max)
  
}

# function to simulate from and fit a model, and return the estimated parameters
fit_model <- function(model = c('logit', 'probit', 'cloglog'),
                      n = 1000,
                      ...) {
  
  model <- match.arg(model)
  
  # simulate x and y
  x <- rnorm(n)
  y <- sim_pres(x, model = model, ...)
  
  # fit a model
  m <- glm(y ~ x, family = binomial(model))
  
  # return the coefficients
  coefs <- coef(m)
  names(coefs) <- c('alpha', 'beta')
  
  coefs
  
}

# models and areas to plot
models <- c('cloglog', 'logit', 'probit')
areas <- 1:3

# define colours for different area sizes
paired <- brewer.pal(12, 'Paired')
cols <- colorRampPalette(paired[2:1])(3)

# and for other plot elements
panel_col <- grey(0.3)
axis_col <- grey(0.4)
lab_col <- grey(0.5)
line_col <- grey(0.7)

# location of panel letter
panel_adj <- 0

# axis label size
axis_size <- 0.8

# line widths
lw1 <- 10

# titles
models_title <- c('complementary log-log', 'logit', 'probit')
areas_title <- paste('area =', areas)

# plot dimensions
xlim <- c(-4, -1.5)
ylim <- c(1.5, 3)

# vertical position of width bars for cloglog panel
wbar_heights <- c(2.9, 2.7)

# device
png('cloglog_box_2.png',
    width = 6000,
    height = 2000,
    pointsize = 80)

par(mfrow = c(1, 3))

for (m in 1:3) {
  
  # set up the plot
  plot.new()
  plot.window(xlim = xlim,
              ylim = ylim)
  
  # add axes and labels
  axis(1,
       col = axis_col)
  title(xlab = expression(alpha),
        cex.lab = axis_size * 1.5,
        col.lab = lab_col,
        line = 3)
  
  # y axis only on first panel
  if (m == 1) {
    axis(2,
         las = 2,
         col = axis_col)
    mtext(text = expression(beta),
          side = 2,
          las = 1,
          cex = axis_size,
          col = lab_col,
          line = 2)
  }
  
  # panel titles (model names)
  mtext(models_title[m],
        side = 3,
        line = 1.2,
        cex = 1.2,
        col = panel_col)
  
  # panel labels (letters)
  mtext(LETTERS[m],
        side = 3,
        cex = 1.5,
        adj = panel_adj,
        line = 1.5,
        col = panel_col)
  
  
  for (a in 3:1) {
    
    # 1000 random draws of the coefficients, fitted to a 1000-observation dataset
    est <- t(replicate(1000, fit_model(models[m], A = areas[a])))
    
    # plot the draws
    points(est,
           pch = 16,
           cex = 0.8,
           col = cols[a])
    
  }
  
  # add lines for the true parameters
  abline(v = -3,
         lwd = lw1,
         lend = 2,
         col = line_col)
  
  abline(h = 2,
         lwd = lw1,
         lend = 2,
         col = line_col)
  
  # add vertical lines and widths for cloglog offsets
  if (models[m] == 'cloglog') {
    
    # vertival lines
    offsets <- log(areas[-1])
    abline(v = -3 + offsets,
           col = line_col,
           lwd = lw1,
           lend = 2,
           lty = 3)
    
    # width bars & labels
    for (i in 1:2) {
      
      # bar
      arrows(x0 = -3,
             y0 = wbar_heights[i],
             x1 = -3 + offsets[i],
             y1 = wbar_heights[i],
             code = 3,
             angle = 90,
             lwd = 8,
             lend = 2,
             col = grey(0.4))
      
      # label
      mid <- -3 + offsets[i] / 2
      text(x = mid,
           y = wbar_heights[i],
           labels = sprintf('ln(%i)',
                            areas[i + 1]),
           pos = 3,
           offset = 0.4,
           col = grey(0.4))
      
    }
    
  }
  
}

dev.off()
