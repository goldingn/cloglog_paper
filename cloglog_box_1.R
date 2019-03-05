# example predicting abundance from presence-absence data generated under a
# Poisson model

# clear workspace and set RNG seed
rm(list = ls())
set.seed(1)

library(RColorBrewer)

# colours
paired <- brewer.pal(12, 'Paired')
col_abs <- paired[1]
col_pres <- paired[2]
line_col <- grey(0.7)
lab_col <- grey(0.5)
link_col <- grey(0.85)
axis_col <- grey(0.4)


# generate fake abundance data
n <- 1000
coef <- c(-1, 2)
df <- data.frame(x = runif(n))
eta <- coef[1] + coef[2] * df$x
lambda <- exp(eta)

sample_presence <- function (lambda) {
  # given vector of lambdas, sample a random presence vector by sampling 
  # abundances then truncating. Ensure there are at least two 1s and two 0s
  OK <- FALSE
  while (!OK) {
    abundance <- rpois(length(lambda), lambda)
    presence <- pmin(abundance, 1)
    OK <- sum(presence) > 2 & sum(1 - presence) > 2
  }
  presence
}

bernoulli_coef <- function (presence, df, link = c('logit', 'probit', 'cloglog')) {
  # given a vector of presence-absence data, and corresponding dataframe of 
  # covariates, return the coefficients from a Bernoulli glm with the stated
  # link function
  link <- match.arg(link)
  glm <- glm(presence ~ ., data = df, family = binomial(link))
  coef(glm)
}

presences <- replicate(1000,
                       sample_presence(lambda),
                       simplify = FALSE)

cloglog_coefs <- sapply(presences, bernoulli_coef, df, 'cloglog')
logit_coefs <- sapply(presences, bernoulli_coef, df, 'logit')
probit_coefs <- sapply(presences, bernoulli_coef, df, 'probit')


# plot figure

# device
png('cloglog_box_1.png',
    width = 2000,
    height = 2000,
    pointsize = 50)

par(mfrow = c(5, 5, 1, 1))

xlim <- c(-1.4, -0.2)
ylim <- c(1, 4.5)

# point & line sizes
cex <- 0.5
lw1 <- 5

plot(t(cloglog_coefs),
     type = 'n',
     xlim = xlim,
     ylim = ylim,
     xlab = '',
     ylab = '',
     axes = FALSE)

axis(1,
     col.axis = axis_col)

axis(2,
     col.axis = axis_col,
     las = 2)

title(xlab = expression(alpha),
      cex.lab = 1.5,
      col.lab = lab_col)

mtext(side = 2,
      text = expression(beta),
      col = lab_col,
      las = 1,
      line = 3,
      cex = 1.5)

lines(x = rep(coef[1], 2),
      y = ylim + c(-10, 0),
      col = line_col,
      lwd = lw1)

lines(x = xlim + c(-10, 0),
      y = rep(coef[2], 2),
      col = line_col,
      lwd = lw1)

points(t(logit_coefs),
       pch = 16,
       cex = cex,
       col = link_col)

points(t(probit_coefs),
       pch = 17,
       cex = cex,
       col = link_col)

points(t(cloglog_coefs),
     pch = 21,
     cex = 0.6,
     bg = col_pres,
     col = grey(0.4),
     lwd = 1)

dev.off()
