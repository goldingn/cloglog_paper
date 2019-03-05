# plotting figure 1 for the cloglog manuscript

# clear workspace
rm(list = ls())

# set up colours
library(RColorBrewer)

paired <- brewer.pal(12, 'Paired')
col_abs <- paired[1]
col_pres <- paired[2]

panel_col <- grey(0.3)
axis_col <- grey(0.4)
lab_col <- grey(0.5)
eqn_col <- grey(0.6)
line_col <- grey(0.7)
link_col <- grey(0.85)

# link plot default line widths
lw1 <- 10
lw2 <- 20

# location of panel letter
panel_adj <- -0.15

# y label position
ylab_line <- 2

# axis label size
axis_size <- 0.8

# size of numeric labels
number_size <- 1.2

# device
png('cloglog_fig_1.png',
    width = 6000,
    height = 2000,
    pointsize = 90)

par(mfrow = c(1, 3),
    mar = c(5, 6, 4, 2),
    oma = c(0, 1.5, 0, 0))

# A) Poisson density
x <- 0:8
dens <- dpois(x, 2)

# colour absence differently & plot
col_pois <- c(col_abs,
              rep(col_pres, length(x) - 1))

bp <- barplot(dens,
              border = NA,
              ylim = c(0, 0.42),
              las = 1,
              col.axis = axis_col,
              col = col_pois)

text(x = 8,
     y = 0.36,
     label = expression(paste('N ~ ',
                              italic(Poisson),
                              (lambda))),
     col = eqn_col,
     cex = 1.3,
     adj = 0.5)

text(x = 8,
     y = 0.33,
     label = expression(paste(lambda,
                              ' = 2')),
     col = eqn_col,
     cex = 1.3,
     adj = 0.5)

text(x = 8,
     y = 0.26,
     label = expression(paste(italic(p),
                              (N > 0) == 0.86)),
     col = col_pres,
     cex = 1.4,
     adj = 0.5)

# label the probability of absence
p_abs <- exp(-2)
text(x = bp[1] - 0.1,
     y = p_abs + 0.02,
     label = round(p_abs, 2),
     cex = number_size,
     col = col_pres)

# label axes etc
title(xlab = expression(N),
      cex.lab = axis_size * 1.5,
      col.lab = lab_col,
      line = 3)

mtext(text = expression(italic(p)(N)),
      side = 2,
      las = 1,
      at = 0.24,
      cex = axis_size,
      col = lab_col,
      line = ylab_line)

axis(side = 1,
     at = bp,
     labels = x,
     lwd = 0,
     col.axis = axis_col)

# add presence/absence dividing line
cutoff <- mean(bp[1:2])
abline(v = cutoff,
       lwd = lw1,
       col = line_col)

text(x = cutoff - 0.4,
     y = 0.34,
     label = 'absence',
     srt = 90,
     col = col_abs,
     cex = 1.2)

text(x = cutoff + 0.4,
     y = 0.34,
     label = 'presence',
     srt = 90,
     col = col_pres,
     cex = 1.2)

# Panel label
mtext('A',
      side = 3,
      cex = 1.5,
      adj = panel_adj,
      line = 1.5,
      col = panel_col)

# B) cloglog function

# sequence of lambda values & cll function
lambda <- seq(0.001, 5, length = 1000)
cloglog <- function (x) 1 - exp(-exp(x))

# set up plot
plot(cloglog(log(lambda)) ~ lambda,
     type = 'n',
     ylim = c(0, 1),
     xlab = '',
     ylab = '',
     las = 1,
     frame.plot = FALSE,
     col.axis = axis_col)

# label axes etc
title(xlab = expression(lambda),
      cex.lab = axis_size * 1.5,
      col.lab = lab_col,
      line = 3)

mtext(text = expression(italic(p)(italic(y) == 1)),
      side = 2,
      las = 1,
      cex = axis_size,
      col = lab_col,
      line = ylab_line)

# add a lookup line and text label
prob <- 1 - exp(-2)

lines(x = c(2, 2),
      y = c(-1, prob),
      lwd = lw1,
      col = line_col)

lines(x = c(-1, 2),
      y = c(prob, prob),
      lwd = lw1,
      col = line_col)

text(x = 0.8,
     y = prob + 0.04,
     label = round(prob, 2),
     cex = number_size,
     col = col_pres)

# the main line
lines(cloglog(log(lambda)) ~ lambda,
      lwd = lw2,
      col = col_pres)

# Panel label
mtext('B',
      side = 3,
      cex = 1.5,
      adj = panel_adj,
      line = 1.5,
      col = panel_col)

# add plot on eta scale
eta <- seq(-3, 3, length = 1000)

# set up plot
plot(cloglog(eta) ~ eta,
     type = 'n',
     ylim = c(0, 1),
     xlab = '',
     ylab = '',
     las = 1,
     frame.plot = FALSE,
     col.axis = axis_col)

lines(x = log(c(2, 2)),
      y = c(-1, prob),
      lwd = lw1,
      col = line_col)

lines(x = c(-5, log(2)),
      y = c(prob, prob),
      lwd = lw1,
      col = line_col)

text(x = -1.2,
     y = prob + 0.04,
     label = round(prob, 2),
     cex = number_size,
     col = col_pres)


# the other links
lines(plogis(eta) ~ eta,
      lwd = lw2,
      col = link_col)

lines(pnorm(eta) ~ eta,
      lwd = lw2,
      lty = 3,
      col = link_col)

# the main line
lines(cloglog(eta) ~ eta,
      lwd = lw2,
      col = col_pres)

# label axes etc
title(xlab = expression(eta == ln(lambda)),
      cex.lab = axis_size * 1.5,
      col.lab = lab_col,
      line = 3)

mtext(text = expression(italic(p)(italic(y) == 1)),
      side = 2,
      las = 1,
      cex = axis_size,
      col = lab_col,
      line = ylab_line)

# Panel label
mtext('C',
      side = 3,
      cex = 1.5,
      adj = panel_adj,
      line = 1.5,
      col = panel_col)

# close device
dev.off()



