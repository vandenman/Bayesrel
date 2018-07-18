#' plot function for a posterior sample, gives posterior and prior distribution and pie plots
#' input is the main reliability estimation object
#'
#' @export
brelPlot <- function(res, top.align = FALSE, name = NULL, blackwhite = FALSE, criteria = TRUE){
  options(warn = -1)
  estimate <- readline(prompt = "Enter the estimate you want to plot: ")
  posi <- grep(estimate, res$estimates, ignore.case = T)
  if (!is.null(name)) {estimate <- name}
  samp <- coda::as.mcmc(unlist(res$bay$samp[posi]))
  prior <- unlist(res$priors[posi])
  if (is.null(prior)) {
    return("You need to rerun the reliability estimation with ‘prior.samp = TRUE‘ for the plot")
  }
  par(cex.main = 1.5, mar = c(4, 1,  3, 1), mgp = c(3.5, 1, 0), cex.lab = 1.5,
      font.lab = 2, cex.axis = 1.8, bty = "n", las = 1)

  hdi <- coda::HPDinterval(samp)
  med <- median(samp) # main alignment variable
  rad <- .06 # radius of pie plots
  peak <- max(density(samp)$y)

  if (med > .65){
    mid <-  med/2
    pos.x <- mid - .11
    pos.y <- peak * .6
  }
  if (med < .35) {
    mid <-  1 - (1 - med)/2
    pos.x <- mid - .11
    pos.y <- peak * .6

  }
  if (med >= .35 && med <= .65) {top.align <- TRUE}
  if (top.align == TRUE) {
    pos.x <- med - .11
    pos.y <- peak * 1.35
    mid <- med #(pos.x + pos.x +.2)/2
    rad <- .04
  }

  cuts <- c(.70, .80, .90)
  main <- paste("Reliability Fit for Distribution Mass of", estimate)
  colos <- c("firebrick", "cadetblue3", "cornflowerblue","navy")
  if (blackwhite){
    colos <- c("gray100", "gray80", "gray50", "gray8")
  }

  dens.prior <- density(prior, from = 0, to = 1, n = 2e3)
  xx0 <- min(which(dens.prior$x <= cuts[1]))
  xx1 <- max(which(dens.prior$x <= cuts[1]))
  xx2 <- max(which(dens.prior$x <= cuts[2]))
  xx3 <- max(which(dens.prior$x <= cuts[3]))
  xx4 <- max(which(dens.prior$x <= 1))

  if (!is.integer(xx0)) xx0 <- 1
  if (!is.integer(xx1)) xx1 <- 1
  if (!is.integer(xx2)) xx2 <- 1
  if (!is.integer(xx3)) xx3 <- 1
  if (!is.integer(xx4)) xx4 <- 1

  dens.post <- density(samp, adjust = 1.75, n = 2e3)
  x0 <- min(which(dens.post$x <= cuts[1]))
  x1 <- max(which(dens.post$x <= cuts[1]))
  x2 <- max(which(dens.post$x <= cuts[2]))
  x3 <- max(which(dens.post$x <= cuts[3]))
  x4 <- max(which(dens.post$x <= 1))

  if (!is.integer(x0)) x0 <- 1
  if (!is.integer(x1)) x1 <- 1
  if (!is.integer(x2)) x2 <- 1
  if (!is.integer(x3)) x3 <- 1
  if (!is.integer(x4)) x4 <- 1

  y1 <- sum(prior <= cuts[1])
  y2 <- sum(prior <= cuts[2])
  y3 <- sum(prior <= cuts[3])
  y4 <- sum(prior <= 1)

  z1 <- sum(samp <= cuts[1])
  z2 <- sum(samp <= cuts[2])
  z3 <- sum(samp <= cuts[3])
  z4 <- sum(samp <= 1)

  pie.prior <- c(y1, y2-y1, y3-y2, y4-y3)
  pie.prior[pie.prior == 0] <- 1e-20
  pie.prior.labels <- as.character(round(pie.prior/(length(prior)*1e-2), 1))
  pie.post <- c(z1, z2-z1, z3-z2, z4-z3)
  pie.post[pie.post == 0] <- 1e-20
  pie.post.labels <- as.character(round(pie.post/(length(samp)*1e-2), 1))
  for (i in 1:4){
    if (as.numeric(pie.prior.labels[i]) > 3) {
      pie.prior.labels[i] <- paste(pie.prior.labels[i], "%")
    } else{
      pie.prior.labels[i] <- ""
    }
    if (as.numeric(pie.post.labels[i]) > 3) {
      pie.post.labels[i] <- paste(pie.post.labels[i], "%")
    } else{
      pie.post.labels[i] <- ""
    }
  }

  # ------------------------------- plotting --------------------------------------

  if (top.align && criteria){
    plot(density(samp, adjust = 1.75), type = "l", axes = F, xlab = NA, ylab = NA,
         xlim = c(0, 1), ylim = c(0,  peak * 1.65),
         lwd = 3, main = "")
    plotShadePrior(dens.prior, xx = c(xx0, xx1, xx2, xx3, xx4), cols = colos, cutoffs = criteria, grey = blackwhite)
    plotShadePost(dens.post, xx = c(x0, x1, x2, x3, x4), cols = colos, cutoffs = criteria, grey = blackwhite)

    lines(density(prior, from = 0, to = 1), lty = 5, lwd = 3)

    axis(side = 1, at = seq(0, 1, by = .2), labels = seq(0, 1, by = .2), cex.axis = 1.25, lwd = 1.5)


    arrows(x0 = hdi[1], y0 = peak, x1 = hdi[2], y1 = peak, angle = 90, length = 0.05,
           code = 3, lwd = 2)
    text(paste("95% HDI: [", round(hdi[1], 3), ", ", round(hdi[2], 3),"]", sep =""), x = sum(hdi)/2, y = peak*1.04, adj = .5, cex = 1.25)
    text(paste("median = ", round(med, 3), sep =""), x = sum(hdi)/2, y = peak*1.09, adj = .5, cex = 1.25)
    title(main, line = -1, cex.main = 2)

    text("poor", x = cuts[1]/2, y = peak*-.025, adj = 0.5, cex = 1.25)
    text("acceptable", x = (cuts[1] + cuts[2])/2, y = peak*-.025, adj = .5, cex = 1.25)
    text("preferable", x = (cuts[2] + cuts[3])/2, y = peak*-.025, adj = .5, cex = 1.25)
    text("desirable?", x = (cuts[3] + 1)/2, y = peak*-.025, adj = .5, cex = 1.25)

    f.prior <- plotrix::floating.pie(xpos = pos.x, ypos = pos.y, x = pie.prior, radius = rad,
                            col = colos, startpos = 0)
    plotrix::pie.labels(pos.x, pos.y, radius = rad + .007, labels = pie.prior.labels, f.prior)
    f.post <- plotrix::floating.pie(xpos = pos.x +.22, ypos = pos.y, x = pie.post, radius = rad,
                           col = colos, startpos = 0)
    plotrix::pie.labels(pos.x +.22, pos.y, radius = rad + .007, labels = pie.post.labels, f.post)
    text("prior", x = pos.x, y = pos.y * 1.12, cex = 1.5)
    text("posterior", x = pos.x + 0.22, y = pos.y * 1.12, cex = 1.5)
    segments(x0 = pos.x -.03, y0 = pos.y * 1.1, x1 = pos.x +.03, y1 = pos.y * 1.1, lwd = 1.75, lty = 5)
    segments(x0 = pos.x + .172, y0 = pos.y * 1.1, x1 = pos.x + .268, y1 = pos.y * 1.1, lwd = 1.75, lty = 1)
    text("Distribution", x = mid, y = pos.y * 1.18, cex = 1.75)
    legend(x = mid, y = pos.y *.925, inset=.02, xjust = 0.5, fill=colos, horiz=TRUE, cex=1.25, bty = "n",
           c("poor","acceptable","preferable", "desirable?"))


  } else {

    plot(density(samp, adjust = 1.75), type = "l", axes = F, xlab = NA, ylab = NA,
         xlim = c(0, 1), ylim = c(0,  peak * 1.2),
         lwd = 3, main = "")
    plotShadePrior(dens.prior, xx = c(xx0, xx1, xx2, xx3, xx4), cols = colos, cutoffs = criteria, grey = blackwhite)
    plotShadePost(dens.post, xx = c(x0, x1, x2, x3, x4), cols = colos, cutoffs = criteria, grey = blackwhite)

    lines(density(prior, from = 0, to = 1), lty = 5, lwd = 3)

    axis(side = 1, at = seq(0, 1, by = .2), labels = seq(0, 1, by = .2), cex.axis = 1.25, lwd = 1.5)


    arrows(x0 = hdi[1], y0 = peak, x1 = hdi[2], y1 = peak, angle = 90, length = 0.05,
           code = 3, lwd = 2)
    text(paste("95% HDI: [", round(hdi[1], 3), ", ", round(hdi[2], 3),"]", sep =""), x = sum(hdi)/2, y = peak*1.04, adj = .5, cex = 1.25)
    text(paste("median = ", round(med, 3), sep =""), x = sum(hdi)/2, y = peak*1.09, adj = .5, cex = 1.25)
    title(main, line = -1, cex.main = 2)

    if (criteria){
      text("poor", x = cuts[1]/2, y = peak*-.025, adj = 0.5, cex = 1.25)
      text("acceptable", x = (cuts[1] + cuts[2])/2, y = peak*-.025, adj = .5, cex = 1.25)
      text("preferable", x = (cuts[2] + cuts[3])/2, y = peak*-.025, adj = .5, cex = 1.25)
      text("desirable?", x = (cuts[3] + 1)/2, y = peak*-.025, adj = .5, cex = 1.25)

      f.prior <- plotrix::floating.pie(xpos = pos.x, ypos = pos.y, x = pie.prior, radius = rad,
                                       col = colos, startpos = 0)
      plotrix::pie.labels(pos.x, pos.y, radius = rad + .007, labels = pie.prior.labels, f.prior)
      f.post <- plotrix::floating.pie(xpos = pos.x +.22, ypos = pos.y, x = pie.post, radius = rad,
                                      col = colos, startpos = 0)
      plotrix::pie.labels(pos.x +.22, pos.y, radius = 0.067, labels = pie.post.labels, f.post)
      text("prior", x = pos.x, y = pos.y * 1.29, cex = 1.5)
      text("posterior", x = pos.x + 0.22, y = pos.y * 1.29, cex = 1.5)
      segments(x0 = pos.x -.03, y0 = pos.y * 1.245, x1 = pos.x +.03, y1 = pos.y * 1.245, lwd = 1.75, lty = 5)
      segments(x0 = pos.x + .172, y0 = pos.y * 1.245, x1 = pos.x + .268, y1 = pos.y * 1.245, lwd = 1.75, lty = 1)
      text("Distribution", x = mid, y = pos.y * 1.41, cex = 1.75)
      legend(x = mid, y = pos.y * 0.77, inset=.02, xjust = 0.5, fill=colos, horiz=TRUE, cex=1.25, bty = "n",
             c("poor","acceptable","preferable", "desirable?"))
    } else {
      legend(x = mid, y = pos.y*0.7, c("posterior", "prior"), lty = c(1, 5), lwd = 2, title = "Distribution", cex = 1.25, bty = "n")
    }

  }
  options(warn = 0)
}


plotShadePost <- function(dens, xx, cols, cutoffs, grey){
  if (cutoffs){
    color_transp <- adjustcolor(cols, alpha.f = .7)
    if (grey) {color_transp <- adjustcolor(cols, alpha.f = .8)}
    with(dens, polygon(x[c(xx[1],xx[1]:xx[2],xx[2])], c(0, y[xx[1]:xx[2]], 0), col = color_transp[1]))
    with(dens, polygon(x[c(xx[2],xx[2]:xx[3],xx[3])], c(0, y[xx[2]:xx[3]], 0), col = color_transp[2]))
    with(dens, polygon(x[c(xx[3],xx[3]:xx[4],xx[4])], c(0, y[xx[3]:xx[4]], 0), col = color_transp[3]))
    with(dens, polygon(x[c(xx[4],xx[4]:xx[5],xx[5])], c(0, y[xx[4]:xx[5]], 0), col = color_transp[4]))
  }
  else {
    color_transp <- adjustcolor(cols[3], alpha.f = .7)
    if (grey) {color_transp <- adjustcolor(cols[3], alpha.f = .8)}
    with(dens, polygon(x[c(xx[1],xx[1]:xx[5],xx[5])], c(0, y[xx[1]:xx[5]], 0), col = color_transp))
  }

}

plotShadePrior <- function(dens, xx, cols, cutoffs, grey){
  if (cutoffs){
    color_transp <- adjustcolor(cols, alpha.f = .5)
    if (grey){color_transp <- adjustcolor(cols, alpha.f = .7)}
    with(dens, polygon(x[c(xx[1],xx[1]:xx[2],xx[2])], c(0, y[xx[1]:xx[2]], 0), col = color_transp[1]))
    with(dens, polygon(x[c(xx[2],xx[2]:xx[3],xx[3])], c(0, y[xx[2]:xx[3]], 0), col = color_transp[2]))
    with(dens, polygon(x[c(xx[3],xx[3]:xx[4],xx[4])], c(0, y[xx[3]:xx[4]], 0), col = color_transp[3]))
    with(dens, polygon(x[c(xx[4],xx[4]:xx[5],xx[5])], c(0, y[xx[4]:xx[5]], 0), col = color_transp[4]))
  }
  else {
    color_transp <- adjustcolor("grey", alpha.f = .5)
    if (grey) {color_transp <- adjustcolor("grey75", alpha.f = .7)}
    with(dens, polygon(x[c(xx[1],xx[1]:xx[5],xx[5])], c(0, y[xx[1]:xx[5]], 0), col = color_transp))
  }
}

