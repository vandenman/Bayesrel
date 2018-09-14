#' plot function for a posterior sample, gives posterior and prior distribution and pie plots
#' input is the main reliability estimation object and the estimate to be plotted
#'
#' @export
plotrel <- function(x, estimate, top.align = FALSE, greek = FALSE, blackwhite = FALSE, criteria = TRUE, cuts = c(.70, .80)){
  options(warn = -1)
  posi <- grep(estimate, x$estimates, ignore.case = T)
  samp <- coda::as.mcmc(unlist(x$bay$samp[posi]))
  prior <- unlist(x$priors[posi])
  if (is.null(prior)) {
    return("You need to rerun the reliability estimation with ‘prior.samp = TRUE‘ for the plot")
  }
  par(cex.main = 1.5, mar = c(4, 1,  3, 1), mgp = c(3.5, 1, 0), cex.lab = 1.5,
      font.lab = 2, cex.axis = 1.8, bty = "n", las = 1)

  hdi <- coda::HPDinterval(samp)
  med <- median(samp) # main alignment variable
  rad <- .05 # radius of pie plots
  peak <- max(density(samp)$y)

  if (med > .65){
    mid <-  med/2
    pos.x <- mid - .1
    pos.y <- peak * .6
  }
  if (med < .35) {
    mid <-  1 - (1 - med)/2
    pos.x <- mid - .1
    pos.y <- peak * .6

  }
  if (med >= .35 && med <= .65) {top.align <- TRUE}
  if (!criteria) {top.align <- F}
  if (top.align == TRUE) {
    pos.x <- med - .1
    pos.y <- peak * 1.35
    mid <- med #(pos.x + pos.x +.2)/2
    rad <- .04
  }

  main <- paste("Reliability Fit for Distribution Mass of", estimate)
  if (greek) {
    if (estimate == "lambda2"){
      main <- substitute(paste("Posterior Plot for ", lambda [2]))
    }
    if (estimate == "lambda4"){
      main <- substitute(paste("Posterior Plot for ", lambda [4]))
    }
    if (estimate == "lambda6"){
      main <- substitute(paste("Posterior Plot for ", lambda [6]))
    }
    if (estimate == "alpha" || estimate == "omega" || estimate == "glb"){
      estimate <- as.symbol(estimate)
      main <- substitute(paste("Posterior Plot for ", estimate), list(estimate = estimate))
    }
  }
  colos <- c("firebrick", "cornflowerblue","navy")
  if (blackwhite){
    colos <- c("gray100", "gray70", "gray10")
  }

  dens.prior <- density(prior, from = 0, to = 1, n = 2e3)
  xx0 <- min(which(dens.prior$x <= cuts[1]))
  xx1 <- max(which(dens.prior$x <= cuts[1]))
  xx2 <- max(which(dens.prior$x <= cuts[2]))
  xx3 <- max(which(dens.prior$x <= 1))

  if (!is.integer(xx0)) xx0 <- 1
  if (!is.integer(xx1)) xx1 <- 1
  if (!is.integer(xx2)) xx2 <- 1
  if (!is.integer(xx3)) xx3 <- 1

  dens.post <- density(samp, adjust = 1.75, n = 2e3)
  x0 <- min(which(dens.post$x <= cuts[1]))
  x1 <- max(which(dens.post$x <= cuts[1]))
  x2 <- max(which(dens.post$x <= cuts[2]))
  x3 <- max(which(dens.post$x <= 1))

  if (!is.integer(x0)) x0 <- 1
  if (!is.integer(x1)) x1 <- 1
  if (!is.integer(x2)) x2 <- 1
  if (!is.integer(x3)) x3 <- 1

  y1 <- sum(prior <= cuts[1])
  y2 <- sum(prior <= cuts[2])
  y3 <- sum(prior <= 1)

  z1 <- sum(samp <= cuts[1])
  z2 <- sum(samp <= cuts[2])
  z3 <- sum(samp <= 1)

  pie.prior <- c(y1, y2-y1, y3-y2)
  pie.prior[pie.prior == 0] <- 1e-20
  pie.prior.labels <- as.character(round(pie.prior/(length(prior)*1e-2), 1))
  pie.post <- c(z1, z2-z1, z3-z2)
  pie.post[pie.post == 0] <- 1e-20
  pie.post.labels <- as.character(round(pie.post/(length(samp)*1e-2), 1))
  for (i in 1:3){
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
    plotShadePrior(dens.prior, xx = c(xx0, xx1, xx2, xx3), cols = colos, criteria = criteria, blackwhite = blackwhite)
    plotShadePost(dens.post, xx = c(x0, x1, x2, x3), cols = colos, criteria = criteria, blackwhite = blackwhite)

    lines(density(prior, from = 0, to = 1), lty = 5, lwd = 3)

    axis(side = 1, at = seq(0, 1, by = .2), labels = seq(0, 1, by = .2), cex.axis = 1.25, lwd = 1.5)


    arrows(x0 = hdi[1], y0 = peak, x1 = hdi[2], y1 = peak, angle = 90, length = 0.05,
           code = 3, lwd = 2)
    text(paste("95% HDI: [", round(hdi[1], 3), ", ", round(hdi[2], 3),"]", sep =""), x = sum(hdi)/2, y = peak*1.04, adj = .5, cex = 1.25)
    text(paste("median = ", round(med, 3), sep =""), x = sum(hdi)/2, y = peak*1.09, adj = .5, cex = 1.25)
    title(main, line = -1, cex.main = 2)

    text("insufficient", x = cuts[1]/2, y = peak*-.025, adj = 0.5, cex = 1.25)
    text("sufficient", x = (cuts[1] + cuts[2])/2, y = peak*-.025, adj = .5, cex = 1.25)
    text("good", x = (cuts[2] + 1)/2, y = peak*-.025, adj = .5, cex = 1.25)

    f.prior <- plotrix::floating.pie(xpos = pos.x, ypos = pos.y, x = pie.prior, radius = rad,
                            col = colos, startpos = 0)
    plotrix::pie.labels(pos.x, pos.y, radius = rad + .006, labels = pie.prior.labels, f.prior)
    f.post <- plotrix::floating.pie(xpos = pos.x +.2, ypos = pos.y, x = pie.post, radius = rad,
                           col = colos, startpos = 0)
    plotrix::pie.labels(pos.x +.2, pos.y, radius = rad + .006, labels = pie.post.labels, f.post)
    text("prior", x = pos.x, y = pos.y * 1.125, cex = 1.5)
    text("posterior", x = pos.x + 0.2, y = pos.y * 1.125, cex = 1.5)
    segments(x0 = pos.x -.03, y0 = pos.y * 1.105, x1 = pos.x +.03, y1 = pos.y * 1.105, lwd = 1.75, lty = 5)
    segments(x0 = pos.x + .15, y0 = pos.y * 1.105, x1 = pos.x + .25, y1 = pos.y * 1.105, lwd = 1.75, lty = 1)
    #text("Distribution", x = mid, y = pos.y * 1.18, cex = 1.75)
    legend(x = mid, y = pos.y *.92, inset=.02, xjust = 0.5, fill=colos, horiz=TRUE, cex=1.25, bty = "n",
           c("insufficicent", "sufficient", "good"))


  } else {

    plot(density(samp, adjust = 1.75), type = "l", axes = F, xlab = NA, ylab = NA,
         xlim = c(0, 1), ylim = c(0,  peak * 1.2),
         lwd = 3, main = "")
    plotShadePrior(dens.prior, xx = c(xx0, xx1, xx2, xx3), cols = colos, criteria = criteria, blackwhite = blackwhite)
    plotShadePost(dens.post, xx = c(x0, x1, x2, x3), cols = colos, criteria = criteria, blackwhite = blackwhite)

    lines(density(prior, from = 0, to = 1), lty = 5, lwd = 3)

    axis(side = 1, at = seq(0, 1, by = .2), labels = seq(0, 1, by = .2), cex.axis = 1.25, lwd = 1.5)


    arrows(x0 = hdi[1], y0 = peak, x1 = hdi[2], y1 = peak, angle = 90, length = 0.05,
           code = 3, lwd = 2)
    text(paste("95% HDI: [", round(hdi[1], 3), ", ", round(hdi[2], 3),"]", sep =""), x = sum(hdi)/2, y = peak*1.04, adj = .5, cex = 1.25)
    text(paste("median = ", round(med, 3), sep =""), x = sum(hdi)/2, y = peak*1.09, adj = .5, cex = 1.25)
    title(main, line = -1, cex.main = 2)

    if (criteria){
      text("insufficient", x = cuts[1]/2, y = peak*-.025, adj = 0.5, cex = 1.25)
      text("sufficient", x = (cuts[1] + cuts[2])/2, y = peak*-.025, adj = .5, cex = 1.25)
      text("good", x = (cuts[2] + 1)/2, y = peak*-.025, adj = .5, cex = 1.25)

      f.prior <- plotrix::floating.pie(xpos = pos.x, ypos = pos.y, x = pie.prior, radius = rad,
                                       col = colos, startpos = 0)
      plotrix::pie.labels(pos.x, pos.y, radius = rad + .006, labels = pie.prior.labels, f.prior)
      f.post <- plotrix::floating.pie(xpos = pos.x +.2, ypos = pos.y, x = pie.post, radius = rad,
                                      col = colos, startpos = 0)
      plotrix::pie.labels(pos.x +.2, pos.y, radius = rad + .006, labels = pie.post.labels, f.post)
      text("prior", x = pos.x, y = pos.y * 1.27, cex = 1.5)
      text("posterior", x = pos.x + 0.2, y = pos.y * 1.27, cex = 1.5)
      segments(x0 = pos.x -.03, y0 = pos.y * 1.23, x1 = pos.x +.03, y1 = pos.y * 1.23, lwd = 1.75, lty = 5)
      segments(x0 = pos.x + .15, y0 = pos.y * 1.23, x1 = pos.x + .25, y1 = pos.y * 1.23, lwd = 1.75, lty = 1)
      #text("Distribution", x = mid, y = pos.y * 1.41, cex = 1.75)
      legend(x = mid, y = pos.y * 0.77, inset=.02, xjust = 0.5, fill=colos, horiz=TRUE, cex=1.25, # bty = "n",
             c("insufficient","sufficient", "good"))
    } else {
      legend(x = mid, y = pos.y*0.7, c("posterior", "prior"), lty = c(1, 5), lwd = 2, title = "Distribution", cex = 1.25, bty = "n")
    }

  }
  options(warn = 0)
}


plotShadePost <- function(dens, xx, cols, criteria, blackwhite){
  if (criteria){
    color_transp <- adjustcolor(cols, alpha.f = .7)
    if (blackwhite) {color_transp <- adjustcolor(cols, alpha.f = .8)}
    with(dens, polygon(x[c(xx[1],xx[1]:xx[2],xx[2])], c(0, y[xx[1]:xx[2]], 0), col = color_transp[1]))
    with(dens, polygon(x[c(xx[2],xx[2]:xx[3],xx[3])], c(0, y[xx[2]:xx[3]], 0), col = color_transp[2]))
    with(dens, polygon(x[c(xx[3],xx[3]:xx[4],xx[4])], c(0, y[xx[3]:xx[4]], 0), col = color_transp[3]))
  }
  else {
    color_transp <- adjustcolor(cols[3], alpha.f = .7)
    if (blackwhite) {color_transp <- adjustcolor(cols[3], alpha.f = .8)}
    with(dens, polygon(x[c(xx[1],xx[1]:xx[4],xx[4])], c(0, y[xx[1]:xx[4]], 0), col = color_transp))
  }

}

plotShadePrior <- function(dens, xx, cols, criteria, blackwhite){
  if (criteria){
    color_transp <- adjustcolor(cols, alpha.f = .5)
    if (blackwhite){color_transp <- adjustcolor(cols, alpha.f = .7)}
    with(dens, polygon(x[c(xx[1],xx[1]:xx[2],xx[2])], c(0, y[xx[1]:xx[2]], 0), col = color_transp[1]))
    with(dens, polygon(x[c(xx[2],xx[2]:xx[3],xx[3])], c(0, y[xx[2]:xx[3]], 0), col = color_transp[2]))
    with(dens, polygon(x[c(xx[3],xx[3]:xx[4],xx[4])], c(0, y[xx[3]:xx[4]], 0), col = color_transp[3]))
  }
  else {
    color_transp <- adjustcolor(cols[2], alpha.f = .5)
    if (blackwhite) {color_transp <- adjustcolor(cols[2], alpha.f = .7)}
    with(dens, polygon(x[c(xx[1],xx[1]:xx[4],xx[4])], c(0, y[xx[1]:xx[4]], 0), col = color_transp))
  }
}


#' plotting function of the posterior for the if-item-dropped statistics
#'
#' @export
plotifitem_one <- function(x, estimate, item.pos, criteria = TRUE, blackwhite = FALSE, top.align = FALSE, greek = FALSE, cuts = c(.70, .80)){
  options(warn = -1)
  posi <- grep(estimate, x$estimates, ignore.case = T)

  samp <- coda::as.mcmc(unlist(x$bay$samp[posi]))
  ifitems <- as.matrix(as.data.frame(x$bay$ifitem$samp[posi]))
  if(dim(ifitems) == 0) {
    return("You need to rerun the reliability estimation with ‘if.item.dropped = TRUE‘ for the plot")
  }

  #item.pos <- as.numeric(readline(prompt = "Which item should hypothetically be dropped: "))
  item <- coda::as.mcmc(ifitems[item.pos, ])

  par(cex.main = 1.5, mar = c(4, 1,  3, 1), mgp = c(3.5, 1, 0), cex.lab = 1.5,
      font.lab = 2, cex.axis = 1.8, bty = "n", las = 1)

  hdi <- coda::HPDinterval(samp)
  med <- median(samp) # main alignment variable
  rad <- .05 # radius of pie plots
  peak <- max(density(samp)$y)
  peak2 <- max(density(item)$y)
  hdi2 <- coda::HPDinterval(item)
  med2 <- median(item)
  midd <- min(med, med2)

  if (midd > .65){
    mid <-  midd/2
    pos.x <- mid - .1
    pos.y <- peak * .6
  }
  if (midd < .35) {
    mid <-  1 - (1 - midd)/2
    pos.x <- mid - .1
    pos.y <- peak * .6

  }
  if (midd >= .35 && midd <= .65) {top.align <- TRUE}
  if (!criteria) {top.align <- F}

  if (top.align == TRUE) {
    pos.x <- (med + med2)/2 -.1
    pos.y <- peak * 1.35
    mid <- (med + med2)/2
    rad <- .04
  }

  main <- paste("Posterior Plot for", estimate, "and If - Item =", item.pos, " - Dropped")
  if (greek) {
    if (estimate == "lambda2"){
      main <- substitute(paste("Posterior Plot for ", lambda [2], " and If - Item = ", item.pos, " - Dropped"), list(item.pos = item.pos))
    }
    if (estimate == "lambda4"){
      main <- substitute(paste("Posterior Plot for ", lambda [4], " and If - Item = ", item.pos, " - Dropped"), list(item.pos = item.pos))
    }
    if (estimate == "lambda6"){
      main <- substitute(paste("Posterior Plot for ", lambda [6], " and If - Item = ", item.pos, " - Dropped"), list(item.pos = item.pos))
    }
    if (estimate == "alpha" || estimate == "omega" || estimate == "glb"){
      estimate <- as.symbol(estimate)
      main <- substitute(paste("Posterior Plot for ", estimate, " and If - Item = ", item.pos, " - Dropped"), list(estimate = estimate, item.pos = item.pos))
    }
  }
  colos <- c("firebrick", "cornflowerblue","navy")
  if (blackwhite){
    colos <- c("gray100", "gray70", "gray10")
  }

  dens.post <- density(samp, adjust = 1.75, n = 2e3)
  x0 <- min(which(dens.post$x <= cuts[1]))
  x1 <- max(which(dens.post$x <= cuts[1]))
  x2 <- max(which(dens.post$x <= cuts[2]))
  x3 <- max(which(dens.post$x <= 1))

  if (!is.integer(x0)) x0 <- 1
  if (!is.integer(x1)) x1 <- 1
  if (!is.integer(x2)) x2 <- 1
  if (!is.integer(x3)) x3 <- 1

  z1 <- sum(samp <= cuts[1])
  z2 <- sum(samp <= cuts[2])
  z3 <- sum(samp <= 1)

  dens.item <- density(item, adjust = 1.75, n = 2e3)
  xx0 <- min(which(dens.item$x <= cuts[1]))
  xx1 <- max(which(dens.item$x <= cuts[1]))
  xx2 <- max(which(dens.item$x <= cuts[2]))
  xx3 <- max(which(dens.item$x <= 1))

  if (!is.integer(xx0)) xx0 <- 1
  if (!is.integer(xx1)) xx1 <- 1
  if (!is.integer(xx2)) xx2 <- 1
  if (!is.integer(xx3)) xx3 <- 1

  y1 <- sum(item <= cuts[1])
  y2 <- sum(item <= cuts[2])
  y3 <- sum(item <= 1)


  pie.item <- c(y1, y2-y1, y3-y2)
  pie.item[pie.item == 0] <- 1e-20
  pie.item.labels <- as.character(round(pie.item/(length(item)*1e-2), 1))
  pie.post <- c(z1, z2-z1, z3-z2)
  pie.post[pie.post == 0] <- 1e-20
  pie.post.labels <- as.character(round(pie.post/(length(samp)*1e-2), 1))

  for (i in 1:3){
    if (as.numeric(pie.item.labels[i]) > 3) {
      pie.item.labels[i] <- paste(pie.item.labels[i], "%")
    } else{
      pie.item.labels[i] <- ""
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
         xlim = c(0, 1), ylim = c(0,  peak * 1.64),
         lwd = 3, main = "")
    plotShadePost(dens.post, xx = c(x0, x1, x2, x3),
                              cols = colos, criteria = criteria, blackwhite = blackwhite)
    plotShadePrior(dens.item, xx = c(xx0, xx1, xx2, xx3),
                             cols = colos, criteria = criteria, blackwhite = blackwhite)

    lines(density(item, adjust = 1.75), lty = 4, lwd = 3)

    axis(side = 1, at = seq(0, 1, by = .2), labels = seq(0, 1, by = .2), cex.axis = 1.25, lwd = 1.5)

    text("insufficient", x = cuts[1]/2, y = peak*-.025, adj = 0.5, cex = 1.25)
    text("sufficient", x = (cuts[1] + cuts[2])/2, y = peak*-.025, adj = .5, cex = 1.25)
    text("good", x = (cuts[2] + 1)/2, y = peak*-.025, adj = .5, cex = 1.25)

    arrows(x0 = hdi[1], y0 = peak, x1 = hdi[2], y1 = peak, angle = 90, length = 0.05,
           code = 3, lwd = 1.5)
    text(substitute(paste("95% HDI"[original]: "[", cc1, ", ", cc2,"]"), list(cc1 = round(hdi[1], 3), cc2 = round(hdi[2], 3))),
         x = pos.x/3, y = peak*1.04, adj = 0)
    text(labels = substitute(paste("median" [original], " = ", cc), list(cc = round(med, 3))),
         x = pos.x/3, y = peak*1.09, adj = 0)

    arrows(x0 = hdi2[1], y0 = peak2, x1 = hdi2[2], y1 = peak2, angle = 90, length = 0.05,
           code = 3, lwd = 1.5, lty = 4)
    text(substitute(paste("95% HDI"[ifitem]: "[", cc1, ", ", cc2,"]"), list(cc1 = round(hdi2[1], 3), cc2 = round(hdi2[2], 3))),
         x = pos.x/3, y = peak*0.94, adj = 0)
    text(labels = substitute(paste("median" [ifitem], " = ", cc), list(cc = round(med2, 3))),
         x = pos.x/3, y = peak*.99, adj = 0)
    title(main, line = -1, cex.main = 1.75)

    f.item <- plotrix::floating.pie(xpos = pos.x, ypos = pos.y, x = pie.item, radius = rad,
                                    col = colos, startpos = 0)
    plotrix::pie.labels(pos.x, pos.y, radius = rad + .006, labels = pie.item.labels, f.item)
    f.post <- plotrix::floating.pie(xpos = pos.x +.2, ypos = pos.y, x = pie.post, radius = rad,
                                    col = colos, startpos = 0)
    plotrix::pie.labels(pos.x +.2, pos.y, radius = rad + .006, labels = pie.post.labels, f.post)
    text("if-item-dropped", x = pos.x, y = pos.y * 1.125, cex = 1.5)
    text("original", x = pos.x + 0.2, y = pos.y * 1.125, cex = 1.5)
    segments(x0 = pos.x -.075, y0 = pos.y * 1.105, x1 = pos.x +.075, y1 = pos.y * 1.105, lwd = 1.75, lty = 4)
    segments(x0 = pos.x + .16, y0 = pos.y * 1.105, x1 = pos.x + .24, y1 = pos.y * 1.105, lwd = 1.75, lty = 1)
    #text("Distribution", x = mid, y = pos.y * 1.18, cex = 1.75)
    legend(x = mid, y = pos.y *.92, inset=.02, xjust = 0.5, fill=colos, horiz=TRUE, cex=1.25, bty = "n",
           c("insufficient","sufficient", "good"))

  } else {

    plot(density(samp, adjust = 1.75), type = "l", axes = F, xlab = NA, ylab = NA,
         xlim = c(0, 1), ylim = c(0,  peak * 1.2),
         lwd = 3, main = "")
    plotShadePost(dens.post, xx = c(x0, x1, x2, x3),
                              cols = colos, criteria = criteria, blackwhite = blackwhite)
    plotShadePrior(dens.item, xx = c(xx0, xx1, xx2, xx3),
                             cols = colos, criteria = criteria, blackwhite = blackwhite)

    lines(density(item, adjust = 1.75), lty = 4, lwd = 3)

    axis(side = 1, at = seq(0, 1, by = .2), labels = seq(0, 1, by = .2), cex.axis = 1.25, lwd = 1.5)

    arrows(x0 = hdi[1], y0 = peak, x1 = hdi[2], y1 = peak, angle = 90, length = 0.05,
           code = 3, lwd = 1.5)
    text(substitute(paste("95% HDI"[original]: "[", cc1, ", ", cc2,"]"), list(cc1 = round(hdi[1], 3), cc2 = round(hdi[2], 3))),
         x = pos.x, y = peak*1.04, adj = 0)
    text(labels = substitute(paste("median" [original], " = ", cc), list(cc = round(med, 3))), x = pos.x, y = peak*1.09, adj = 0)

    arrows(x0 = hdi2[1], y0 = peak2, x1 = hdi2[2], y1 = peak2, angle = 90, length = 0.05,
           code = 3, lwd = 1.5, lty = 4)
    text(substitute(paste("95% HDI"[ifitem]: "[", cc1, ", ", cc2,"]"), list(cc1 = round(hdi2[1], 3), cc2 = round(hdi2[2], 3))),
         x = pos.x, y = peak*0.94, adj = 0)
    text(labels = substitute(paste("median" [ifitem], " = ", cc), list(cc = round(med2, 3))), x = pos.x, y = peak*.99, adj = 0)
    title(main, line = -1, cex.main = 1.75)
    if (criteria){
      text("insufficient", x = cuts[1]/2, y = peak*-.025, adj = 0.5, cex = 1.25)
      text("sufficient", x = (cuts[1] + cuts[2])/2, y = peak*-.025, adj = .5, cex = 1.25)
      text("good", x = (cuts[2] + 1)/2, y = peak*-.025, adj = .5, cex = 1.25)
      f.post <- plotrix::floating.pie(xpos = pos.x, ypos = pos.y, x = pie.post, radius = rad,
                                      col = colos, startpos = 0)
      plotrix::pie.labels(pos.x, pos.y, radius = rad + 0.006, labels = pie.post.labels, f.post)
      f.item <- plotrix::floating.pie(xpos = pos.x +.2, ypos = pos.y, x = pie.item, radius = rad,
                                      col = colos, startpos = 0)
      plotrix::pie.labels(pos.x +.2, pos.y, radius = rad + .006, labels = pie.item.labels, f.item)

      text("original", x = pos.x, y = pos.y * 1.29, cex = 1.5)
      text("if-item-dropped", x = pos.x + 0.2, y = pos.y * 1.29, cex = 1.5)
      segments(x0 = pos.x -.04, y0 = pos.y * 1.245, x1 = pos.x +.04, y1 = pos.y * 1.245, lwd = 1.75, lty = 1)
      segments(x0 = pos.x + .13, y0 = pos.y * 1.245, x1 = pos.x + .27, y1 = pos.y * 1.245, lwd = 1.75, lty = 4)
      #text("Distribution", x = mid, y = pos.y * 1.41, cex = 1.75)
      legend(x = mid, y = pos.y * 0.77, inset=.02, xjust = .5, fill=colos, horiz=TRUE, cex=1.25, # bty = "n",
             c("insufficient","sufficient", "good"))
    } else {
      legend(x = mid, y = pos.y*0.7, c("original", "if-item-dropped"), lty = c(1, 4), lwd = 2, title = "Distribution", cex = 1.25, bty = "n")
    }
  }
  options(warn = 0)
}

#' plots posterior distributions of chosen estimate and the item-dropped cases in one plot
#' @export
plotifitem_all <- function(x, estimate, ordering = FALSE){
  n.row <- length(unlist(x$bay$ifitem$est[1]))
  posi <- grep(estimate, x$estimates, ignore.case = T)
  main <- paste("If-Item-Dropped Posterior Plot for", estimate)

  dat <- as.data.frame(as.matrix(unlist(x$bay$samp[posi])))
  colnames(dat) <- "value"
  dat$var <- "original"
  dat$colo <- "1"

  dat.del <- t(as.matrix(as.data.frame(x$bay$ifitem$samp[posi])))

  names <- NULL
  for(i in 1:(n.row)){
    names[i] <- paste0("x", i)
  }

  for (i in 1:n.row){
    tmp <- as.data.frame(dat.del[, i])
    colnames(tmp) <- "value"
    tmp$var <- names[i]
    tmp$colo <- "2"
    dat <- rbind(dat, tmp)
  }
  dat$var <- factor(dat$var, levels = unique(dat$var))
  if (ordering){
    est <- as.data.frame(unlist(x$bay$ifitem$est[posi]))
    est[n.row + 1, ] <- 1
    colnames(est) <- "value"
    est$name <- c(names, "original")
    est <- est[order(est$value, decreasing = T), ]
    dat$var <- factor(dat$var, levels = c(est$name))
  }

  ggplot2::ggplot(dat, ggplot2::aes(x = value, y = var, fill = colo)) +
    ggridges::stat_density_ridges(quantile_lines = TRUE, quantiles = c(0.025, 0.5, 0.975), alpha = .85, show.legend = F) +
    ggplot2::theme_linedraw() +
    ggplot2::theme(strip.background = ggplot2::element_rect(fill = "white"),
                   strip.text = ggplot2::element_text(colour = "black")) +
    ggplot2::xlab("\n Reliability") +
    ggplot2::ylab("Item Dropped") +
    ggplot2::scale_y_discrete(expand = ggplot2::expand_scale(add = c(0.25, 1.5))) +
    ggplot2::ggtitle(main) +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, vjust = 4, size = 20),
                   axis.title = ggplot2::element_text(size = 16),
                   axis.text = ggplot2::element_text(size = 12),
                   plot.margin = ggplot2::unit(c(1,1,1,1), "cm"))

}
