mywormplot <- function (m = NULL, residuals = NULL, 
                        age = NA, n.inter = 1, 
                        y.limits = c(-1, 1)) 
{
  if (inherits(m, "gamlss")) 
    residuals <- residuals(m)
  mm <- tibble::tibble(x = seq(-4, 4, l = 1000), yu = 1.96 * 
                         sqrt(stats::pnorm(.data$x) * 
                                (1 - stats::pnorm(.data$x))/length(residuals))/stats::dnorm(.data$x),
                       yl = -.data$yu)
  mm <- reshape2::melt(mm, id.var = "x")
  tmp <- data.frame(residuals = residuals, age = age)
  if (n.inter > 1) {
    if (all(is.na(age)))
      stop("intervals only possible of a vector of ages  is given")
    tmp$ag <- cut(tmp$age, n.inter)
  }
  else {
    tmp$ag <- "all ages"
  }
  tmp <- dplyr::group_by(tmp, .data$ag) %>% tidyr::nest()
  tmp <- dplyr::mutate(tmp, qq = purrr::map(.data$data, function(x) {
    qq <- as.data.frame(stats::qqnorm(x$residuals, plot.it = F))
    qq$y <- qq$y - qq$x
    qq
  }))
  tmp <- tidyr::unnest(tmp, .data$qq)
  ggplot2::ggplot(tmp, ggplot2::aes_string(x = "x", y = "y")) + 
    ggplot2::geom_point(shape = 21, size = 3, colour = "#2171B5") + 
    ggplot2::geom_line(data = mm, inherit.aes = F, 
                       ggplot2::aes_string(x = "x", 
                                           y = "value", group = "variable"), 
                       linetype = 2, size=1.2, colour = "forestgreen") + 
    ggplot2::geom_vline(xintercept = 0, colour = "firebrick3", 
                        linetype = 4, size=1.2) + ggplot2::scale_y_continuous(limits = 
                                                                                y.limits) + 
    # ggplot2::facet_wrap(~ag) + 
    ggplot2::labs(x = "Unit normal quantiles", 
                  y = "Deviation") 
  #+ ggplot2::theme_minimal() + ggplot2::theme(legend.position = "none", )
  #                                                                                                         axis.title = ggplot2::element_text(size = 13, colour = "black"), 
  #                                                                                                         axis.text = ggplot2::element_text(size = 13, colour = "black"), 
  #                                                                                                         title = ggplot2::element_text(colour = "black"))
}