# Figures of SSDFGS(). They need the optional package ggplot2. Colors are from the
# Okabe-Ito palette (vermillion #D55E00 for the marked points, blue #0072B2 for the fit).

.plot_operating_curve = function(OC.fit, anno){
  nt = r.score = NULL  # column names used inside aes()
  ggplot2::ggplot() +
    ggplot2::geom_line(ggplot2::aes(x = nt, y = r.score), OC.fit) +
    ggplot2::geom_point(ggplot2::aes(x = nt, y = r.score), anno, color = "#D55E00") +
    ggplot2::annotate("segment", x = -Inf, xend = anno$nt, y = anno$r.score, yend = anno$r.score, color = "#D55E00", linetype = "dashed") +
    ggplot2::annotate("segment", x = anno$nt, xend = anno$nt, y = anno$r.score, yend = -Inf, color = "#D55E00", linetype = "dashed") +
    ggplot2::geom_text(ggplot2::aes(x = nt, y = r.score, label = paste0(" (", nt, ", ", round(r.score, digits = 2), ")")), anno, color = "#D55E00", hjust = "left", vjust = "top") +
    ggplot2::labs(x = expression('Training set size (n'[t]*')'), y = expression('RErs(n'[t]*')'), title = "Operating curve")
}

.plot_growth_curve = function(par, GC.fit){
  nt = r.score = n = r = NULL  # column names used inside aes()
  gc.title = bquote("Logistic growth curve: " ~ y == frac(.(round(par$alpha, 3)), 1 + exp(.(round(par$beta, 3)) - .(round(par$gamma, 3)) * x)))
  ggplot2::ggplot() +
    ggplot2::geom_point(ggplot2::aes(x = n, y = r), data = par$sim, size = 0.8, alpha = 0.5) +
    ggplot2::geom_line(ggplot2::aes(x = nt, y = r.score), data = GC.fit, color = "#0072B2") +
    ggplot2::labs(title = gc.title, x = expression('Training set size (n'[t]*')'), y = "r-score")
}
