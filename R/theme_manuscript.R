# ======================================================================
# theme_manuscript.R
#
# Purpose  : House style for the figures of the paper, for print. No grid and no panel border; axis lines capped at the outer breaks; words in Latin Modern Roman,
#            which matches the Computer Modern text of the manuscript; axis titles and strip labels typeset by LaTeX
#            (xdvir), so math in a figure looks as it does in the text; white background; curves labeled directly where
#            they are apart, and otherwise a legend inside the plot area, never beside it.
# Colors   : series in the order purple, red, then Okabe-Ito blue, bluish green, orange, reddish purple, which keeps
#            neighbouring series separable for color-blind readers and in grayscale; ink for data that are not a
#            series (points, reference curves); muted gray for reference lines.
# Provides : theme_manuscript(), ms_col, ms_series, ms_serif, ms_minus(), ms_pow10(), legend_inside(), ms_save()
# Requires : ggplot2 >= 4.0, systemfonts, ragg; xdvir and a TeX installation for typeset labels (falls back to plain
#            text without them)
# Use      : source("R/theme_manuscript.R"); ggplot(...) + theme_manuscript(); ms_save("figures/x.pdf", p, 7.2, 3.3)
# ======================================================================
suppressPackageStartupMessages(library(ggplot2))

ms_col <- c(ink = "#1F1F1F", note = "#4A4A4A", muted = "#8A8A8A", rule = "#BDBDBD", paper = "#FFFFFF",
            purple = "#6A2C91", red = "#C2453E", blue = "#0072B2", green = "#009E73", orange = "#E69F00",
            pink = "#CC79A7")
ms_series <- unname(ms_col[c("purple", "red", "blue", "green", "orange", "pink")])

# Latin Modern ships with every TeX installation; register it for ragg so that words match the manuscript text
.ms_font_file <- function(f) {
  p <- tryCatch(suppressWarnings(system2("kpsewhich", f, stdout = TRUE, stderr = FALSE)), error = function(e) "")
  if (length(p) && nzchar(p[1]) && file.exists(p[1])) p[1] else NA_character_
}
.ms_register <- function(name, files) {
  if (name %in% systemfonts::registry_fonts()$family) return(name)
  f <- vapply(files, .ms_font_file, "")
  if (anyNA(f)) return("serif")
  systemfonts::register_font(name, f[1], f[2], f[3], f[4])
  name
}
ms_serif <- .ms_register("Latin Modern Roman",
                         c("lmroman10-regular.otf", "lmroman10-bold.otf", "lmroman10-italic.otf", "lmroman10-bolditalic.otf"))
ms_latex <- requireNamespace("xdvir", quietly = TRUE) && nzchar(Sys.which("latex"))
ms_text_size <- 0.85 * 9 / .pt            # in-panel labels at 85% of the 9 pt base size, in ggplot2 size units (mm)
if (ms_latex) {
  suppressPackageStartupMessages(library(xdvir))
  update_geom_defaults(xdvir::GeomLatex, list(family = ms_serif, colour = ms_col[["ink"]], size = ms_text_size))
  # xdvir keeps typeset labels in memory only, so every new R session reran TeX for every label (seconds each);
  # keep them on disk too, under xdvir's own key (text, engine, packages)
  local({
    ns <- asNamespace("xdvir")
    if (!isTRUE(attr(get("getDVI", ns), "ms_cache"))) {
      orig <- get("getDVI", ns); key_of <- get("typesetKey", ns)
      md5 <- function(x) {                       # base-R md5 of a character vector (no extra package)
        f <- tempfile(); on.exit(unlink(f))
        writeLines(enc2utf8(as.character(x)), f, useBytes = TRUE)
        unname(tools::md5sum(f))
      }
      dir <- file.path(tools::R_user_dir("h-index-manuscript", "cache"), "xdvir")
      dir.create(dir, recursive = TRUE, showWarnings = FALSE)
      cached <- function(tex, engine, packages, texFile, sig = TRUE) {
        f <- file.path(dir, paste0(md5(key_of(tex, engine, packages)), ".rds"))
        if (file.exists(f)) return(readRDS(f))
        dvi <- orig(tex, engine, packages, texFile, sig)
        saveRDS(dvi, f)
        dvi
      }
      attr(cached, "ms_cache") <- TRUE
      unlockBinding("getDVI", ns); assign("getDVI", cached, envir = ns); lockBinding("getDVI", ns)
    }
  })
}

theme_manuscript <- function(base_size = 9, family = ms_serif, latex = ms_latex) {
  k <- ms_col; ink <- k[["ink"]]
  tex <- if (latex) function(...) xdvir::element_latex(..., width = unit(30, "in")) else element_text
  t <- theme_minimal(base_size = base_size, base_family = family, ink = ink, paper = k[["paper"]], accent = k[["red"]]) +
    theme(geom = element_geom(linewidth = 0.7, pointsize = 1.2, family = family, fontsize = 0.85 * base_size),
          palette.colour.discrete = ms_series, palette.fill.discrete = ms_series,
          panel.grid = element_blank(),
          axis.line = element_line(colour = ink, linewidth = 0.35),
          axis.ticks = element_line(colour = ink, linewidth = 0.35), axis.ticks.length = unit(3, "pt"),
          axis.title = element_text(colour = ink, size = rel(1)),
          axis.title.x = tex(margin = margin(t = 4)), axis.title.y = tex(angle = 90, margin = margin(r = 4)),
          axis.text = tex(colour = alpha(ink, 0.9), size = rel(0.9)),
          # ggplot2 4 sets the bottom and left tick labels itself, so they take the typeset element explicitly
          axis.text.x.bottom = tex(margin = margin(t = 1.8), vjust = 1),
          axis.text.y.left = tex(margin = margin(r = 1.8), hjust = 1),
          strip.text = tex(colour = ink, size = rel(0.95), hjust = 0, margin = margin(b = 3)),
          panel.spacing = unit(10, "pt"),
          legend.position = "none", legend.title = element_blank(), legend.text = tex(size = rel(0.9)),
          legend.key.width = unit(18, "pt"), legend.key.height = unit(9, "pt"), legend.spacing.y = unit(1, "pt"),
          plot.title = tex(colour = ink, size = rel(1), hjust = 0, margin = margin(b = 4)),
          plot.title.position = "plot",
          plot.margin = margin(4, 8, 4, 4))
  list(t, guides(x = guide_axis(cap = "both"), y = guide_axis(cap = "both")))
}

# tick labels with a true minus sign instead of a hyphen (math-mode minus when TeX typesets the labels, U+2212
# otherwise): scale_y_continuous(labels = ms_minus)
ms_minus <- function(x) {
  s <- format(x, trim = TRUE, drop0trailing = FALSE)
  if (ms_latex) ifelse(grepl("^-", s), paste0("$", s, "$"), s) else sub("^-", "−", s)
}

# power-of-ten tick labels: scale_x_log10(breaks = 10^(2:6), labels = ms_pow10(2:6))
ms_pow10 <- function(k) if (ms_latex) sprintf("$10^{%d}$", k) else parse(text = paste0("10^", k))

# a legend that labels on the curves cannot replace goes inside the plot area
legend_inside <- function(x = 0.98, y = 0.02) {
  theme(legend.position = "inside", legend.position.inside = c(x, y), legend.justification = c(x, y))
}

# vector PDF through cairo for the manuscript (text stays text: TeX-typeset labels embed Latin Modern with a Unicode
# map, so words and numbers in a figure can be selected and searched); a .png file name gives a print-resolution PNG
# through ragg instead
ms_save <- function(file, plot, width, height, dpi = 600) {
  if (grepl("[.]pdf$", file)) {
    ggsave(file, plot, width = width, height = height, device = grDevices::cairo_pdf, bg = ms_col[["paper"]])
  } else {
    ggsave(file, plot, width = width, height = height, dpi = dpi, device = ragg::agg_png, bg = ms_col[["paper"]])
  }
}
