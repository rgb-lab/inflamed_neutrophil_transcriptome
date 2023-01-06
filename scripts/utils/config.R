##########
# CONFIG #
##########
#
# This file contains global configuration options for things like plotting.
#
##########

###########
# IMPORTS #
###########
if (!require(MASS)) install.packages("MASS")
if (!require(viridis)) install.packages("viridis")
if (!require(ggplot2)) install.packages("ggplot2")
if (!require(RColorBrewer)) install.packages("RColorBrewer")
if (!require(ggsci)) install.packages("ggsci")

if (!require(BiocManager)) install.packages("BiocManager")
if (!require(circlize)) BiocManager::install("circlize")

###################
# GRAPHING CONFIG #
###################

# THEME
theme_rgb <- function () { 
  theme_classic(base_family="Helvetica", base_size=6) %+replace% 
    theme(
      axis.text.x = element_text(size = 6, color = "black"),
      axis.text.y = element_text(size = 6, color = "black"),
      #	axis.text.x = element_text(size=5, color = "black"),
      text = element_text(size = 6, color = "black"),
      #	axis.title = element_text(size = 6),
      axis.title = element_blank(),
      axis.ticks=element_line(color = "black", size=0.25),
      panel.border = element_rect(colour = "black", fill=NA, size=0.5),
      axis.line = element_blank(),
      strip.text = element_text(size = 6, color = "black", vjust = 1),
      plot.background = element_blank(),
      #	strip.text = element_blank(),
      strip.background = element_blank())
}

# PALETTES
# discrete
# exp_cond_pal <- c(
#   "U" = "#556FB1",
#   "GI" = "#FF725D",
#   "GL" = "#FFC743"
#   )
exp_cond_pal <- c(
  "U" = "#150578",
  "GI" = "#FFA400",
  "GL" = "#D90368"
  )

exp_cond_long_pal <- c(
  "unstim" = "#150578",
  "GI" = "#FFA400",
  "GL" = "#D90368"
)

ice_blue_pal <- c(
  "#03045e",
  "#023e8a",
  "#0077b6",
  "#0096c7",
  "#00b4d8",
  "#48cae4"#,
  # "#90e0ef",
  # "#ade8f4",
  # "#caf0f8"
)

species_pal <- brewer.pal(6, "RdGy")[5:6]
names(species_pal) <- c("Hs", "Mm")

species_pal_long <- brewer.pal(6, "RdGy")[5:6]
names(species_pal_long) <- c("Homo sapiens", "Mus musculus")

species_pal_en <- brewer.pal(6, "RdGy")[5:6]
names(species_pal_en) <- c("Human", "Mouse")

hc_infl_pal <- brewer.pal(2, "Pastel1")[c(3, 1)]
names(hc_infl_pal) <- c("HC", "INFL")

lineage_pal <- pal_lancet()(6)
lineage_pal <- c(lineage_pal[5:6], lineage_pal[1:2], lineage_pal[3:4])
names(lineage_pal) <- c("B cells",
                        "Dendritic cells",
                        "Monocytes",
                        "Neutrophils",
                        "NK cells",
                        "T cells")

lineage_pal_short <- lineage_pal
names(lineage_pal_short) <- ifelse(
  names(lineage_pal) == "B cells",
  "B",
  ifelse(names(lineage_pal) == "Dendritic cells",
         "DC",
         ifelse(names(lineage_pal) == "Monocytes",
                "Mono",
                ifelse (names(lineage_pal) == "Neutrophils",
                        "GN",
                        ifelse (names(lineage_pal) == "NK cells",
                                "NK",
                                "T")
                        )
                )
         )
  )

model_pal <- c("in vitro" = "#fca311", "in vivo" = "#e5e5e5")

# continuous div
# for complex heatmap
generate_rd_bu_colfun <- function (min = -1, max = 1) {
  size = max - min
  ramp <- colorRampPalette(rev(brewer.pal(n = 11, name = "RdBu")))(25)
  colfun <- colorRamp2(seq(min, max, size/24), ramp)
  
  return (colfun)
}
generate_pi_y_g_colfun <- function (min = -1, max = 1) {
  size = max - min
  ramp <- colorRampPalette(rev(brewer.pal(n = 11, name = "PiYG")))(25)
  colfun <- colorRamp2(seq(min, max, size/24), ramp)
  
  return (colfun)
}

generate_blues_w_oranges_colfun <- function (min = -1, max = 1) {
  size = max - min
  ramp <- colorRampPalette(
    c(
      rev(brewer.pal(n = 9, name = "Blues")),
      "#FFFFFF",
      brewer.pal(n = 9, name = "Oranges")
      )
    )(25)
  colfun <- colorRamp2(seq(min, max, size/24), ramp)
  
  return (colfun)
}

# continuous
# for complex heatmap
generate_inferno_colfun <- function(min = 0, max = 1) {
  size <- max-min
  col_fun <- colorRamp2(seq(min, max, size/24), inferno(25))
  
  return(col_fun)
}

generate_blues_colfun <- function (min = 0, max = 1) {
  size = max - min
  ramp <- colorRampPalette(brewer.pal(n = 9, name = "Blues"))(25)
  colfun <- colorRamp2(seq(min, max, size/24), ramp)
  
  return (colfun)
}
generate_greens_colfun <- function (min = 0, max = 1) {
  size = max - min
  ramp <- colorRampPalette(brewer.pal(n = 9, name = "Greens"))(25)
  colfun <- colorRamp2(seq(min, max, size/24), ramp)
  
  return (colfun)
}
generate_or_rd_colfun <- function (min = 0, max = 1) {
  size = max - min
  ramp <- colorRampPalette(rev(brewer.pal(9, "OrRd")[4:9]))(25)
  colfun <- colorRamp2(seq(min, max, size/24), ramp)
  
  return (colfun)
}

# continuous
# for complex heatmap
# TODO: implement general colfun generator
# generate_brewer_colfun <- function (min = -1, max = 1, pal = "RdBu") {
#   
# }

# Custom enhancedVolcano (based on CRAN v1.14.0) function to introduce new vars "labAlpha" and "labPad"
enhancedVolcanoFAR <- function (toptable, lab, x, y, selectLab = NULL, xlim = c(min(toptable[[x]], 
                                                                                    na.rm = TRUE) - 1.5, max(toptable[[x]], na.rm = TRUE) + 
                                                                                  1.5), ylim = c(0, max(-log10(toptable[[y]]), na.rm = TRUE) + 
                                                                                                   5), xlab = bquote(~Log[2] ~ "fold change"), ylab = bquote(~-Log[10] ~ 
                                                                                                                                                               italic(P)), axisLabSize = 18, title = "Volcano plot", subtitle = bquote(italic(EnhancedVolcano)), 
                                caption = paste0("total = ", nrow(toptable), " variables"), 
                                titleLabSize = 18, subtitleLabSize = 14, captionLabSize = 14, 
                                pCutoff = 1e-05, pCutoffCol = y, FCcutoff = 1, FCcutoff_label = 0, cutoffLineType = "longdash", 
                                cutoffLineCol = "black", cutoffLineWidth = 0.4, pointSize = 2, 
                                labSize = 5, labCol = "black", labAlpha = 1, labPad = 0.25, labFace = "plain", boxedLabels = FALSE, 
                                parseLabels = FALSE, shape = 19, shapeCustom = NULL, col = c("grey30", 
                                                                                             "forestgreen", "royalblue", "red2"), colCustom = NULL, 
                                colAlpha = 1/2, colGradient = NULL, colGradientBreaks = c(pCutoff, 
                                                                                          1), colGradientLabels = c("0", "1.0"), colGradientLimits = c(0, 
                                                                                                                                                       1), legendLabels = c("NS", expression(Log[2] ~ FC), 
                                                                                                                                                                            "p-value", expression(p - value ~ and ~ log[2] ~ FC)), 
                                legendPosition = "top", legendLabSize = 14, legendIconSize = 5, 
                                legendDropLevels = TRUE, encircle = NULL, encircleCol = "black", 
                                encircleFill = "pink", encircleAlpha = 3/4, encircleSize = 2.5, 
                                shade = NULL, shadeFill = "grey", shadeAlpha = 1/2, shadeSize = 0.01, 
                                shadeBins = 2, drawConnectors = FALSE, widthConnectors = 0.5, 
                                typeConnectors = "closed", endsConnectors = "first", lengthConnectors = unit(0.01, 
                                                                                                             "npc"), colConnectors = "grey10", max.overlaps = 15, 
                                maxoverlapsConnectors = NULL, min.segment.length = 0, directionConnectors = "both", 
                                arrowheads = TRUE, hline = NULL, hlineType = "longdash", 
                                hlineCol = "black", hlineWidth = 0.4, vline = NULL, vlineType = "longdash", 
                                vlineCol = "black", vlineWidth = 0.4, gridlines.major = TRUE, 
                                gridlines.minor = TRUE, border = "partial", borderWidth = 0.8, 
                                borderColour = "black", raster = FALSE) 
{
  if (!is.numeric(toptable[[x]])) {
    stop(paste(x, " is not numeric!", sep = ""))
  }
  if (!is.numeric(toptable[[pCutoffCol]])) {
    stop(paste(y, " is not numeric!", sep = ""))
  }
  if (raster) {
    has_ggrastr <- !is(try(find.package("ggrastr"), silent = TRUE), 
                       "try-error")
    if (has_ggrastr) {
      geom_point <- ggrastr::geom_point_rast
    }
    else {
      warning("raster disabled, required package \"ggrastr\" not installed")
    }
  }
  if (!is.null(maxoverlapsConnectors)) {
    max.overlaps <- maxoverlapsConnectors
  }
  i <- xvals <- yvals <- Sig <- NULL
  toptable <- as.data.frame(toptable)
  toptable$Sig <- "NS"
  toptable$Sig[(abs(toptable[[x]]) > FCcutoff)] <- "FC"
  toptable$Sig[(toptable[[pCutoffCol]] < pCutoff)] <- "P"
  toptable$Sig[(toptable[[pCutoffCol]] < pCutoff) & (abs(toptable[[x]]) > 
                                                       FCcutoff)] <- "FC_P"
  toptable$Sig <- factor(toptable$Sig, levels = c("NS", "FC", 
                                                  "P", "FC_P"))
  if (pCutoffCol != y) {
    pCutoff = max(toptable[which(toptable[pCutoffCol] <= 
                                   pCutoff), y])
  }
  if (min(toptable[[y]], na.rm = TRUE) == 0) {
    warning(paste("One or more p-values is 0.", "Converting to 10^-1 * current", 
                  "lowest non-zero p-value..."), call. = FALSE)
    toptable[which(toptable[[y]] == 0), y] <- min(toptable[which(toptable[[y]] != 
                                                                   0), y], na.rm = TRUE) * 10^-1
  }
  toptable$lab <- lab
  toptable$xvals <- toptable[[x]]
  toptable$yvals <- toptable[[y]]
  if (!is.null(selectLab)) {
    names.new <- rep(NA, length(toptable$lab))
    indices <- which(toptable$lab %in% selectLab)
    names.new[indices] <- toptable$lab[indices]
    toptable$lab <- names.new
  }
  th <- theme_bw(base_size = 24) + theme(legend.background = element_rect(), 
                                         plot.title = element_text(angle = 0, size = titleLabSize, 
                                                                   face = "bold", vjust = 1), plot.subtitle = element_text(angle = 0, 
                                                                                                                           size = subtitleLabSize, face = "plain", vjust = 1), 
                                         plot.caption = element_text(angle = 0, size = captionLabSize, 
                                                                     face = "plain", vjust = 1), axis.text.x = element_text(angle = 0, 
                                                                                                                            size = axisLabSize, vjust = 1), axis.text.y = element_text(angle = 0, 
                                                                                                                                                                                       size = axisLabSize, vjust = 0.5), axis.title = element_text(size = axisLabSize), 
                                         legend.position = legendPosition, legend.key = element_blank(), 
                                         legend.key.size = unit(0.5, "cm"), legend.text = element_text(size = legendLabSize), 
                                         title = element_text(size = legendLabSize), legend.title = element_blank())
  if (!is.null(colCustom) & !is.null(shapeCustom)) {
    plot <- ggplot(toptable, aes(x = xvals, y = -log10(yvals))) + 
      th + guides(colour = guide_legend(order = 1, override.aes = list(size = legendIconSize)), 
                  shape = guide_legend(order = 2, override.aes = list(size = legendIconSize))) + 
      geom_point(aes(color = factor(names(colCustom)), 
                     shape = factor(names(shapeCustom))), alpha = colAlpha, 
                 size = pointSize, na.rm = TRUE) + scale_color_manual(values = colCustom) + 
      scale_shape_manual(values = shapeCustom)
  }
  else if (!is.null(colCustom) & is.null(shapeCustom) & length(shape) == 
           1) {
    plot <- ggplot(toptable, aes(x = xvals, y = -log10(yvals))) + 
      th + guides(colour = guide_legend(order = 1, override.aes = list(size = legendIconSize)), 
                  shape = guide_legend(order = 2, override.aes = list(size = legendIconSize))) + 
      geom_point(aes(color = factor(names(colCustom))), 
                 alpha = colAlpha, shape = shape, size = pointSize, 
                 na.rm = TRUE) + scale_color_manual(values = colCustom) + 
      scale_shape_manual(guide = TRUE)
  }
  else if (!is.null(colCustom) & is.null(shapeCustom) & length(shape) == 
           4) {
    plot <- ggplot(toptable, aes(x = xvals, y = -log10(yvals))) + 
      th + guides(colour = guide_legend(order = 1, override.aes = list(size = legendIconSize)), 
                  shape = guide_legend(order = 2, override.aes = list(size = legendIconSize))) + 
      geom_point(aes(color = factor(names(colCustom)), 
                     shape = Sig), alpha = colAlpha, size = pointSize, 
                 na.rm = TRUE) + scale_color_manual(values = colCustom) + 
      scale_shape_manual(values = c(NS = shape[1], FC = shape[2], 
                                    P = shape[3], FC_P = shape[4]), labels = c(NS = legendLabels[1], 
                                                                               FC = legendLabels[2], P = legendLabels[3], FC_P = legendLabels[4]), 
                         guide = TRUE, drop = legendDropLevels)
  }
  else if (is.null(colCustom) & !is.null(shapeCustom)) {
    if (is.null(colGradient)) {
      plot <- ggplot(toptable, aes(x = xvals, y = -log10(yvals))) + 
        th + guides(colour = guide_legend(order = 1, 
                                          override.aes = list(size = legendIconSize)), 
                    shape = guide_legend(order = 2, override.aes = list(size = legendIconSize))) + 
        geom_point(aes(color = Sig, shape = factor(names(shapeCustom))), 
                   alpha = colAlpha, size = pointSize, na.rm = TRUE) + 
        scale_color_manual(values = c(NS = col[1], FC = col[2], 
                                      P = col[3], FC_P = col[4]), labels = c(NS = legendLabels[1], 
                                                                             FC = legendLabels[2], P = legendLabels[3], 
                                                                             FC_P = legendLabels[4]), drop = legendDropLevels) + 
        scale_shape_manual(values = shapeCustom)
    }
    else {
      plot <- ggplot(toptable, aes(x = xvals, y = -log10(yvals))) + 
        th + guides(shape = guide_legend(order = 2, 
                                         override.aes = list(size = legendIconSize))) + 
        geom_point(aes(color = Sig, shape = factor(names(shapeCustom))), 
                   alpha = colAlpha, size = pointSize, na.rm = TRUE) + 
        scale_colour_gradient(low = colGradient[1], 
                              high = colGradient[2], limits = colGradientLimits, 
                              breaks = colGradientBreaks, labels = colGradientLabels)
      scale_shape_manual(values = shapeCustom)
    }
  }
  else if (is.null(colCustom) & is.null(shapeCustom) & length(shape) == 
           1) {
    if (is.null(colGradient)) {
      plot <- ggplot(toptable, aes(x = xvals, y = -log10(yvals))) + 
        th + guides(colour = guide_legend(order = 1, 
                                          override.aes = list(shape = shape, size = legendIconSize))) + 
        geom_point(aes(color = Sig), alpha = colAlpha, 
                   shape = shape, size = pointSize, na.rm = TRUE) + 
        scale_color_manual(values = c(NS = col[1], FC = col[2], 
                                      P = col[3], FC_P = col[4]), labels = c(NS = legendLabels[1], 
                                                                             FC = legendLabels[2], P = legendLabels[3], 
                                                                             FC_P = legendLabels[4]), drop = legendDropLevels)
    }
    else {
      plot <- ggplot(toptable, aes(x = xvals, y = -log10(yvals))) + 
        th + geom_point(aes(color = yvals), alpha = colAlpha, 
                        shape = shape, size = pointSize, na.rm = TRUE) + 
        scale_colour_gradient(low = colGradient[1], 
                              high = colGradient[2], limits = colGradientLimits, 
                              breaks = colGradientBreaks, labels = colGradientLabels)
    }
  }
  else if (is.null(colCustom) & is.null(shapeCustom) & length(shape) == 
           4) {
    if (is.null(colGradient)) {
      plot <- ggplot(toptable, aes(x = xvals, y = -log10(yvals))) + 
        th + guides(colour = guide_legend(order = 1, 
                                          override.aes = list(shape = c(NS = shape[1], 
                                                                        FC = shape[2], P = shape[3], FC_P = shape[4]), 
                                                              size = legendIconSize))) + geom_point(aes(color = Sig, 
                                                                                                        shape = Sig), alpha = colAlpha, size = pointSize, 
                                                                                                    na.rm = TRUE) + scale_color_manual(values = c(NS = col[1], 
                                                                                                                                                  FC = col[2], P = col[3], FC_P = col[4]), labels = c(NS = legendLabels[1], 
                                                                                                                                                                                                      FC = legendLabels[2], P = legendLabels[3], FC_P = legendLabels[4]), 
                                                                                                                                       drop = legendDropLevels) + scale_shape_manual(values = c(NS = shape[1], 
                                                                                                                                                                                                FC = shape[2], P = shape[3], FC_P = shape[4]), 
                                                                                                                                                                                     guide = FALSE, drop = legendDropLevels)
    }
    else {
      plot <- ggplot(toptable, aes(x = xvals, y = -log10(yvals))) + 
        th + geom_point(aes(color = yvals, shape = Sig), 
                        alpha = colAlpha, size = pointSize, na.rm = TRUE) + 
        scale_colour_gradient(low = colGradient[1], 
                              high = colGradient[2], limits = colGradientLimits, 
                              breaks = colGradientBreaks, labels = colGradientLabels) + 
        scale_shape_manual(values = c(NS = shape[1], 
                                      FC = shape[2], P = shape[3], FC_P = shape[4]), 
                           guide = FALSE, drop = legendDropLevels)
    }
  }
  plot <- plot + xlab(xlab) + ylab(ylab) + xlim(xlim[1], xlim[2]) + 
    ylim(ylim[1], ylim[2]) + geom_vline(xintercept = c(-FCcutoff, 
                                                       FCcutoff), linetype = cutoffLineType, colour = cutoffLineCol, 
                                        size = cutoffLineWidth) + geom_hline(yintercept = -log10(pCutoff), 
                                                                             linetype = cutoffLineType, colour = cutoffLineCol, size = cutoffLineWidth)
  plot <- plot + labs(title = title, subtitle = subtitle, 
                      caption = caption)
  if (!is.null(vline)) {
    plot <- plot + geom_vline(xintercept = vline, linetype = vlineType, 
                              colour = vlineCol, size = vlineWidth)
  }
  if (!is.null(hline)) {
    plot <- plot + geom_hline(yintercept = -log10(hline), 
                              linetype = hlineType, colour = hlineCol, size = hlineWidth)
  }
  if (border == "full") {
    plot <- plot + theme(panel.border = element_rect(colour = borderColour, 
                                                     fill = NA, size = borderWidth))
  }
  else if (border == "partial") {
    plot <- plot + theme(axis.line = element_line(size = borderWidth, 
                                                  colour = borderColour), panel.border = element_blank(), 
                         panel.background = element_blank())
  }
  else {
    stop("Unrecognised value passed to 'border'. Must be 'full' or 'partial'")
  }
  if (gridlines.major) {
    plot <- plot + theme(panel.grid.major = element_line())
  }
  else {
    plot <- plot + theme(panel.grid.major = element_blank())
  }
  if (gridlines.minor) {
    plot <- plot + theme(panel.grid.minor = element_line())
  }
  else {
    plot <- plot + theme(panel.grid.minor = element_blank())
  }
  if (!boxedLabels) {
    if (drawConnectors && is.null(selectLab)) {
      if (arrowheads) {
        arr <- arrow(length = lengthConnectors, type = typeConnectors, 
                     ends = endsConnectors)
      }
      else {
        arr <- NULL
      }
      plot <- plot + geom_text_repel(data = subset(toptable, 
                                                   toptable[[y]] < pCutoff & abs(toptable[[x]]) > 
                                                     FCcutoff_label), aes(label = subset(toptable, toptable[[y]] < 
                                                                                     pCutoff & abs(toptable[[x]]) > FCcutoff_label)[["lab"]]), 
                                     xlim = c(NA, NA), ylim = c(NA, NA), size = labSize, 
                                     segment.color = colConnectors, segment.size = widthConnectors, 
                                     arrow = arr, colour = labCol, alpha = labAlpha, label.padding = unit(labPad, "lines"), fontface = labFace, 
                                     parse = parseLabels, na.rm = TRUE, direction = directionConnectors, 
                                     max.overlaps = max.overlaps, min.segment.length = min.segment.length)
    }
    else if (drawConnectors && !is.null(selectLab)) {
      if (arrowheads) {
        arr <- arrow(length = lengthConnectors, type = typeConnectors, 
                     ends = endsConnectors)
      }
      else {
        arr <- NULL
      }
      plot <- plot + geom_text_repel(data = subset(toptable, 
                                                   !is.na(toptable[["lab"]])), aes(label = subset(toptable, 
                                                                                                  !is.na(toptable[["lab"]]))[["lab"]]), xlim = c(NA, 
                                                                                                                                                 NA), ylim = c(NA, NA), size = labSize, segment.color = colConnectors, 
                                     segment.size = widthConnectors, arrow = arr, 
                                     colour = labCol,  alpha = labAlpha, label.padding = unit(labPad, "lines"), fontface = labFace, parse = parseLabels, 
                                     na.rm = TRUE, direction = directionConnectors, 
                                     max.overlaps = max.overlaps, min.segment.length = min.segment.length)
    }
    else if (!drawConnectors && !is.null(selectLab)) {
      plot <- plot + geom_text(data = subset(toptable, 
                                             !is.na(toptable[["lab"]])), aes(label = subset(toptable, 
                                                                                            !is.na(toptable[["lab"]]))[["lab"]]), size = labSize, 
                               check_overlap = TRUE, colour = labCol, alpha = labAlpha, label.padding = unit(labPad, "lines"), fontface = labFace, 
                               parse = parseLabels, na.rm = TRUE)
    }
    else if (!drawConnectors && is.null(selectLab)) {
      plot <- plot + geom_text(data = subset(toptable, 
                                             toptable[[y]] < pCutoff & abs(toptable[[x]]) > 
                                               FCcutoff_label), aes(label = subset(toptable, toptable[[y]] < 
                                                                               pCutoff & abs(toptable[[x]]) > FCcutoff_label)[["lab"]]), 
                               size = labSize, check_overlap = TRUE, colour = labCol, alpha = labAlpha, label.padding = unit(labPad, "lines"), 
                               fontface = labFace, parse = parseLabels, na.rm = TRUE)
    }
  }
  else {
    if (drawConnectors && is.null(selectLab)) {
      if (arrowheads) {
        arr <- arrow(length = lengthConnectors, type = typeConnectors, 
                     ends = endsConnectors)
      }
      else {
        arr <- NULL
      }
      plot <- plot + geom_label_repel(data = subset(toptable, 
                                                    toptable[[y]] < pCutoff & abs(toptable[[x]]) > 
                                                      FCcutoff_label), aes(label = subset(toptable, toptable[[y]] < 
                                                                                      pCutoff & abs(toptable[[x]]) > FCcutoff_label)[["lab"]]), 
                                      xlim = c(NA, NA), ylim = c(NA, NA), size = labSize, 
                                      segment.color = colConnectors, segment.size = widthConnectors, 
                                      arrow = arr, colour = labCol, alpha = labAlpha, label.padding = unit(labPad, "lines"), fontface = labFace, 
                                      parse = parseLabels, na.rm = TRUE, direction = directionConnectors, 
                                      max.overlaps = max.overlaps, min.segment.length = min.segment.length)
    }
    else if (drawConnectors && !is.null(selectLab)) {
      if (arrowheads) {
        arr <- arrow(length = lengthConnectors, type = typeConnectors, 
                     ends = endsConnectors)
      }
      else {
        arr <- NULL
      }
      plot <- plot + geom_label_repel(data = subset(toptable, 
                                                    !is.na(toptable[["lab"]])), aes(label = subset(toptable, 
                                                                                                   !is.na(toptable[["lab"]]))[["lab"]]), xlim = c(NA, 
                                                                                                                                                  NA), ylim = c(NA, NA), size = labSize, segment.color = colConnectors, 
                                      segment.size = widthConnectors, arrow = arr, 
                                      colour = labCol, alpha = labAlpha, label.padding = unit(labPad, "lines"), fontface = labFace, parse = parseLabels, 
                                      na.rm = TRUE, direction = directionConnectors, 
                                      max.overlaps = max.overlaps, min.segment.length = min.segment.length)
    }
    else if (!drawConnectors && !is.null(selectLab)) {
      plot <- plot + geom_label(data = subset(toptable, 
                                              !is.na(toptable[["lab"]])), aes(label = subset(toptable, 
                                                                                             !is.na(toptable[["lab"]]))[["lab"]]), size = labSize, 
                                colour = labCol, alpha = labAlpha, label.padding = unit(labPad, "lines"), fontface = labFace, parse = parseLabels, 
                                na.rm = TRUE)
    }
    else if (!drawConnectors && is.null(selectLab)) {
      plot <- plot + geom_label(data = subset(toptable, 
                                              toptable[[y]] < pCutoff & abs(toptable[[x]]) > 
                                                FCcutoff_label), aes(label = subset(toptable, toptable[[y]] < 
                                                                                pCutoff & abs(toptable[[x]]) > FCcutoff_label)[["lab"]]), 
                                size = labSize, colour = labCol, alpha = labAlpha, label.padding = unit(labPad, "lines"), fontface = labFace, 
                                parse = parseLabels, na.rm = TRUE)
    }
  }
  if (!is.null(encircle)) {
    if (is(try(find.package("ggalt"), silent = TRUE), "try-error")) {
      stop("Please install package \"ggalt\" to access the \"encircle\" features")
    }
    plot <- plot + ggalt::geom_encircle(data = subset(toptable, 
                                                      rownames(toptable) %in% encircle), colour = encircleCol, 
                                        fill = encircleFill, alpha = encircleAlpha, size = encircleSize, 
                                        show.legend = FALSE, na.rm = TRUE)
  }
  if (!is.null(shade)) {
    plot <- plot + stat_density2d(data = subset(toptable, 
                                                rownames(toptable) %in% shade), fill = shadeFill, 
                                  alpha = shadeAlpha, geom = "polygon", contour = TRUE, 
                                  size = shadeSize, bins = shadeBins, show.legend = FALSE, 
                                  na.rm = TRUE)
  }
  plot <- plot + coord_cartesian(clip = "off")
  return(plot)
}

# UTILS

theme_set(theme_bw(base_size = 16))

# thanks Kamil Slowikowski
# Get density of points in 2 dimensions.
# @param x A numeric vector.
# @param y A numeric vector.
# @param n Create a square n by n grid to compute density.
# @return The density within each square.
get_density <- function(x, y, ...) {
  dens <- MASS::kde2d(x, y, ...)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}

add_label <- function(figure, label, fontsize = 12) {
  return(arrangeGrob(figure,
                     top = textGrob(label,
                                    x = unit(0, "npc"),
                                    y = unit(1, "npc"),
                                    just=c("left","top"),
                                    gp=gpar(col="black",
                                            fontfamily="Helvetica",
                                            fontface = "bold",
                                            fontsize = fontsize))))
}


