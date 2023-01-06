#!/bin/env Rscript

###################################
###################################
## SCRIPT THAT GENERATES FIGURES ##
###################################
###################################

#
# This script requires all analysis notebooks to have been run previous to
# its execution.
#

if (!require(gridExtra)) install.packages("gridExtra")
if (!require(cowplot)) install.packages("cowplot")
if (!require(magrittr)) install.packages("magrittr")
if (!require(writexl)) install.packages("writexl")
if (!require(gridtext)) remotes::install_github("wilkelab/gridtext") #To make element_markdown italic compatible with R4.2.0 (not yet in CRAN)


source("scripts/utils/figure_builder.R")

################
################
## SUBFIGURES ##
################
################

# ############
# # FIGURE 2 #
# ############
# ############
# 
# if (!dir.exists("figures/figure_2")) dir.create("figures/figure_2")
# 
# #############
# # FIGURE 2A #
# #############
# 
# fig_2a_1 <- build_fig_2a_1()
# fig_2a_2 <- build_fig_2a_2()
# 
# pdf("figures/figure_2/figure_2a.pdf",
#     width = 10,
#     height = 5)
# ggarrange(fig_2a_1,
#           fig_2a_2,
#           common.legend = TRUE,
#           legend = "right") %>%
#   annotate_figure(
#     bottom = text_grob("Number of genes compared",
#                        hjust = 0.5, x = 0.45),
#     left = text_grob("Fraction of reads mapped to top n genes", color = "black", rot = 90)
#   )
# dev.off()
# 
# #############
# # FIGURE 2B #
# #############
# # fig_2b_1 <- build_fig_2b_1()
# 
# pdf("figures/figure_2/figure_2b.pdf",
#     width = 10,
#     height = 5)
# fig_2b_1
# dev.off()
# 
# #############
# # FIGURE 2C #
# #############
# fig_2c_1 <- build_fig_2c_1()
# 
# pdf("figures/figure_2/figure_2c_1.pdf",
#     width = 10,
#     height = 10*9/16)
# fig_2c_1
# dev.off()
# 
# 
# #############
# # FIGURE 2C #
# #############
# fig_2c_2 <- build_fig_2c_2()
# 
# pdf("figures/figure_2/figure_2c_2.pdf",
#     width = 10,
#     height = 10*9/16)
# fig_2c_2
# dev.off()
# 
# #############
# # FIGURE 2D #
# #############
# fig_2d_2 <- build_fig_2d_2()
# 
# pdf("figures/figure_2/figure_2d.pdf",
#     width = 10*0.75,
#     height = 7*0.75)
# fig_2d_2
# dev.off()
# 
# #############
# # FIGURE 2E #
# #############
# fig_2e <- build_fig_2e()
# pdf("figures/figure_2/figure_2e.pdf",
#     width = 15,
#     height = 15*9/16)
# fig_2e
# dev.off()


############
# FIGURE 3 #
############
############

# if (!dir.exists("figures/figure_3")) dir.create("figures/figure_3")

#############
# FIGURE 3A #
# #############
# # all genes
# fig_3a_1 <- build_fig_3a_1()
# # lineage specific genes (from haemopedia)
# fig_3a_2 <- build_fig_3a_2()
# # transcription factors
# fig_3a_3 <- build_fig_3a_3()
# 
# pdf("figures/figure_3/figure_3a.pdf",
#     width = 15,
#     height = 5)
# ggarrange(
#   fig_3a_1,
#   fig_3a_2,
#   fig_3a_3,
#   ncol = 3
# )
# dev.off()

#############
# FIGURE 3B #
#############
# fig_3b <- build_fig_3b()
# pdf("figures/figure_3/figure_3b.pdf",
#     width = 5,
#     height = 5)
# fig_3b
# dev.off()

#############
# FIGURE 3C #
#############
# fig_3c <- build_fig_3c()
# pdf("figures/figure_3/figure_3c.pdf",
#     width = 5,
#     height = 5)
# fig_3c
# dev.off()
# 

#############
# FIGURE 3D #
#############
# fig_3d <- build_fig_3d()
# pdf("figures/figure_3/figure_3d.pdf",
#     width = 5,
#     height = 5)
# fig_3d
# dev.off()

#############
# FIGURE 3E #
#############
# fig_3e <- build_fig_3e()
# pdf("figures/figure_3/figure_3e.pdf",
#     width = 5,
#     height = 5)
# fig_3e
# dev.off()


############
# FIGURE 4 #
############
############

# if (!dir.exists("figures/figure_4")) dir.create("figures/figure_4")

#############
# FIGURE 4C #
#############
# pdf("figures/figure_4/figure_4c.pdf",
#     width = 10,
#     height = 10*9/16)
# 
# 
# draw_fig_4c_2()
# 
# dev.off()

# 
# #############
# # FIGURE 4D #
# #############
# fig_4d_1 <- build_fig_4d_1()
# pdf("figures/figure_4/figure_4d.pdf",
#     width = 10,
#     height = 10*9/20)
# 
# fig_4d_1
# 
# dev.off()
# 
# #############
# # FIGURE 4E #
# #############
# fig_4e_1 <- build_fig_4e_1()
# fig_4e_2 <- build_fig_4e_2()
# 
# pdf("figures/figure_4/figure_4e.pdf",
#     width = 12,
#     height = 12*9/16)
# 
# fig_4e_1
# fig_4e_2
# 
# dev.off()
# 
# #############
# # FIGURE 4F #
# #############
# fig_4f_1 <- build_fig_4f_1()
# fig_4f_2 <- build_fig_4f_2()
# 
# pdf("figures/figure_4/figure_4f.pdf",
#     width = 12,
#     height = 12*9/16)
# 
# ggarrange(
#   fig_4f_1,
#   fig_4f_2,
#   ncol = 2
# )
# 
# dev.off()
# 
# #############
# # FIGURE 4G AND H #
# #############
# fig_4g <- build_fig_4g()
# fig_4h <- build_fig_4h()
# 
# pdf("figures/figure_4/figure_4gh.pdf",
#     width = 12,
#     height = 12*9/16)
# ggarrange(
#   fig_4g,
#   fig_4h,
#   ncol = 2
# )
# 
# dev.off()

############
# FIGURE  #
############
############

# if (!dir.exists("figures/figure_5")) dir.create("figures/figure_5")


# #############
# # FIGURE 5A #
# #############
# fig_5a_1 <- build_fig_5a_1()
# fig_5a_2 <- build_fig_5a_2()
# 
# pdf("figures/figure_5/figure_5a.pdf",
#     width = 14,
#     height = 7)
# 
# ggarrange(
#   fig_5a_1,
#   fig_5a_2,
#   ncol = 2
# )
# 
# dev.off()
# 
# #############
# # FIGURE 5C #
# #############
# fig_5c <- build_fig_5c()
# 
# pdf("figures/delete/regulators_3-7_hm.pdf",
#     width = 13,
#     height = 13*9/20)
# 
# fig_5c
# 
# dev.off()


################
################
# MAIN FIGURES #
################
################

############
# FIGURE 1 #
############

# plot panels

fig_1a1 <- plot_grid(
  ggdraw() + 
    draw_image("figures/figure_1a.png"),
  labels = c("A"),
  label_size = 8)

fig_1a <- build_fig_1a_3()
fig_1a_lgd <- get_legend(build_fig_1a_3(return_color_lgd = TRUE))

fig_1d <- build_fig_1d_1()
# fig_2d <- build_fig_2d_2()
fig_1d_grob <- grid.grabExpr(draw(fig_1d),
                             height = 9.375*1/(1+1+0.75),
                             width = 6.75/2)
fig_1d_legend <- build_fig_1d_1(return_lgd = TRUE)
# fig_1d_legend <- build_fig_1d_2(return_lgd = TRUE)

fig_1e <- build_fig_1e()
fig_1e_grob <- grid.grabExpr(draw(fig_1e))
fig_1e_legend <- build_fig_1e(return_lgd = TRUE)

fig_1f <- build_fig_1f()

# plot legends

legends <- grid.grabExpr({
  #grid.rect()  # border
  draw(fig_1d_legend, x = unit(0.25, "npc"), y = unit(0.8, "npc"), just = c("center", "center"))
  draw(fig_1e_legend, x = unit(0.75, "npc"), y = unit(0.8, "npc"), just = c("center", "center"))
})


fig_1_legends <- plot_grid(
  NULL,
  legends,
  fig_1a_lgd,
  NULL,
  ncol = 1,
  rel_heights = c(0.5,1,0.1,0.5)
) #+ theme(plot.background = element_rect(color = "black", size = 0.5))


fig_1_row_1 <- plot_grid(
  fig_1a1,
  fig_1a, labels = c("A","B"),
  label_size = 8,
  rel_widths = c(2,1)
)

fig_1_row_2 <- plot_grid(
  fig_1d_grob,
  fig_1e_grob,
  nrow = 1,
  labels = c("C", "D"),
  label_size = 8,
  rel_widths = c(1,1)
)

fig_1_row_3 <- plot_grid(
  fig_1f,
  fig_1_legends,
  nrow = 1, ncol = 2,
  rel_widths = c(2, 1),
  rel_heights = c(1,1),
  labels = c("E", ""),
  label_size = 8
)

fig_1_cowplot <- plot_grid(
  fig_1_row_1,
  fig_1_row_2,
  fig_1_row_3,
  ncol = 1,
  rel_heights = c(0.75,1,1)
)

save_plot("figures/Figure 1.pdf",
          fig_1_cowplot,
          ncol = 4,
          nrow = 3,
          base_width = 6.75/4,
          base_height = 9.375/3
          )

# layout_matrix <-
#   rbind(
#     c(1, 2, 5, 5),
#     c(1, 2, 5, 5),
#     c(3, 4, 5, 5),
#     c(3, 4, 5, 5),
#     c(6, 6, 7, 7),
#     c(6, 6, 7, 7),
#     c(6, 6, 7, 7),
#     c(6, 6, 7, 7),
#     c(8, 8, 8, 9),
#     c(8, 8, 8, 9),
#     c(8, 8, 8, 9),
#     c(8, 8, 8, 9)
#   )
# 
# pdf("figures/figure_1.pdf",
#     height = 20*0.9,
#     width = 15*0.9)
# grid.arrange(
#   fig_1a_1 %>% add_label("A"),
#   fig_1a_2 %>% add_label(""),
#   fig_1b_1 %>% add_label("B"),
#   fig_1b_2 %>% add_label(""),
#   fig_1c_1 %>% add_label("C"),
#   fig_1d_grob %>% add_label("D"),
#   fig_1e_grob %>% add_label("E"),
#   fig_1f %>% add_label("F"),
#   grid.text("Add another legend here?"),
#   layout_matrix = layout_matrix
# )
# dev.off()

############
# FIGURE 2 #
############

# all genes
fig_2a_1 <- build_fig_2a_1()
# lineage specific genes (from haemopedia)
fig_2a_2 <- build_fig_2a_2()
# transcription factors
fig_2a_3 <- build_fig_2a_3()

# known neutrophil genes
fig_2b <- build_fig_2b()

# fig_2c <- build_fig_2c()
fig_2c <- build_fig_2c_updated()

fig_2c_grob <- grid.grabExpr(draw(fig_2c),
                             # set the height according to the actual pdf size
                             # and multiply by relative space occupied by heatmap
                             height = 9.375*4/5,
                             width = 6.75,
                             wrap.grobs = T)

# layout_matrix <- rbind(
#   c(1,2,3,4),
#   c(5,5,5,5),
#   c(5,5,5,5),
#   c(5,5,5,5),
#   c(5,5,5,5)
# )

# pdf("figures/Figure 2.pdf",
#     width = 6.75,
#     height = 9.375)
# grid.arrange(
#   fig_2a_1 %>% add_label("A"),
#   fig_2a_2 %>% add_label(""),
#   fig_2a_3 %>% add_label(""),
#   fig_2b %>% add_label("B"),
#   fig_2c_grob,
#   # arrangeGrob(
#   # fig_2c_grob %>% add_label("C"),
#   # ncol = 1),
#   layout_matrix = layout_matrix
# )
# dev.off()


pdf("figures/Figure 2.pdf", height = 9.375, width = 6.75)
plot_grid(
  plot_grid(fig_2a_1,
            fig_2a_2,
            fig_2a_3,
            fig_2b,
            ncol = 4,
            labels = c("A", "", "", "B"),
            label_size = 8
  ),
  fig_2c_grob,
  ncol=1,
  nrow = 2,
  labels = c("", "C"),
  label_size = 8,
  rel_heights = c(1, 4)
)
dev.off()

# pdf("figures/Figure 2 rownames.pdf",
#     width = 6.75,
#     height = 9.375)
# build_fig_2c_updated()
# dev.off()


############
# FIGURE 3 #
############

# table
fig_3a_2 <- build_fig_3a_2()
fig_3a_2_grob <- grid.grabExpr(draw(fig_3a_2))

# volcano ! currently removed
# fig_3b <- build_fig_3b(return_lgd = FALSE)
# fig_3b_legend <- build_fig_3b(return_lgd = TRUE)



# old S2 -> NEW 3b
fig_3b_outer <-build_fig_s8a()
fig_3b_inset_1 <-build_fig_s2e()
fig_3b_inset_2 <-build_fig_s2f()

fig_3b_outer_legend <- get_legend(fig_3b_outer +
                                    theme(legend.title = element_text(face = "bold", size = 8),
                                          legend.text = element_text(size = 6)) +
                                    guides(color = guide_legend(
                                      override.aes = aes(size = 2.5)
                                    )))
fig_3b_inset_legend <- get_legend(fig_3b_inset_1 +
                                   theme(legend.key.height = unit(2, "mm"),
                                         legend.key.width = unit(4, "mm"),
                                         legend.title = element_text(size = 8),
                                         legend.text = element_text(size = 6)))

fig_3b <- ggdraw(fig_3b_outer +
                   theme(legend.position = "none")) +
  draw_plot(fig_3b_inset_1 +
              theme(legend.position = "none"),
            x = .25,
            y = .3,
            width = .35,
            height = .65) +
  draw_plot(fig_3b_inset_2 +
              theme(legend.position = "none"),
            x = .6,
            y = .3,
            width = .35,
            height = .65) +
  draw_plot_label(
    label = c("", ""),
    x = c(.25, .6),
    y = c(.95, .95),
    size = c(8, 8)
  )

fig_3b_legends <- plot_grid(
  fig_3b_outer_legend,
  fig_3b_inset_legend,
  ncol = 2
)


# lfc heatmap no scaling
fig_3c <- build_fig_3c_1()
fig_3c_grob <- grid.grabExpr(draw(fig_3c),
                             width = 6.75*1,
                             height = 9.375*(1/4))

# expr heatmap
fig_3d_3<- build_fig_3d_3()
fig_3d_3_grob <- grid.grabExpr(draw(fig_3d_3),
                               width = 6.75*(2/3),
                               height = 9.375*(1/4))
# grid.arrange(fig_3d_grob)
# pathway heatmap
fig_3e_2 <- build_fig_3e_2()
fig_3e_2_grob <- grid.grabExpr(draw(fig_3e_2),
                               width = 6.75*(1/3),
                               height = 9.375*(1/4))


# Neutrotime and expr. ! currently removed
# fig_3f <- build_fig_3f()
# fig_3f_1 <- build_fig_3f_1()
# fig_3f_2 <- build_fig_3f_2()



# layout_matrix <- rbind(
#   c(1,1,1,1,2,2),
#   c(1,1,1,1,2,2),
#   c(3,3,3,3,3,3),
#   c(3,3,3,3,3,3),
#   c(4,4,4,5,5,5),
#   c(4,4,4,5,5,5),
#   c(4,4,4,5,5,5),
#   c(4,4,4,5,5,5)

fig_3_cowplot <- plot_grid(
  fig_3a_2_grob,
  plot_grid(
    fig_3b,
    fig_3b_legends,
    nrow = 2,
    rel_heights = c(8, 1)
  ),
  fig_3c_grob,
  plot_grid(
    fig_3d_3_grob,
    fig_3e_2_grob,
    ncol = 2,
    rel_widths = c(2, 1),
    labels = c("D", "E"),
    label_size = 8
    ),
  labels = c("A", "B", "C", ""),
  label_size = 8,
  nrow = 4,
  rel_heights = c(0.8, 1, 1, 1),
  align = "v",
  axis = "lr"
)

save_plot(
  "figures/Figure 3.pdf",
  fig_3_cowplot,
  ncol = 2,
  nrow = 4,
  base_height = 9.375/4,
  base_width = 6.75/2
)


# OLD
# 
# layout_matrix <- rbind(
#   c(1,1,1,1,1,1),
#   c(1,1,1,1,1,1),
#   c(2,2,2,2,2,2),
#   c(2,2,2,2,2,2),
#   c(3,3,3,4,4,4),
#   c(3,3,3,4,4,4),
#   c(3,3,3,4,4,4),
#   c(3,3,3,4,4,4)
# )
# 
# # pdf("figures/figure_3.pdf",
# #     width = 25,
# #     height = 12)
# # fig3b_grob <- fig_3b +
# #   annotation_custom(
# #     grob = fig_3b_legend,
# #     ymin = -700,
# #     xmin = -6
# #   )
# 
# pdf("figures/Figure 3.pdf",
#     width = 6.75,
#     height = 9.375)
# 
# grid.arrange(
#   fig_3a_2_grob %>% add_label("A"),
#   # fig_4b_grob %>% add_label("B"),
#   # fig3b_grob %>% add_label("B"),
#   # placeholder_grob,
#   fig_3c_grob %>% add_label("C"),
#   fig_3d_grob %>% add_label("D"),
#   fig_3e_2_grob %>% add_label("E"),
#   # fig_3f %>% add_label("F"),
#   # fig_3f_1 %>% add_label("F"),
#   # fig_3f_2 %>% add_label(""),
#   # fig_3g %>% add_label("G"),
#   # fig_3h %>% add_label("H"),
#   layout_matrix = layout_matrix
# )
# dev.off()


############
# FIGURE 4 #
############
fig_4a <- build_fig_4a()
fig_4b<- build_fig_4b()
fig_4c <- build_fig_4c()
fig_4d <- build_fig_4d()

# w <- convertWidth(unit(1, "npc")*(1/1), "inch", valueOnly = TRUE)
# h <- convertHeight(unit(1, "npc")*(2/3), "inch", valueOnly = TRUE)

# for grid version
# fig_4d_grob <- grid.grabExpr(draw(fig_4d),
#                              # set the height according to the actual pdf size
#                              # and multiply by relative space occupied by heatmap
#                              # add a correction factor to adjust for add_label
#                              height = 6.75*1.34*2/3,
#                              width = 9.375)

# for cowplot version (no correction factor)
fig_4d_grob <- grid.grabExpr(draw(fig_4d),
                             # set the height according to the actual pdf size
                             # and multiply by relative space occupied by heatmap
                             height = 9.375*1/2,
                             width = 6.75
)

# get around some very annoying alignment problems
aligned_plots <- cowplot::align_plots(fig_4a, fig_4c + theme(legend.position = "none"), align = "v", axis = "l")

top_block <- plot_grid(
  plot_grid(aligned_plots[[1]],
            fig_4b,
            ncol = 2,
            labels = c("A", "B"),
            label_size = 8,
            rel_widths = c(1, 1),
            rel_heights = c(1)),
  aligned_plots[[2]],
  #fig_4d_grob,
  ncol = 1,
  # rel_heights = c(1, 1, 2),
  labels = c("", "C"),
  label_size = 8
  # align = "v",
  # axis = "r"
)


fig_4_cowplot <- plot_grid(
  plot_grid(
    top_block,
    plot_grid(NULL,
              get_legend(fig_4c),
              ncol = 1),
    ncol = 2,
    rel_widths = c(0.8, 0.2),
    label_size = 8
  ),
  fig_4d_grob,
  labels = c("", "D"),
  label_size = 8,
  nrow = 2,
  rel_heights = c(1, 1),
  align = "v",
  axis = "lr"
)

save_plot(
  "figures/Figure 4.pdf",
  fig_4_cowplot,
  ncol = 2,
  nrow = 3,
  base_height = 9.375/3,
  base_width = 6.75/2
)


# layout_matrix <-
#   rbind(
#     c(1,1,2,2,3,3,3),
#     c(1,1,2,2,3,3,3),
#     c(4,4,4,4,4,4,4),
#     c(4,4,4,4,4,4,4),
#     c(4,4,4,4,4,4,4),
#     c(4,4,4,4,4,4,4)
#   )

# cairo_pdf("figures/Figure 4.pdf",
#     height = 9.375,
#     width = 6.75)
# 
# grid.arrange(
#   fig_4a %>% add_label("A"),
#   fig_4b %>% add_label("B"),
#   fig_4c %>% add_label("C"),
#   fig_4d_grob %>% add_label("D"),
#   layout_matrix = layout_matrix
# )
# dev.off()

############
# FIGURE 5 #
############
# fig_5b <- build_fig_5b(return_lgd = FALSE)
# fig_5b_grob <- grid.grabExpr(draw(fig_5b))
# fig_5b_legend <- build_fig_5b(return_lgd = TRUE)
# 
# fig_5a <- build_fig_5a()
# fig_5a_legend <- build_fig_5a(return_lgd = TRUE)
# 
# fig_5a_full <- arrangeGrob(fig_5a, fig_5a_legend, ncol = 1, heights = c(5,1))
# 
# fig_s6a <- build_fig_S6a()
# fig_s6a_1 <- fig_s6a[[1]]
# 
# fig_5b <- build_fig_5b()
# 
# fig_5c <- build_fig_5c(ratio = 0.5)
# fig_5c_legend <- build_fig_5c(return_lgd = TRUE)
# 
# fig_5e <- build_fig_5e()
# 
# fig_5d <- build_fig_5d()
# fig_5d_grob <- grid.grabExpr(draw(fig_5d))
# 
# 
# # fig_5_row2 <- arrangeGrob(fig_5c %>% add_label("C"), fig_5d %>% add_label("D"), nrow = 1, widths = c(1,3))
# # fig_5_row2 <- arrangeGrob(fig_5c %>% add_label("C"), fig_5d_1 %>% add_label("D"), fig_5d_2 %>% add_label(""), fig_5d_3 %>% add_label(""), fig_5d_4 %>% add_label(""), ncol = 3)
# 
# fig_5_row2 <- arrangeGrob(fig_5c %>% add_label("C"), fig_5d_grob %>% add_label("D"), nrow = 1, widths = c(1,1))
# 
# 
# fig_5e <- build_fig_5e()
# fig_5e_full <- arrangeGrob(
#   fig_5e[[1]],
#   fig_5e[[2]],
#   fig_5e[[3]],
#   fig_5e[[4]]
#   
# )
# 
# 
# grDevices::cairo_pdf("figures/Figure 5.pdf",
#                      width = 6.75,
#                      height = 9.375)
# layout_matrix <- rbind(
#   c(1,1,2,2),
#   c(3,3,3,3),
#   c(4,4,4,4),
#   c(4,4,4,4)
# )
# grid.arrange(
#   fig_5a_full %>% add_label("A"),
#   fig_5b %>% add_label("B"),
#   fig_5_row2,
#   fig_5e_full %>% add_label("E"),
#   
#   layout_matrix = layout_matrix
# )
# dev.off()


##################
# Reworked Fig 5 #
##################

# 5A: kept from original
fig_5a <- build_fig_5a()
fig_5a_legend <- build_fig_5a(return_lgd = TRUE)

# 5B taken from supp: S7A heatmap
fig_s6a <- build_fig_S6a(
  width = 6.75/3,
  height = 9.375/6
)
fig_s6a_2 <- fig_s6a[[1]]

# 5C: moved old 5B to C
fig_5b <- build_fig_5b()


# 5D: rows for each comparison; mixing supp and current 5
# col 1
# fig_5d_c1r1 <- build_fig_5c() +
#   theme(title = element_blank()) +
#   ggtitle("")
# fig_5d_c1r2 <- build_fig_s6c_1()
# fig_5d_c1r3 <- build_fig_s6c_2()
# fig_5d_c1r4 <- build_fig_s6c_3()
# fig_5d_c1r5 <- build_fig_s6c_4()

fig_5d_c1r1 <- build_fig_5d_by_row(1)[[1]]
fig_5d_c1r2 <- build_fig_5d_by_row(2)[[1]]
fig_5d_c1r3 <- build_fig_5d_by_row(3)[[1]]
fig_5d_c1r4 <- build_fig_5d_by_row(4)[[1]]
fig_5d_c1r5 <- build_fig_5d_by_row(5)[[1]]


fig_5d_c2r1 <- build_fig_5d_by_row(1)[[2]]
fig_5d_c2r2 <- build_fig_5d_by_row(2)[[2]]
fig_5d_c2r3 <- build_fig_5d_by_row(3)[[2]]
fig_5d_c2r4 <- build_fig_5d_by_row(4)[[2]]
fig_5d_c2r5 <- build_fig_5d_by_row(5)[[2]]


# col 2
fig_s6b <- build_fig_s6b()

fig_5d_c3r1 <- fig_s6b[[1]]
fig_5d_c3r2 <- fig_s6b[[2]]
fig_5d_c3r3 <- fig_s6b[[3]]
fig_5d_c3r4 <- fig_s6b[[4]]
fig_5d_c3r5 <- fig_s6b[[5]]

fig_5_cowplot <- plot_grid(
  plot_grid(
    fig_5a,
    fig_s6a_2,
    fig_5b,
    ncol = 3,
    rel_widths = c(1, 1, 1),
    labels = c("A", "B", "C"),
    label_size = 8
    ),
  plot_grid(
    fig_5d_c1r1,
    fig_5d_c2r1,
    fig_5d_c3r1,
    align = "h",
    axis = "tb",
    ncol = 3,
    rel_widths = c(1, 1, 1)
  ),
  plot_grid(
    fig_5d_c1r2,
    fig_5d_c2r2,
    fig_5d_c3r2,
    align = "h",
    axis = "tb",
    ncol = 3,
    rel_widths = c(1, 1, 1)
  ),
  plot_grid(
    fig_5d_c1r3,
    fig_5d_c2r3,
    fig_5d_c3r3,
    align = "h",
    axis = "tb",
    ncol = 3,
    rel_widths = c(1, 1, 1)
  ),
  plot_grid(
    fig_5d_c1r4,
    fig_5d_c2r4,
    fig_5d_c3r4,
    align = "h",
    axis = "tb",
    ncol = 3,
    rel_widths = c(1, 1, 1)
  ),
  plot_grid(
    fig_5d_c1r5,
    fig_5d_c2r5,
    fig_5d_c3r5,
    align = "h",
    axis = "tb",
    ncol = 3,
    rel_widths = c(1, 1, 1)
  ),
  nrow = 6,
  rel_heights = c(
    2,
    1,
    1,
    1,
    1,
    1
    ),
  labels = c(
    "",
    "D",
    "",
    "",
    "",
    ""
    ),
  label_size = 8
)

save_plot("figures/Figure 5.pdf",
          fig_5_cowplot,
          ncol = 3,
          nrow = 6,
          base_width = 6.75/3,
          base_height = 9.375/6)

############
# FIGURE 6 #
############

fig_6a <-build_fig_6a()
fig_6a_grob <- grid.grabExpr(draw(fig_6a),
                             width = 1.2/2.2*6.75)

# fig_6b_grob <- grid.text("Placeholder")

fig_6a2 <- ggdraw() + 
draw_image("figures/figure_5b.png")

fig_6_row_1 <- plot_grid(
  fig_6a_grob, 
  fig_6a2,
  labels = c("A", "B"),
  label_size = 8,
  rel_widths = c(1.2,1))

fig_6b_1 <- build_fig_6b_1()
fig_6b_2 <- build_fig_6b_2()

fig_6b <- plot_grid(
  fig_6b_1,
  fig_6b_2,
  nrow = 2
)

fig_6c_conditions <- build_fig_6c_conditions()
fig_6c_species <- build_fig_6c_species()

fig_6c_arrange <- plot_grid(
  fig_6c_conditions +
    theme(legend.position = "none"),
  fig_6c_species +
    theme(legend.position = "none"),
  nrow = 2,
  rel_heights = c(1, 1)
)





fig_6d <- build_fig_6d()

fig_6c_conditions_legend <- get_legend(fig_6c_conditions +
                                         theme(
                                           legend.title = element_text(size = 8),
                                           legend.text = element_text(size = 6)
                                         )
                                       )
fig_6c_species_legend <- get_legend(fig_6c_species +
                                      theme(
                                        legend.title = element_text(size = 8),
                                        legend.text = element_text(size = 6)
                                      )
                                    )
fig_6d_legend <- get_legend(fig_6d +
                              theme(
                                legend.title = element_text(size = 8),
                                legend.text = element_text(size = 6)
                              )
                            )

fig_6cd_legends <- plot_grid(
  fig_6c_conditions_legend,
  fig_6c_species_legend,
  fig_6d_legend,
  ncol = 3,
  rel_widths = c(1, 0.66, 1)
)

fig_6 <- plot_grid(
  # fig_6a_grob,
  # fig_6a2,
  fig_6_row_1,
  fig_6b,
  plot_grid(
    fig_6c_arrange,
    fig_6d +
      theme(legend.position = "none"),
    # plot_grid(
    #   fig_6d +
    #     theme(legend.position = "none"),
    #   fig_6cd_legends,
    #   rel_heights = c(6,1),
    #   nrow = 2
    # ),
    ncol = 2,
    rel_widths = c(1,3),
    labels = c("D", "E"),
    label_size = 8
  ),
  fig_6cd_legends,
  ncol = 1,
  # rel_heights = c(1.31,0.9,1,1.25,1.3),
  #rel_heights = c(1.2,0.8,1.25,1.2),
  rel_heights = c(1,1,1,0.11),
  labels = c("", "C", ""),
  label_size = 8
)


save_plot("figures/Figure 6.pdf",
          fig_6,
          ncol = 1,
          nrow = 5,
          base_width = 6.75/1,
          base_height = 9.375/5)

# layout_matrix <- rbind(
#   c(1,1),
#   c(2,2),
#   c(3,3),
#   c(4,4),
#   c(5,6),
#   c(7,8)
# )
# 
# plot_grid(
#   fig_6a_grob,
#   fig_6b_grob,
#   fig_6c_1,
#   fig_6c_2,
#   plot_grid(fig_6d_1,
#             fig_6d_2,
#             ncol = 2),
#   plot_grid(fig_6e_1,
#             fig_6e_2,
#             ncol = 2),
#   ncol = 1,
#   labels = c("A", "B", "C", NULL, NULL, NULL)
# )
# 
# 
# pdf("figures/figure_6.pdf",
#     height = 15,
#     width = 10)
# 
# grid.arrange(
#   fig_6a_grob,
#   fig_6b_grob,
#   fig_6c_1,
#   fig_6c_2,
#   fig_6d_1,
#   fig_6d_2,
#   fig_6e_1,
#   fig_6e_2,
#   layout_matrix = layout_matrix
# )
# dev.off()


############
# FIGURE 7 #
############

### 3D plots (Fig. 7A) must be manually inserted via Illustrator

fig_7A <- build_fig_7a()
fig_7B <- build_fig_7b()
fig_7C <- build_fig_7c()

pdf("figures/Figure 7.pdf", height = 9.375, width = 6.75)
cowplot::plot_grid(
  fig_7A,
  fig_7B,
  fig_7C,
  ncol=1,
  nrow = 3,
  labels = c("A", "B", "C"),
  label_size = 8,
  rel_widths = c(1,1,1),
  rel_heights = c(0.9, 0.7, 1),
  align = "hv",
  axis = "rlbt",
  scale = c(0.9, 0.9, 1)
)
dev.off()


#########################
#########################
# SUPPLEMENTARY FIGURES #
#########################
#########################


#############
# FIGURE S1 #
#############

fig_s1_1 <- build_fig_s1_1()

pdf("figures/Figure S1.pdf", width = 6.75, height = 9.375)
fig_s1_1
dev.off()

# fig_s1_3 <- build_fig_s1_3()

# pdf("figures/Figure S5.pdf", width = 6.75, height = 9.375)
# fig_s5_3
# dev.off()


#############
# FIGURE S2 #
#############

fig_s2a <- build_fig_s2a()
fig_s2b <- build_fig_s2b()
fig_s2c <- build_fig_s2c()

# old S8
# fig_s2d <-build_fig_s8a()
# fig_s2e <-build_fig_s2e()
# fig_s2f <-build_fig_s2f()
# 
# fig_s2d_legend <- get_legend(fig_s2d +
#                                theme(legend.title = element_text(face = "bold")))
# fig_s2ef_legend <- get_legend(fig_s2e +
#                                 theme(legend.key.height = unit(4, "mm")))
# 
# fig_s2def_inset <- ggdraw(fig_s2d +
#                             theme(legend.position = "none")) +
#   draw_plot(fig_s2e +
#               theme(legend.position = "none"),
#             x = .25,
#             y = .3,
#             width = .35,
#             height = .65) +
#   draw_plot(fig_s2f +
#               theme(legend.position = "none"),
#             x = .6,
#             y = .3,
#             width = .35,
#             height = .65) +
#   draw_plot_label(
#     label = c("", ""),
#     x = c(.25, .6),
#     y = c(.95, .95),
#     size = c(8, 8)
#   )
# 
# legends <- plot_grid(
#   fig_s2d_legend,
#   fig_s2ef_legend,
#   ncol = 2
# )

fig_s2 <- plot_grid(
  fig_s2a,
  # fig_s2def_inset,
  # legends,
  plot_grid(
    fig_s2b,
    fig_s2c,
    labels = c("B", "C"),
    label_size = 8,
    rel_widths = c(1, 2),
    ncol = 2,
    # TODO: fix alignment
    align = "h",
    axis = "tb"
    ),
  labels = c("A", ""),
  label_size = 8,
  nrow = 2,
  rel_heights = c(1, 1)
)

save_plot("figures/Figure S2.pdf",
          fig_s2,
          ncol = 2,
          nrow = 2,
          base_width = 6.75/2,
          base_height = 9.375/4
)

# OLD 
# layout_matrix <- rbind(
#   c(1,1),
#   c(2,3)
# )
# pdf("figures/Figure S2.pdf",
#     width = 6.75,
#     height = 9.375)
# 
# grid.arrange(
#   fig_s2a %>% add_label("A"),
#   fig_s2b %>% add_label("B"),
#   fig_s2c %>% add_label("C"),
#   layout_matrix = layout_matrix
# )
# dev.off()



#############
# FIGURE S3 #
#############

fig_s3a <- build_fig_s3a()
fig_s3b <- build_fig_s3b()
fig_s3c <- build_fig_s3c()

# pdf("figures/Figure S6.pdf",
#     height = 9.375,
#     width = 6.75)

fig_s3 <- plot_grid(
  fig_s3a,
  plot_grid(
    fig_s3b,
    fig_s3c,
    labels = c("B", "C"),
    label_size = 8,
    ncol = 2,
    rel_widths = c(1, 3),
    align = "h",
    axis = "tb"
  ),
  labels = c("A", ""),
  label_size = 8,
  nrow = 2,
  rel_heights = c(1, 1)
)

save_plot("figures/Figure S3.pdf",
          fig_s3,
          ncol = 2,
          nrow = 2,
          base_width = 6.75/2,
          base_height = 9.375/2
)


#############
# FIGURE S4 #
#############

fig_s4a <- build_fig_s4a()
fig_s4b <- build_fig_s4b()
fig_s4c <- build_fig_s4c()

# pdf("figures/Figure S6.pdf",
#     height = 9.375,
#     width = 6.75)

fig_s4 <- plot_grid(
  fig_s4a,
  fig_s4b,
  fig_s4c,
  labels = c("A", "B", "C"),
  label_size = 8,
  nrow = 3,
  align = "v",
  axis = "lr"
)

save_plot("figures/Figure S4.pdf",
          fig_s4,
          ncol = 1,
          nrow = 3,
          base_width = 6.75/1,
          base_height = 9.375/3
)


#############
# FIGURE S5 #
#############

fig_s7a <-build_fig_s7a()
fig_s7a_grob <- grid.grabExpr(draw(fig_s7a, merge_legend = TRUE),
                              height = 9.375/2)

fig_s7b <-build_fig_s7b()
fig_s7b_grob <- grid.grabExpr(draw(fig_s7b))

fig_s7 <- plot_grid(
  fig_s7b_grob,
  fig_s7a_grob,
  labels = c("A", "B"),
  label_size = 8,
  nrow = 2,
  rel_heights = c(1, 1),
  align = "v",
  axis = "l"
)

save_plot("figures/Figure S5.pdf",
          fig_s7,
          ncol = 1,
          nrow = 2,
          base_width = 6.75/1,
          base_height = 9.375/2
)


#############
# FIGURE S6 #
#############
fig_s5a <-build_fig_s5a()
fig_s5b <-build_fig_s5b()

fig_s5 <- plot_grid(
  fig_s5a,
  fig_s5b,
  labels = c("A", "B"),
  label_size = 8,
  nrow = 2,
  rel_heights = c(1, 1),
  align = "v",
  axis = "l"
)

save_plot("figures/Figure S6.pdf",
          fig_s5,
          ncol = 1,
          nrow = 2,
          base_width = 6.75/1,
          base_height = 9.375/2
)

# #############
# # FIGURE S7
# #############
# 
# # fig_s6a <- build_fig_S6a()
# # fig_s6a_1 <- fig_s6a[[1]]
# 
# # m <- rbind(
# #   c(NA,1,NA)
# # )
# 
# fig_s6a_2 <- fig_s6a[[2]]
# 
# # Plot full list 6a
# fig_s6a_rownames <- build_fig_S6a(full_list = TRUE)
# pdf("figures/figure_S7a_rownames.pdf", height = 15)
# grid.arrange(fig_s6a_rownames[[1]])
# dev.off()
# 
# fig_s6a_full <- grid.arrange(
#   fig_s6a_2,
#   fig_s6a_1,
#   ncol = 1
# )
# 
# 
# 
# fig_s6b <- build_fig_s6b()
# fig_s6b_full <- 
#   arrangeGrob(
#     fig_s6b[[1]],
#     fig_s6b[[2]],
#     fig_s6b[[3]],
#     fig_s6b[[4]],
#     fig_s6b[[5]],
#     ncol = 3
#   )
# 
# fig_s6_top <- arrangeGrob(fig_s6a_full %>% add_label("A"), fig_s6b_full %>% add_label("B"), nrow = 1, widths = c(1,2))
# 
# 
# fig_s6c_1 <- build_fig_s6c_1(ratio = 0.5)
# fig_s6c_2 <- build_fig_s6c_2(ratio = 0.5)
# fig_s6c_3 <- build_fig_s6c_3(ratio = 0.5)
# fig_s6c_4 <- build_fig_s6c_4(ratio = 0.5)
# 
# fig_s6c <- arrangeGrob(fig_s6c_1, fig_s6c_2, fig_s6c_3, fig_s6c_4, ncol = 2)
# fig_s6c_full <- arrangeGrob(fig_s6c, fig_5c_legend, heights = c(10,1))
# 
# grDevices::cairo_pdf("figures/Figure S7.pdf",
#                      width = 6.75,
#                      height = 9.375)
# layout_matrix <- rbind(
#   c(1,1,1,1),
#   c(2,2,2,2),
#   c(2,2,2,2)
# )
# grid.arrange(
#   fig_s6_top %>% add_label(""),
#   fig_s6c_full %>% add_label("C"),
#   layout_matrix = layout_matrix
# )
# dev.off()

#############
# FIGURE S7 #
#############

fig_s8_1 <- build_fig_s8_1()
fig_s8_2 <- build_fig_s8_2()

fig_s8 <- plot_grid(
  fig_s8_1,
  fig_s8_2,
  nrow = 2
)
pdf("figures/Figure S7.pdf", height = 4, width = 6.75)
fig_s8
dev.off()




#############
# FIGURE SX #
#############

# Draw top 10 motifs from AP_VS_BL_increase motif enrichment
motifs <- plot_grid(
  ggdraw() + draw_image("data/raw_data/zymosan/atac/AP_VS_BL_increase/knownResults/known1.logo.svg"),
  ggdraw() + draw_image("data/raw_data/zymosan/atac/AP_VS_BL_increase/knownResults/known2.logo.svg"),
  ggdraw() + draw_image("data/raw_data/zymosan/atac/AP_VS_BL_increase/knownResults/known3.logo.svg"),
  ggdraw() + draw_image("data/raw_data/zymosan/atac/AP_VS_BL_increase/knownResults/known4.logo.svg"),
  ggdraw() + draw_image("data/raw_data/zymosan/atac/AP_VS_BL_increase/knownResults/known5.logo.svg"),
  ggdraw() + draw_image("data/raw_data/zymosan/atac/AP_VS_BL_increase/knownResults/known6.logo.svg"),
  ggdraw() + draw_image("data/raw_data/zymosan/atac/AP_VS_BL_increase/knownResults/known7.logo.svg"),
  ggdraw() + draw_image("data/raw_data/zymosan/atac/AP_VS_BL_increase/knownResults/known8.logo.svg"),
  ggdraw() + draw_image("data/raw_data/zymosan/atac/AP_VS_BL_increase/knownResults/known9.logo.svg"),
  ggdraw() + draw_image("data/raw_data/zymosan/atac/AP_VS_BL_increase/knownResults/known10.logo.svg"),
  ncol = 1
)

#############
# FIGURE S8 #
#############

# fig_s8a <-build_fig_s8a()
# fig_s8b <-build_fig_s8b()
# 
# fig_s8 <- plot_grid(
#   fig_s8a,
#   fig_s8b,
#   labels = c("A", "B"),
#   label_size = 8,
#   nrow = 2,
#   rel_heights = c(1, 1),
#   align = "v",
#   axis = "l"
# )
# 
# save_plot("figures/Figure S8.pdf",
#           fig_s8,
#           ncol = 1,
#           nrow = 2,
#           base_width = 6.75/1,
#           base_height = 9.375/2
# )

########################
########################
# SUPPLEMENTARY TABLES #
########################
########################

build_table_S5()
