# Draw Quad Venn for Epigenetically Repressed Genes Across Subtypes
# Author: Sean Maden

library(VennDiagram)

draw.quad.venn(area1=69, 
               area2=51, 
               area3=3, 
               area4=0, 
               n12=19, 
               n13=0, 
               n14=0, 
               n23=1, 
               n24=0,
               n34=0, 
               n123=0, 
               n124=0, 
               n134=0, 
               n234=0, 
               n1234=0, 
               category = c("HM (N=69)",
                            "IM (N=51)",
                            "LM (N=3)",
                            "MM (N=0)"), 
               lwd = rep(0, 4), 
               lty = rep(0,4), 
               col = rep(NA,4), 
               fill = c("yellow",
                        "coral",
                        "gray",
                        "lightblue"), 
               alpha = rep(0.5, 4),
               label.col = rep("black", 15), 
               cex = rep(1, 15),
               fontface = rep("plain", 15), 
               fontfamily = rep("",15), 
               cat.pos = c(-15, 15, 0, 0), 
               cat.dist = c(0.22, 0.22, 0.11, 0.11), 
               cat.col = rep("black", 4), 
               cat.cex = rep(1, 4), 
               cat.fontface = rep("plain", 4),
               cat.fontfamily = rep("", 4), 
               cat.just =rep(list(c(0.5, 0.5)), 4), 
               rotation.degree = 0,
               rotation.centre = c(0.5, 0.5), ind = TRUE, 
               cex.prop =NULL, 
               print.mode = "raw", 
               sigdigs = 3, 
               direct.area =FALSE, 
               area.vector = 0)

dev.off()

#
