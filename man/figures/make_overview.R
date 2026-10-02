suppressPackageStartupMessages({ library(spicyR); library(SpatialDatasets); library(ggplot2); library(grid) })
spe <- spe_Ali_2020()
cells <- as.data.frame(SummarizedExperiment::colData(spe))
cells <- cells[cells$ER.Status %in% c("neg", "pos"), ]
cells$ER <- factor(ifelse(cells$ER.Status == "pos", "ER+", "ER-"), levels = c("ER-", "ER+"))
res <- spicy(cells, condition = "ER", subject = "metabricId", r = 25, imageID = "file_id", cellType = "description",
             spatialCoords = c("Location_Center_X", "Location_Center_Y"), from = "HR- Ki67+", to = "T cells")
## left: a 160 x 160 µm window of an ER+ core with several proliferating tumour cells and T cells near them
f <- "HR- Ki67+"; t <- "T cells"
cand <- aggregate(cbind(nf = description == f, nt = description == t) ~ file_id, data = cells[cells$ER == "ER+", ], FUN = sum)
img <- cand$file_id[order(-pmin(cand$nf, cand$nt))][1]
z <- cells[cells$file_id == img, ]; z$x <- z$Location_Center_X; z$y <- z$Location_Center_Y
zf <- z[z$description == f, ]; best <- NULL; bestn <- -1
for (i in seq_len(nrow(zf))) { d <- sqrt((z$x - zf$x[i])^2 + (z$y - zf$y[i])^2); n <- sum(z$description == t & d <= 25)
  if (n > bestn) { bestn <- n; best <- zf[i, ] } }
w <- 80; z <- z[abs(z$x - best$x) <= w & abs(z$y - best$y) <= w, ]
z$role <- ifelse(z$description == f, "proliferating HR- tumour cell (from)", ifelse(z$description == t, "T cell (to)", "other cells"))
z$role <- factor(z$role, levels = c("proliferating HR- tumour cell (from)", "T cell (to)", "other cells"))
circ <- data.frame(x = best$x + 25 * cos(seq(0, 2 * pi, length.out = 200)), y = best$y + 25 * sin(seq(0, 2 * pi, length.out = 200)))
pal <- c("proliferating HR- tumour cell (from)" = "#b3261e", "T cell (to)" = "#1f6fb4", "other cells" = "grey82")
p1 <- ggplot(z, aes(x, y)) +
  geom_point(aes(colour = role, size = role)) +
  geom_path(data = circ, colour = "black", linewidth = 0.5, linetype = 2) +
  annotate("label", x = best$x + 27, y = best$y + 27, label = "r = 25 µm", size = 3.2, label.size = 0, fill = "white") +
  scale_colour_manual(values = pal, name = NULL) + scale_size_manual(values = c(2.6, 2.6, 1.2), guide = "none") +
  coord_equal() + theme_void(base_size = 11) +
  theme(legend.position = "bottom", legend.direction = "vertical", plot.title = element_text(face = "bold", size = 11)) +
  labs(title = "Count the T cells near each tumour cell,\nand compare with chance")
## right: the excess of every patient, by ER status
b <- bind(res); names(b)[ncol(b)] <- "excess"
b$ER <- factor(b$condition, levels = c("ER-", "ER+"))
p2 <- ggplot(b[is.finite(b$excess), ], aes(ER, excess)) +
  geom_hline(yintercept = 0, colour = "grey55", linetype = 2) +
  geom_jitter(width = 0.18, height = 0, alpha = 0.35, size = 1.1, colour = "grey30") +
  geom_boxplot(fill = NA, outlier.shape = NA, width = 0.5, linewidth = 0.6) +
  coord_cartesian(ylim = c(-1.7, 4)) +
  theme_classic(base_size = 11) + theme(plot.title = element_text(face = "bold", size = 11)) +
  labs(x = NULL, y = "Extra T cells per tumour cell\n(beyond chance)", title = "Compare patients between groups",
       subtitle = sprintf("one point per patient (axis cut at 4); p = %.1e", res$cellResults$p_value))
png(file.path(Sys.getenv("OUT"), "spicyR_overview.png"), width = 2000, height = 1000, res = 220)
grid.newpage(); pushViewport(viewport(layout = grid.layout(1, 2, widths = unit(c(1.05, 1), "null"))))
print(p1, vp = viewport(layout.pos.row = 1, layout.pos.col = 1)); print(p2, vp = viewport(layout.pos.row = 1, layout.pos.col = 2))
invisible(dev.off()); cat(img, bestn, nrow(b), res$cellResults$p_value, "\n")
