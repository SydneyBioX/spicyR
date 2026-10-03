suppressPackageStartupMessages({ library(spicyR); library(SpatialDatasets); library(ggplot2); library(grid) })
spe <- spe_Ali_2020()
cells <- as.data.frame(SummarizedExperiment::colData(spe))
cells <- cells[cells$ER.Status %in% c("neg", "pos"), ]
cells$ER <- factor(ifelse(cells$ER.Status == "pos", "ER+", "ER-"), levels = c("ER-", "ER+"))
res <- spicy(cells, condition = "ER", subject = "metabricId", r = 25, imageID = "file_id", cellType = "description",
             spatialCoords = c("Location_Center_X", "Location_Center_Y"))      # every pair, as in the vignette
## left: a 160 x 160 µm window of an ER+ core: the proliferating tumour cells (to), filled when a T cell (from) lies
## within 25 µm, with the circle drawn around two of them. Pair T cells -> HR- Ki67+.
f <- "HR- Ki67+"; t <- "T cells"
cand <- aggregate(cbind(nf = description == f, nt = description == t) ~ file_id, data = cells[cells$ER == "ER+", ], FUN = sum)
img <- cand$file_id[order(-pmin(cand$nf, cand$nt))][1]
z <- cells[cells$file_id == img, ]; z$x <- z$Location_Center_X; z$y <- z$Location_Center_Y
zf <- z[z$description == f, ]; zt <- z[z$description == t, ]
zf$near <- vapply(seq_len(nrow(zf)), function(i) any(sqrt((zt$x - zf$x[i])^2 + (zt$y - zf$y[i])^2) <= 25), TRUE)
w <- 80
# a window with 8 to 16 tumour cells, as evenly split as possible between those with and without a T cell nearby
inner <- function(i) abs(zf$x - zf$x[i]) <= w - 25 & abs(zf$y - zf$y[i]) <= w - 25
score <- vapply(seq_len(nrow(zf)), function(i) { k <- inner(i)
  if (sum(k) < 8 || sum(k) > 16) -Inf else min(sum(zf$near[k]), sum(!zf$near[k])) - 0.01 * sum(k) }, 0)
zf <- zf[inner(which.max(score)), ]
best <- data.frame(x = mean(range(zf$x)), y = mean(range(zf$y)))     # centre the window on those tumour cells
z <- z[abs(z$x - best$x) <= w & abs(z$y - best$y) <= w & z$description != f, ]
z$role <- ifelse(z$description == t, "T cell (from)", "other cells")
zf$role <- ifelse(zf$near, "tumour cell (to), a T cell within 25 µm", "tumour cell (to), none within 25 µm")
lev <- c("tumour cell (to), a T cell within 25 µm", "tumour cell (to), none within 25 µm", "T cell (from)", "other cells")
pts <- rbind(z[, c("x", "y", "role")], zf[, c("x", "y", "role")]); pts$role <- factor(pts$role, levels = lev)
pts <- pts[order(pts$role, decreasing = TRUE), ]
d0 <- (zf$x - best$x)^2 + (zf$y - best$y)^2      # circles around the most central tumour cell of each kind
show <- c(which(zf$near)[which.min(d0[zf$near])], which(!zf$near)[which.min(d0[!zf$near])])
circ <- do.call(rbind, lapply(show, function(i) data.frame(id = i,
  x = zf$x[i] + 25 * cos(seq(0, 2 * pi, length.out = 120)), y = zf$y[i] + 25 * sin(seq(0, 2 * pi, length.out = 120)))))
lab <- zf[show[1], ]
p1 <- ggplot(pts, aes(x, y)) +
  geom_path(data = circ, aes(group = id), colour = "black", linewidth = 0.5, linetype = 2) +
  geom_point(aes(colour = role, fill = role, size = role), shape = 21, stroke = 0.9) +
  annotate("label", x = lab$x + 20, y = lab$y + 27, label = "r = 25 µm", size = 3.2, fill = "white") +
  scale_colour_manual(values = c("#b3261e", "#b3261e", "#1f6fb4", "grey82"), name = NULL, drop = FALSE) +
  scale_fill_manual(values = c("#b3261e", "white", "#1f6fb4", "grey82"), name = NULL, drop = FALSE) +
  scale_size_manual(values = c(2.8, 2.8, 2.4, 1.2), guide = "none") +
  coord_equal() + theme_void(base_size = 11) +
  theme(legend.position = "bottom", legend.direction = "vertical", plot.title = element_text(face = "bold", size = 11)) +
  labs(title = "Does each tumour cell have a T cell\nwithin 25 µm? Compare with chance")
## right: the effect of every patient (the extra fraction of tumour cells next to T cells), by ER status
b <- bind(res, pairName = "T cells__HR- Ki67+"); names(b)[ncol(b)] <- "effect"
b$ER <- factor(b$condition, levels = c("ER-", "ER+"))
p2 <- ggplot(b[is.finite(b$effect), ], aes(ER, effect)) +
  geom_hline(yintercept = 0, colour = "grey55", linetype = 2) +
  geom_jitter(width = 0.18, height = 0, alpha = 0.35, size = 1.1, colour = "grey30") +
  geom_boxplot(fill = NA, outlier.shape = NA, width = 0.5, linewidth = 0.6) +
  theme_classic(base_size = 11) + theme(plot.title = element_text(face = "bold", size = 11)) +
  labs(x = NULL, y = "Extra fraction of tumour cells\nnext to T cells (beyond chance)", title = "Compare patients between groups",
       subtitle = "one point per patient")
png(file.path(Sys.getenv("OUT"), "spicyR_overview.png"), width = 2000, height = 1000, res = 220)
grid.newpage(); pushViewport(viewport(layout = grid.layout(1, 2, widths = unit(c(1.05, 1), "null"))))
print(p1, vp = viewport(layout.pos.row = 1, layout.pos.col = 1)); print(p2, vp = viewport(layout.pos.row = 1, layout.pos.col = 2))
invisible(dev.off()); cat(img, nrow(zf), sum(zf$near), nrow(b), "\n")
