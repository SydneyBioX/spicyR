model_data_gen <- function(spe, ref, target, rad, window = "convex", weights = NULL) {
  coords <- spatialCoords(spe)
  cell_type <- colData(spe)$cellType
  image_id <- colData(spe)$imageID
  condition <- colData(spe)$condition
  subject <- colData(spe)$subject
  #cell_id <- paste0(cell_type, seq_len(nrow(coords)))
  
  result_list <- list()
  unique_imgs <- unique(image_id)
  img_lookup <- setNames(seq_along(unique(image_id)), as.character(unique(image_id)))
  
  for (img in unique_imgs) {
    weight_img <- if (!is.null(weights)) weights[which(unique_imgs == img)] else NA
    idx_img <- which(image_id == img)
    subject_img <- subject[idx_img][1]
    condition_img <- condition[idx_img][1]
    coords_img <- coords[idx_img, ]
    types_img <- cell_type[idx_img]
    type_counts <- ave(seq_along(types_img), types_img, FUN = seq)
    cell_id_img <- paste0(types_img, "_", img, "_", type_counts)
    
    img_wt_index <- img_lookup[as.character(img)]
    weight_img <- if (!is.null(weights)) weights[img_wt_index] else NA
    
    
    idx_ref <- which(types_img == ref)
    idx_target <- which(types_img == target)
    
    if (length(intersect(idx_ref, idx_target)) > 0) {
      warning(glue::glue("Ref and target overlap in image {img}"))
    }
    
    # Set spatial window
    if (window == "rectangle") {
      x_range <- range(coords_img[, 1])
      y_range <- range(coords_img[, 2])
      win <- owin(xrange = x_range, yrange = y_range)
    } else if (window == "convex") {
      win <- convexhull.xy(coords_img[, 1], coords_img[, 2])
    } else if (window == "concave") {
      hull_coords <- concaveman::concaveman(coords_img)
      win <- owin(poly = list(x = hull_coords[, 1], y = hull_coords[, 2]))
    } else {
      stop("Invalid value for `window`. Use 'rectangle', 'convex', or 'concave'.")
    }
    
    ppp <- spatstat.geom::ppp(x = coords_img[, 1], y = coords_img[, 2], window = win)
    area_img <- spatstat.geom::area.owin(win)
    dens <- (length(idx_target) / area_img) * (pi * rad^2)
      
    if (length(idx_ref) > 0 && length(idx_target) > 0) {
      D <- rdist(
      as.matrix(coords_img[idx_ref, c("x", "y"), drop = FALSE]),
      as.matrix(coords_img[idx_target, c("x", "y"), drop = FALSE]))
      counts_target_for_each_ref <- rowSums(D <= rad)
      
      df_img <- data.frame(
        cellID = as.factor(cell_id_img[idx_ref]),
        imageID = as.factor(img),
        subject = as.factor(subject_img),
        n = counts_target_for_each_ref,
        condition = as.factor(condition_img),
        density = dens,
        weights = weight_img
      )
      
      result_list[[img]] <- df_img
    } else {
      cat(sprintf("No reference or target cells found for image %s\n", img))
    }
  }
  
  final_table <- bind_rows(result_list)
  return(final_table)
}

plotImage <- function(spe, image_id, type1, type2) {
  fspe <- spe[, colData(spe)$imageID == image_id]
  coords <- as.data.frame(fspe@int_colData$spatialCoords)
  celltypes <- fspe@colData$cellType
  
  # Filter for only the two selected cell types
  keep_idx <- which(celltypes %in% c(type1, type2))
  coords <- coords[keep_idx, ]
  coords$celltype <- factor(celltypes[keep_idx])
  
  ggplot(coords, aes(x = x, y = y, color = celltype)) +
    geom_point(size = 1) +
    coord_fixed() +
    #scale_y_reverse() +  # Optional: flip y-axis if image coordinate system
    theme_minimal() +
    labs(
      title = paste("Cells in", image_id, "for", type1, "and", type2),
      x = "x", y = "y"
    )
}

count_cells_per_image <- function(spe, cell_types) {
  # Safety check
  if (!"imageID" %in% colnames(colData(spe))) {
    stop("The 'imageID' column is missing from colData(spe)")
  }
  
  if (!"cellType" %in% colnames(colData(spe))) {
    stop("The 'cellType' column is missing from colData(spe)")
  }
  
  # Extract relevant metadata
  meta <- colData(spe) %>%
    as.data.frame() %>%
    filter(cellType %in% cell_types)
  
  # Count cells
  counts <- meta %>%
    group_by(imageID, cellType) %>%
    summarise(n_cells = n(), .groups = "drop") %>%
    pivot_wider(names_from = cellType, values_from = n_cells, values_fill = 0) %>%
    arrange(imageID)
  
  return(counts)
}


mem_stats_summary <- function(spe, rad, model) {
  coords <- SpatialExperiment::spatialCoords(spe)
  image_ids <- colData(spe)$imageID
  cell_types <- colData(spe)$cellType
  unique_images <- unique(image_ids)
  
  ppp_list <- list()
  
  for (img in unique_images) {
    idx <- which(image_ids == img)
    coords_img <- coords[idx, ]
    marks_img <- factor(cell_types[idx])  # ensure it's a factor
    
    if (nrow(coords_img) < 2 || nlevels(marks_img) < 2) next
    
    win <- spatstat.geom::owin(
      xrange = range(coords_img[, 1]),
      yrange = range(coords_img[, 2])
    )
    
    ppp_list[[img]] <- spatstat.geom::ppp(
      x = coords_img[, 1],
      y = coords_img[, 2],
      window = win,
      marks = marks_img
    )
  }
  
  # Compute Kcross and Lcross
  k_results <- lapply(ppp_list, function(pp) {
    if (spatstat.geom::npoints(pp) >= 2) 
      spatstat.explore::Kcross(pp, i = "A", j = "B", correction = "Ripley") 
    else NULL
  })
  
  l_results <- lapply(ppp_list, function(pp) {
    if (spatstat.geom::npoints(pp) >= 2) 
      spatstat.explore::Lcross(pp, i = "A", j = "B", correction = "Ripley")
    else NULL
  })
  
  # Extract K/L at radius `rad`
  extract_iso_at_r <- function(x, target_r) {
    if (is.null(x)) return(NA)
    idx <- which.min(abs(x$r - target_r))
    return(x$iso[idx])
  }
  
  k_at_r <- unlist(lapply(k_results, extract_iso_at_r, target_r = rad))
  l_at_r <- unlist(lapply(l_results, extract_iso_at_r, target_r = rad))
  
  # Correct against null expectations
  k_corrected <- k_at_r - (pi * rad^2)
  l_corrected <- l_at_r - rad
  
  # Vector of image IDs in your desired order (e.g. from l_results)
  ordered_ids <- names(l_results)
  # Reorder ranef(pm_s)$imageID to match that order
  ranef <- ranef(model)[["cond"]][["imageID"]][ordered_ids, , drop = FALSE]
  
  comp_data <- data.frame(
    imageID = names(ppp_list),
    K = k_corrected,
    L = l_corrected,
    R = ranef
  )
  
  return(comp_data)
}

simulate_AB <- function(structure, 
                        A_count, 
                        B_count, 
                        window, 
                        width = NA, 
                        sigma = NA){
  
  xrange <- spatstat.geom::as.rectangle(window)$xrange
  yrange <- spatstat.geom::as.rectangle(window)$yrange
  xwidth <- diff(xrange)
  yheight <- diff(yrange)
  
  if (structure == "Random") {
    A <- spatstat.random::rpoispp(A_count / (xwidth * yheight), win = window)
    B <- spatstat.random::rpoispp(B_count / (xwidth * yheight), win = window)
  } 
  
  else if (structure == "DependentCluster") {
    if (is.null(sigma) || sigma <= 0 || is.na(sigma)) {
      stop(glue::glue("Invalid sigma value: {sigma}"))
    }
    
    A <- spatstat.random::rpoispp(A_count / (xwidth * yheight), win = window)
    aDens <- spatstat.explore::density.ppp(A, sigma = sigma, kernel = "disc")
    aDens$v <- pmax(aDens$v, 0) * B_count / A_count
    B <- spatstat.random::rpoispp(aDens)
  } 
  
  else if (structure == "Overlap") {
    if (is.null(width)) width <- 0.1
    gap <- width * xwidth / 2
    
    A_win <- spatstat.geom::owin(xrange = c(xrange[1], xrange[2] - gap),
                                 yrange = yrange)
    B_win <- spatstat.geom::owin(xrange = c(xrange[1] + gap, xrange[2]),
                                 yrange = yrange)
    
    A <- spatstat.random::rpoispp(A_count/spatstat.geom::area.owin(A_win), win = A_win)
    B <- spatstat.random::rpoispp(B_count/spatstat.geom::area.owin(B_win), win = B_win)
  } 
  
  else if (structure == "Corner") {
    if (is.null(width)) width <- 0
    
    corner_size <- 1 - width
    A_win <- spatstat.geom::owin(
      xrange = c(xrange[1], xrange[1] + corner_size * xwidth),
      yrange = c(yrange[1], yrange[1] + corner_size * yheight)
    )
    B_win <- spatstat.geom::owin(
      xrange = c(xrange[2] - corner_size * xwidth, xrange[2]),
      yrange = c(yrange[2] - corner_size * yheight, yrange[2])
    )
    
    #cat(glue::glue( "[DEBUG] width: {width}, corner_size: {corner_size}, A_area: {spatstat.geom::area.owin(A_win)}, B_area: {spatstat.geom::area.owin(B_win)}\n" ))
    
    A <- spatstat.random::rpoispp(A_count/spatstat.geom::area.owin(A_win), win = A_win)
    B <- spatstat.random::rpoispp(B_count/spatstat.geom::area.owin(B_win), win = B_win)
    #cat(glue::glue("[DEBUG] A points: {A$n}, B points: {B$n}\n"))
    
  }  else {
    stop("Unsupported structure type: ", structure)
  }
  
  
  return(list(A = A, B = B))
}

build_spe <- function(i, 
                      counts, 
                      nPatients, 
                      nIm, 
                      window, 
                      structure1, 
                      structure2, 
                      width1 = NA, 
                      width2 = NA, 
                      lambda = NA, 
                      delta = NA) {
  
  set.seed(i)
  
  g1 <- rpois(nPatients/2, lambda)
  g2 <- rpois(nPatients/2, lambda + delta)
  adjustSigma <- c(g1, g2) + 1
  #cat("[DEBUG] lambda =", lambda, " | delta =", delta, "\n")
  #cat("[DEBUG] adjustSigma (length =", length(adjustSigma), "):\n")
  #print(adjustSigma)
  
  # Build table of (patient, image) combinations and assign groups
  combo_df <- expand.grid(
    patient = 1:nPatients,
    image = 1:nIm
  ) |>
    dplyr::mutate(
      condition = ifelse(patient <= nPatients / 2, "Group1", "Group2"),
      structure = ifelse(condition == "Group1", structure1, structure2),
      width = ifelse(condition == "Group1", width1, width2)
    )
  
  # Simulate per (patient, image)
  result_list <- purrr::pmap(combo_df, function(patient, image, condition, structure, width) {
    # DEBUG: check type of patient  
    #cat("[DEBUG] class(patient):", class(patient),
    #    "| value:", patient, 
    #    "| typeof:", typeof(patient), "\n")
    
    sCount1 <- sample(counts, 1)
    sCount2 <- sample(counts, 1)
    
    # Determine sigma only if needed
    sigma <- if (structure == "DependentCluster") {
      adjustSigma[patient]
    } else {
      NULL
    }
    
    AB <- simulate_AB(
      structure = structure,
      A_count = sCount1,
      B_count = sCount2,
      window = window,
      width = width,
      sigma = sigma
    )
    
    a <- AB$A
    b <- AB$B
    
    list(points = data.frame(
      x = c(a$x, b$x),
      y = c(a$y, b$y),
      cellType = rep(c("A", "B"), c(a$n, b$n)),
      imageID = rep(paste(patient, image, sep = "_"), a$n + b$n)
    ),
    subject = patient,
    condition = condition)
  })
  
  # Combine results into final data frame
  points_df <- dplyr::bind_rows(purrr::map(result_list, "points"))
  
  phenoData <- dplyr::bind_rows(purrr::map(result_list, function(r) {
    data.frame(imageID = unique(r$points$imageID), 
               subject = r$subject, 
               condition = r$condition)
  }))
  
  cellExp <- dplyr::left_join(points_df, phenoData, by = "imageID")
  
  # Create SpatialExperiment object
  setClass("ExpData", contains = "environment")
  
  spe <- SpatialExperiment::SpatialExperiment(
    assays = list(counts = matrix(0, nrow = 1, ncol = nrow(cellExp))),  # dummy assay
    spatialCoords = cbind(x = cellExp$x, y = cellExp$y),
    colData = cellExp[, c("cellType", "imageID", "subject", "condition")]
  )
  #print(i)
  return(spe)
}