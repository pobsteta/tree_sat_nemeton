#!/usr/bin/env Rscript
# ==============================================================================
# TreeSatAI-Time-Series — Visualisation
# Profils phénologiques, cartes de confusion, importance des variables
# ==============================================================================

# Configuration chargée via le package
# ggplot2 et patchwork sont des dépendances du package (Imports)

# --- Palette de couleurs adaptative ------------------------------------------

# Génère une palette nommée par species_name (fonctionne avec noms latins,
# français, ou genres). Réutilise SPECIES_COLORS si les noms correspondent,
# sinon génère une palette qualitative.
.get_species_palette <- function(species_names) {
  species_names <- sort(unique(species_names))
  n <- length(species_names)

  # Tenter de réutiliser SPECIES_COLORS (clés = noms latins)
  if (all(species_names %in% names(SPECIES_COLORS))) {
    return(SPECIES_COLORS[species_names])
  }

  # Tenter correspondance via noms français
  french_match <- match(species_names, SPECIES$french)
  if (!any(is.na(french_match))) {
    cols <- SPECIES_COLORS[SPECIES$latin[french_match]]
    names(cols) <- species_names
    return(cols)
  }

  # Fallback : palette qualitative large
  if (n <= 12) {
    pal <- RColorBrewer::brewer.pal(max(3, n), "Set3")[seq_len(n)]
  } else {
    pal <- grDevices::hcl.colors(n, palette = "Dark 3")
  }
  names(pal) <- species_names
  pal
}

# --- Profils phénologiques par espèce ----------------------------------------

#' Tracé des profils phénologiques NDVI pour toutes les espèces
#' @param ts_long Dataset long avec colonnes date, species_name, NDVI
#' @param index_name Nom de l'indice à tracer (défaut "NDVI")
#' @param save_path Chemin pour sauvegarder le graphique (NULL = pas de sauvegarde)
#' @return Objet ggplot
plot_phenological_profiles <- function(ts_long, index_name = "NDVI",
                                        save_path = NULL) {
  # Calcul du profil moyen par espèce et date
  profiles <- ts_long |>
    dplyr::group_by(species_name, date) |>
    dplyr::summarise(
      mean_val = mean(.data[[index_name]], na.rm = TRUE),
      sd_val   = sd(.data[[index_name]], na.rm = TRUE),
      .groups = "drop"
    ) |>
    dplyr::mutate(doy = as.numeric(format(date, "%j")))

  sp_colors <- .get_species_palette(profiles$species_name)
  n_species <- length(sp_colors)

  p <- ggplot(profiles, aes(x = doy, y = mean_val, color = species_name)) +
    geom_ribbon(aes(ymin = mean_val - sd_val, ymax = mean_val + sd_val,
                    fill = species_name), alpha = 0.1, color = NA) +
    geom_line(linewidth = 0.8) +
    scale_color_manual(values = sp_colors, name = "Esp\u00e8ce") +
    scale_fill_manual(values = sp_colors, guide = "none") +
    scale_x_continuous(
      breaks = c(1, 32, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335),
      labels = c("Jan", "F\u00e9v", "Mar", "Avr", "Mai", "Jun",
                 "Jul", "Ao\u00fb", "Sep", "Oct", "Nov", "D\u00e9c")
    ) +
    labs(
      title = paste("Profils ph\u00e9nologiques \u2014", index_name),
      subtitle = glue::glue("Moyenne \u00b1 \u00e9cart-type par esp\u00e8ce ({n_species} classes)"),
      x = "Mois",
      y = index_name
    ) +
    theme_minimal(base_size = VIS_PARAMS$font_size) +
    theme(
      legend.position = "right",
      legend.text = element_text(face = "italic", size = 7),
      legend.key.height = unit(0.4, "cm"),
      plot.title = element_text(face = "bold")
    )

  if (!is.null(save_path)) {
    ggsave(save_path, p, width = VIS_PARAMS$width_cm, height = VIS_PARAMS$height_cm,
           units = "cm", dpi = VIS_PARAMS$dpi)
    log_msg("Profils ph\u00e9nologiques sauvegard\u00e9s : {save_path}", level = "success")
  }

  p
}

#' Profils par type : feuillus vs résineux
#' @param ts_long Dataset long
#' @param save_path Chemin de sauvegarde
#' @return Objet ggplot
plot_profiles_by_type <- function(ts_long, save_path = NULL) {
  # Ajouter le type d'essence — tenter la jointure, sinon deviner
  species_info <- SPECIES |> dplyr::select(french, type, phenologie)
  ts_typed <- ts_long |>
    dplyr::left_join(species_info, by = c("species_name" = "french"))

  # Si la jointure n'a rien donné (noms de genres), deviner le type
  if (all(is.na(ts_typed$type))) {
    coniferes <- c("Pin sylvestre", "\u00c9pic\u00e9a commun", "Sapin pectin\u00e9",
                    "Douglas", "M\u00e9l\u00e8ze d'Europe")
    ts_typed$type <- ifelse(ts_typed$species_name %in% coniferes, "r\u00e9sineux", "feuillu")
    caducs <- c("H\u00eatre", "Ch\u00eane p\u00e9doncul\u00e9", "Bouleau verruqueux",
                "Fr\u00eane commun", "\u00c9rable sycomore", "Charme", "Peupliers",
                "Aulne", "Tilleul", "Prunus", "M\u00e9l\u00e8ze d'Europe", "Cleared")
    ts_typed$phenologie <- ifelse(ts_typed$species_name %in% caducs, "caducifoli\u00e9", "sempervirent")
  }

  profiles <- ts_typed |>
    dplyr::group_by(species_name, type, phenologie, date) |>
    dplyr::summarise(mean_ndvi = mean(NDVI, na.rm = TRUE), .groups = "drop") |>
    dplyr::mutate(doy = as.numeric(format(date, "%j")))

  sp_colors <- .get_species_palette(profiles$species_name)

  p <- ggplot(profiles, aes(x = doy, y = mean_ndvi, color = species_name)) +
    geom_line(linewidth = 0.7, alpha = 0.8) +
    facet_grid(type ~ phenologie, scales = "free_y") +
    scale_color_manual(values = sp_colors, name = "Esp\u00e8ce") +
    scale_x_continuous(
      breaks = c(1, 91, 182, 274),
      labels = c("Jan", "Avr", "Jul", "Oct")
    ) +
    labs(
      title = "Profils ph\u00e9nologiques par type et r\u00e9gime foliaire",
      x = "Mois", y = "NDVI"
    ) +
    theme_minimal(base_size = VIS_PARAMS$font_size) +
    theme(
      legend.position = "bottom",
      legend.text = element_text(face = "italic", size = 6),
      strip.text = element_text(face = "bold")
    ) +
    guides(color = guide_legend(ncol = 5))

  if (!is.null(save_path)) {
    ggsave(save_path, p, width = 28, height = 20, units = "cm", dpi = VIS_PARAMS$dpi)
    log_msg("Profils par type sauvegard\u00e9s : {save_path}", level = "success")
  }

  p
}

# --- Matrice de confusion ----------------------------------------------------

#' Tracé d'une matrice de confusion avec ggplot2
#' @param conf_matrix Objet table ou matrix de confusion
#' @param class_names Noms des classes
#' @param title Titre du graphique
#' @param save_path Chemin de sauvegarde
#' @return Objet ggplot
plot_confusion_matrix <- function(conf_matrix, class_names = NULL,
                                   title = "Matrice de confusion",
                                   save_path = NULL) {
  if (is.null(class_names)) {
    class_names <- sort(unique(c(rownames(conf_matrix), colnames(conf_matrix))))
  }

  # S'assurer que la matrice couvre toutes les classes (certaines peuvent manquer)
  full_mat <- matrix(0L, nrow = length(class_names), ncol = length(class_names),
                     dimnames = list(class_names, class_names))
  common_rows <- intersect(rownames(conf_matrix), class_names)
  common_cols <- intersect(colnames(conf_matrix), class_names)
  full_mat[common_rows, common_cols] <- conf_matrix[common_rows, common_cols]

  # Normaliser par ligne (recall par classe)
  row_sums <- rowSums(full_mat)
  conf_norm <- sweep(full_mat, 1, ifelse(row_sums == 0, 1, row_sums), "/")
  conf_norm[is.nan(conf_norm)] <- 0

  # Conversion en data.frame long
  conf_df <- expand.grid(
    Reference  = class_names,
    Prediction = class_names,
    stringsAsFactors = FALSE
  )
  conf_df$Count    <- as.vector(full_mat)
  conf_df$Percent  <- as.vector(conf_norm) * 100

  # Labels
  conf_df$label <- ifelse(
    conf_df$Count > 0,
    paste0(conf_df$Count, "\n(", round(conf_df$Percent, 1), "%)"),
    ""
  )

  # OA sur les classes présentes
  diag_sum <- sum(diag(full_mat))
  total    <- sum(full_mat)
  oa_pct   <- if (total > 0) round(diag_sum / total * 100, 1) else 0

  p <- ggplot(conf_df, aes(x = Prediction, y = Reference, fill = Percent)) +
    geom_tile(color = "white", linewidth = 0.5) +
    geom_text(aes(label = label), size = 2.2) +
    scale_fill_gradientn(
      colors = c("white", "#deebf7", "#9ecae1", "#3182bd", "#08519c"),
      limits = c(0, 100),
      name = "Rappel (%)"
    ) +
    scale_x_discrete(position = "top") +
    labs(
      title = title,
      subtitle = paste0("OA = ", oa_pct, "%"),
      x = "Pr\u00e9diction",
      y = "R\u00e9f\u00e9rence"
    ) +
    theme_minimal(base_size = 8) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 0, face = "italic", size = 6),
      axis.text.y = element_text(face = "italic", size = 6),
      plot.title = element_text(face = "bold"),
      legend.position = "right"
    ) +
    coord_fixed()

  if (!is.null(save_path)) {
    ggsave(save_path, p, width = 25, height = 22, units = "cm", dpi = VIS_PARAMS$dpi)
    log_msg("Matrice de confusion sauvegardée : {save_path}", level = "success")
  }

  p
}

# --- Importance des variables ------------------------------------------------

#' Tracé de l'importance des variables (Random Forest)
#' @param importance_df data.frame avec colonnes variable, importance
#' @param top_n Nombre de variables à afficher
#' @param title Titre
#' @param save_path Chemin de sauvegarde
#' @return Objet ggplot
plot_variable_importance <- function(importance_df, top_n = 30,
                                      title = "Importance des variables",
                                      save_path = NULL) {
  # Top N
  imp_top <- importance_df |>
    dplyr::arrange(dplyr::desc(importance)) |>
    dplyr::slice_head(n = top_n)

  # Détecter le type de feature
  imp_top <- imp_top |>
    dplyr::mutate(
      feature_type = dplyr::case_when(
        grepl("^pheno_|_pheno_", variable)   ~ "Phénologie",
        grepl("^NDVI|^EVI|^NDWI|^CRI|^NBR", variable) ~ "Indice spectral",
        grepl("^fourier_|_fourier_", variable) ~ "Fourier",
        grepl("^B\\d|^B8A", variable)          ~ "Bande spectrale",
        TRUE                                    ~ "Autre"
      )
    )

  type_colors <- c(
    "Phénologie"       = "#e41a1c",
    "Indice spectral"  = "#377eb8",
    "Fourier"          = "#4daf4a",
    "Bande spectrale"  = "#984ea3",
    "Autre"            = "#999999"
  )

  p <- ggplot(imp_top, aes(x = reorder(variable, importance),
                            y = importance, fill = feature_type)) +
    geom_col(alpha = 0.85) +
    coord_flip() +
    scale_fill_manual(values = type_colors, name = "Type de feature") +
    labs(
      title = title,
      subtitle = paste0("Top ", top_n, " variables (Random Forest MDA)"),
      x = NULL,
      y = "Mean Decrease Accuracy"
    ) +
    theme_minimal(base_size = VIS_PARAMS$font_size) +
    theme(
      axis.text.y = element_text(size = 7),
      plot.title = element_text(face = "bold"),
      legend.position = "bottom"
    )

  if (!is.null(save_path)) {
    ggsave(save_path, p, width = VIS_PARAMS$width_cm, height = 20,
           units = "cm", dpi = VIS_PARAMS$dpi)
    log_msg("Importance des variables sauvegardée : {save_path}", level = "success")
  }

  p
}

# --- Heatmap des signatures spectrales ---------------------------------------

#' Heatmap des profils spectraux moyens par espèce
#' @param feature_matrix Matrice de features (1 ligne par parcelle)
#' @param save_path Chemin de sauvegarde
#' @return Objet pheatmap (invisible)
plot_spectral_heatmap <- function(feature_matrix, save_path = NULL) {
  # Extraire les colonnes NDVI temporelles
  ndvi_cols <- grep("^NDVI_d\\d+$", names(feature_matrix), value = TRUE)

  if (length(ndvi_cols) == 0) {
    cli::cli_alert_warning("Aucune colonne NDVI temporelle trouvée")
    return(invisible(NULL))
  }

  # Moyenne par espèce
  heatmap_data <- feature_matrix |>
    dplyr::group_by(species_name) |>
    dplyr::summarise(dplyr::across(dplyr::all_of(ndvi_cols), \(x) mean(x, na.rm = TRUE)),
                     .groups = "drop") |>
    tibble::column_to_rownames("species_name")

  # Noms des colonnes en DOY
  doy_labels <- gsub("NDVI_d", "DOY ", ndvi_cols)
  colnames(heatmap_data) <- doy_labels

  # Normaliser par ligne (z-score)
  heatmap_scaled <- t(scale(t(as.matrix(heatmap_data))))

  if (!is.null(save_path)) {
    png(save_path, width = VIS_PARAMS$width_cm * 40, height = VIS_PARAMS$height_cm * 40,
        res = VIS_PARAMS$dpi)
  }

  result <- pheatmap::pheatmap(
    heatmap_scaled,
    cluster_cols = FALSE,
    cluster_rows = TRUE,
    color = viridis::viridis(100),
    main = glue::glue("Signatures ph\u00e9nologiques NDVI \u2014 {nrow(heatmap_data)} classes (z-score)"),
    fontsize = 8,
    fontsize_row = 7,
    angle_col = 45,
    border_color = NA
  )

  if (!is.null(save_path)) {
    dev.off()
    log_msg("Heatmap sauvegardée : {save_path}", level = "success")
  }

  invisible(result)
}

# --- Comparaison phénologique caduc vs persistant ----------------------------

#' Graphique comparatif caduc vs sempervirent
#' @param ts_long Dataset long
#' @param save_path Chemin de sauvegarde
#' @return Objet ggplot
plot_deciduous_vs_evergreen <- function(ts_long, save_path = NULL) {
  ts_typed <- ts_long |>
    dplyr::left_join(
      SPECIES |> dplyr::select(french, phenologie),
      by = c("species_name" = "french")
    )

  # Si la jointure n'a rien donné (noms de genres), deviner la phénologie
  if (all(is.na(ts_typed$phenologie))) {
    caducs <- c("H\u00eatre", "Ch\u00eane p\u00e9doncul\u00e9", "Bouleau verruqueux",
                "Fr\u00eane commun", "\u00c9rable sycomore", "Charme", "Peupliers",
                "Aulne", "Tilleul", "Prunus", "M\u00e9l\u00e8ze d'Europe", "Cleared")
    ts_typed$phenologie <- ifelse(ts_typed$species_name %in% caducs,
                                   "caducifoli\u00e9", "sempervirent")
  }

  ts_typed <- ts_typed |>
    dplyr::mutate(doy = as.numeric(format(date, "%j")))

  # Profils moyens par type
  profiles <- ts_typed |>
    dplyr::group_by(phenologie, doy) |>
    dplyr::summarise(
      mean_ndvi = mean(NDVI, na.rm = TRUE),
      q25 = quantile(NDVI, 0.25, na.rm = TRUE),
      q75 = quantile(NDVI, 0.75, na.rm = TRUE),
      .groups = "drop"
    )

  phenol_colors <- c("caducifolié" = "#e66101", "sempervirent" = "#1a9641")

  p <- ggplot(profiles, aes(x = doy, y = mean_ndvi, color = phenologie)) +
    geom_ribbon(aes(ymin = q25, ymax = q75, fill = phenologie),
                alpha = 0.15, color = NA) +
    geom_line(linewidth = 1.2) +
    scale_color_manual(values = phenol_colors,
                       labels = c("Caducifolié", "Sempervirent"),
                       name = "Type phénologique") +
    scale_fill_manual(values = phenol_colors, guide = "none") +
    scale_x_continuous(
      breaks = c(1, 32, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335),
      labels = c("J", "F", "M", "A", "M", "J", "J", "A", "S", "O", "N", "D")
    ) +
    annotate("rect", xmin = 60, xmax = 150, ymin = -Inf, ymax = Inf,
             alpha = 0.05, fill = "#4daf4a") +
    annotate("text", x = 105, y = Inf, label = "Verdissement",
             vjust = 1.5, size = 3, color = "#4daf4a") +
    annotate("rect", xmin = 250, xmax = 335, ymin = -Inf, ymax = Inf,
             alpha = 0.05, fill = "#e66101") +
    annotate("text", x = 292, y = Inf, label = "Sénescence",
             vjust = 1.5, size = 3, color = "#e66101") +
    labs(
      title = "Caducifolié vs Sempervirent — Signatures NDVI annuelles",
      subtitle = "Médiane et intervalle interquartile (Q25-Q75)",
      x = "Mois", y = "NDVI"
    ) +
    theme_minimal(base_size = VIS_PARAMS$font_size) +
    theme(
      plot.title = element_text(face = "bold"),
      legend.position = "top"
    )

  if (!is.null(save_path)) {
    ggsave(save_path, p, width = VIS_PARAMS$width_cm, height = 12,
           units = "cm", dpi = VIS_PARAMS$dpi)
    log_msg("Comparaison caduc/persistant sauvegardée : {save_path}", level = "success")
  }

  p
}

# --- Métriques par espèce (barplot) ------------------------------------------

#' Barplot des métriques de classification par espèce
#' @param metrics_df data.frame avec colonnes species, precision, recall, f1
#' @param save_path Chemin de sauvegarde
#' @return Objet ggplot
plot_species_metrics <- function(metrics_df, save_path = NULL) {
  metrics_long <- metrics_df |>
    tidyr::pivot_longer(cols = c(precision, recall, f1),
                        names_to = "metric", values_to = "value") |>
    dplyr::mutate(metric = factor(metric,
      levels = c("precision", "recall", "f1"),
      labels = c("Précision", "Rappel", "F1-Score")))

  p <- ggplot(metrics_long, aes(x = reorder(species, -value), y = value, fill = metric)) +
    geom_col(position = position_dodge(width = 0.8), alpha = 0.85, width = 0.7) +
    scale_fill_manual(values = c("#2166ac", "#b2182b", "#4daf4a"), name = "Métrique") +
    geom_hline(yintercept = 0.8, linetype = "dashed", color = "grey50", alpha = 0.7) +
    labs(
      title = "Performance de classification par espèce",
      subtitle = "Précision, Rappel et F1-Score",
      x = NULL,
      y = "Score"
    ) +
    theme_minimal(base_size = VIS_PARAMS$font_size) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, face = "italic", size = 7),
      plot.title = element_text(face = "bold"),
      legend.position = "top"
    ) +
    coord_cartesian(ylim = c(0, 1))

  if (!is.null(save_path)) {
    ggsave(save_path, p, width = VIS_PARAMS$width_cm, height = 12,
           units = "cm", dpi = VIS_PARAMS$dpi)
    log_msg("Métriques par espèce sauvegardées : {save_path}", level = "success")
  }

  p
}

# --- Dashboard récapitulatif ------------------------------------------------

#' Création d'un dashboard récapitulatif (4 panneaux)
#' @param plots Liste de ggplots (pheno_profiles, confusion, importance, metrics)
#' @param save_path Chemin de sauvegarde
#' @return Objet patchwork
plot_dashboard <- function(plots, save_path = NULL) {
  dashboard <- (plots$profiles + plots$confusion) /
               (plots$importance + plots$metrics) +
    patchwork::plot_annotation(
      title = "TreeSatAI-Time-Series — Tableau de bord de classification",
      subtitle = glue::glue("{nrow(SPECIES_GROUPS_INFO)} classes foresti\u00e8res \u2014 Sentinel-1/2 s\u00e9ries temporelles"),
      theme = theme(
        plot.title = element_text(size = 14, face = "bold"),
        plot.subtitle = element_text(size = 10)
      )
    )

  if (!is.null(save_path)) {
    ggsave(save_path, dashboard, width = 40, height = 30, units = "cm", dpi = VIS_PARAMS$dpi)
    log_msg("Dashboard sauvegardé : {save_path}", level = "success")
  }

  dashboard
}

# --- Rapport cartographique PDF -----------------------------------------------

#' Construction d'un composite True Color Sentinel-2 (RGB)
#'
#' S\u00e9lectionne une date estivale (proche du 1er juillet) et empile
#' B04 (Rouge), B03 (Vert), B02 (Bleu) pour produire un fond satellite.
#'
#' @param cube_list Liste de SpatRasters S2 (une par date)
#' @param dates Vecteur de dates correspondant aux \u00e9l\u00e9ments de cube_list
#' @return SpatRaster 3 bandes (Red, Green, Blue) ou NULL si \u00e9chec
#' @export
build_s2_rgb_composite <- function(cube_list, dates = NULL) {
  if (is.null(cube_list) || length(cube_list) == 0) return(NULL)

  # Choisir une date estivale (DOY 152-243 = juin-ao\u00fbt)
  if (!is.null(dates) && length(dates) == length(cube_list)) {
    doy <- as.numeric(format(as.Date(dates), "%j"))
    summer_idx <- which(doy >= 152 & doy <= 243)
    if (length(summer_idx) > 0) {
      best <- summer_idx[which.min(abs(doy[summer_idx] - 182))]
    } else {
      best <- ceiling(length(cube_list) / 2)
    }
  } else {
    best <- ceiling(length(cube_list) / 2)
  }

  scene <- cube_list[[best]]
  bn <- names(scene)
  b04 <- grep("^B04", bn, value = TRUE)[1]
  b03 <- grep("^B03", bn, value = TRUE)[1]
  b02 <- grep("^B02", bn, value = TRUE)[1]

  if (is.na(b04) || is.na(b03) || is.na(b02)) {
    log_msg("Bandes RGB introuvables dans le cube S2", level = "warning")
    return(NULL)
  }

  rgb <- c(scene[[b04]], scene[[b03]], scene[[b02]])
  names(rgb) <- c("Red", "Green", "Blue")
  rgb
}

#' Palette de couleurs pour les classes d'essences
#' @param class_names Vecteur de noms de classes
#' @return Vecteur nomm\u00e9 de couleurs hex
#' @keywords internal
.match_species_colors <- function(class_names) {
  n <- length(class_names)
  fallback <- grDevices::hcl.colors(n, "Dark 3")
  colors <- character(n)
  for (i in seq_len(n)) {
    nm <- class_names[i]
    if (nm %in% names(SPECIES_GROUP_COLORS)) {
      colors[i] <- SPECIES_GROUP_COLORS[[nm]]
    } else {
      idx <- match(nm, SPECIES$french)
      if (!is.na(idx)) {
        colors[i] <- SPECIES_COLORS[idx]
      } else {
        colors[i] <- fallback[i]
      }
    }
  }
  names(colors) <- class_names
  colors
}

#' G\u00e9n\u00e9ration d'un rapport cartographique PDF multi-pages
#'
#' Produit un PDF A4 paysage contenant :
#' \itemize{
#'   \item Carte des essences foresti\u00e8res (fond satellite optionnel,
#'         coupes/vides transparents)
#'   \item Entropie de Shannon (diversit\u00e9 / taux de m\u00e9lange)
#'   \item Carte de confiance du mod\u00e8le
#'   \item Carte feuillus / r\u00e9sineux
#'   \item Richesse sp\u00e9cifique (nombre d'essences par pixel)
#'   \item Composition foresti\u00e8re (diagramme en barres)
#'   \item Confiance moyenne par essence
#' }
#'
#' @param rasters Liste issue de build_species_raster() :
#'   species, confidence, shannon, probas, presence, legend
#' @param statistics data.frame issu de compute_map_statistics()
#' @param output_dir R\u00e9pertoire de sortie
#' @param s2_rgb SpatRaster 3 bandes (R, G, B) pour le fond satellite (NULL = pas de fond)
#' @param aoi sf object \u2014 contour de la zone d'int\u00e9r\u00eat (NULL = pas de contour)
#' @return Chemin vers le fichier PDF (invisible)
#' @export
generate_prediction_report_pdf <- function(rasters, statistics, output_dir,
                                            s2_rgb = NULL, aoi = NULL) {
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  pdf_path <- file.path(output_dir, "rapport_cartographique.pdf")
  log_msg("G\u00e9n\u00e9ration du rapport cartographique PDF...")

  # --- M\u00e9tadonn\u00e9es ---
  legend_df <- rasters$legend
  class_names <- legend_df$species
  n_classes <- length(class_names)
  species_colors <- .match_species_colors(class_names)

  # Identifier Coupe/Vide / Cleared
  cleared_idx <- which(class_names %in% c("Coupe/Vide", "Cleared"))

  # Raster d'essences : Coupe/Vide \u2192 NA (transparent)
  r_species_clean <- rasters$species
  if (length(cleared_idx) > 0) {
    rcl_na <- cbind(cleared_idx, rep(NA_real_, length(cleared_idx)))
    r_species_clean <- terra::classify(r_species_clean, rcl_na)
  }
  # Supprimer les niveaux cat\u00e9goriels pour un trac\u00e9 num\u00e9rique propre
  if (!is.null(terra::levels(r_species_clean)[[1]])) {
    levels(r_species_clean) <- NULL
  }

  # Couleurs/noms pour les classes pr\u00e9sentes dans le raster
  present_vals <- sort(unique(terra::values(r_species_clean)))
  present_vals <- present_vals[!is.na(present_vals)]
  plot_colors <- species_colors[present_vals]
  plot_names  <- class_names[present_vals]

  # Projeter l'AOI
  aoi_proj <- NULL
  if (!is.null(aoi)) {
    aoi_proj <- sf::st_transform(aoi, terra::crs(rasters$species))
  }

  # --- Ouvrir le PDF (A4 paysage) ---
  grDevices::pdf(pdf_path, width = 11.69, height = 8.27,
                  title = "TreeSatAI Nemeton \u2014 Rapport cartographique")
  on.exit(grDevices::dev.off(), add = TRUE)

  # ================================================================
  # PAGE 1 : Carte des essences foresti\u00e8res
  # ================================================================
  tryCatch({
    if (!is.null(s2_rgb)) {
      # Fond satellite Sentinel-2 True Color
      terra::plotRGB(s2_rgb, stretch = "lin",
                      mar = c(3.5, 2, 3.5, 9),
                      axes = TRUE,
                      main = "Carte des essences foresti\u00e8res")
      terra::plot(r_species_clean, type = "classes",
                  col = plot_colors,
                  alpha = 0.7,
                  legend = FALSE, axes = FALSE,
                  add = TRUE)
    } else {
      par(mar = c(3.5, 2, 3.5, 9))
      terra::plot(r_species_clean, type = "classes",
                  col = plot_colors,
                  main = "Carte des essences foresti\u00e8res",
                  legend = FALSE, axes = TRUE)
    }
    if (!is.null(aoi_proj)) {
      plot(sf::st_geometry(aoi_proj), add = TRUE,
           border = "white", lwd = 2, lty = 2)
    }
    # L\u00e9gende manuelle
    par(xpd = TRUE)
    legend("right", inset = c(-0.13, 0),
           legend = plot_names, fill = plot_colors,
           cex = 0.6, bty = "n",
           title = "Essences", title.font = 2, border = NA)
    par(xpd = FALSE)
    bg_label <- if (!is.null(s2_rgb)) "Fond : Sentinel-2 True Color" else ""
    mtext(paste0(bg_label, " | Coupes/Vides = transparents"),
          side = 1, line = 2.2, cex = 0.7, adj = 1, font = 3)
  }, error = function(e) {
    plot.new()
    text(0.5, 0.5, paste("Erreur carte essences :", e$message), cex = 0.8)
  })

  # ================================================================
  # PAGE 2 : Entropie de Shannon
  # ================================================================
  if (!is.null(rasters$shannon)) {
    tryCatch({
      par(mar = c(3.5, 2, 3.5, 6))
      terra::plot(rasters$shannon,
                  col = viridis::viridis(100),
                  main = "Entropie de Shannon \u2014 Taux de m\u00e9lange des essences",
                  axes = TRUE,
                  plg = list(title = "Shannon\n(0\u20131)",
                             title.cex = 0.8, cex = 0.7))
      if (!is.null(aoi_proj)) {
        plot(sf::st_geometry(aoi_proj), add = TRUE,
             border = "white", lwd = 2, lty = 2)
      }
      mtext(paste0("0 = peuplement pur (100% une essence) | ",
                    "1 = m\u00e9lange maximal (toutes essences \u00e9quiprobables)"),
            side = 1, line = 2.2, cex = 0.7, font = 3)
    }, error = function(e) {
      plot.new()
      text(0.5, 0.5, paste("Erreur carte Shannon :", e$message), cex = 0.8)
    })
  }

  # ================================================================
  # PAGE 3 : Carte de confiance
  # ================================================================
  tryCatch({
    par(mar = c(3.5, 2, 3.5, 6))
    conf_pal <- grDevices::colorRampPalette(
      c("#d73027", "#fc8d59", "#fee08b", "#d9ef8b", "#91cf60", "#1a9850")
    )(100)
    terra::plot(rasters$confidence,
                col = conf_pal,
                main = "Carte de confiance \u2014 Probabilit\u00e9 de l'essence pr\u00e9dite",
                axes = TRUE,
                range = c(0, 1),
                plg = list(title = "Probabilit\u00e9",
                           title.cex = 0.8, cex = 0.7))
    if (!is.null(aoi_proj)) {
      plot(sf::st_geometry(aoi_proj), add = TRUE,
           border = "black", lwd = 2, lty = 2)
    }
    mtext(paste0("Rouge = faible confiance (zones de m\u00e9lange/confusion) | ",
                  "Vert = forte confiance (peuplement monosp\u00e9cifique)"),
          side = 1, line = 2.2, cex = 0.65, font = 3)
  }, error = function(e) {
    plot.new()
    text(0.5, 0.5, paste("Erreur carte confiance :", e$message), cex = 0.8)
  })

  # ================================================================
  # PAGE 4 : Feuillus / R\u00e9sineux
  # ================================================================
  tryCatch({
    type_info <- data.frame(species = class_names,
                             code = seq_along(class_names),
                             stringsAsFactors = FALSE)
    type_source <- if (all(class_names %in% SPECIES_GROUPS_INFO$group)) {
      SPECIES_GROUPS_INFO[, c("group", "type")]
    } else {
      data.frame(group = SPECIES$french, type = SPECIES$type,
                 stringsAsFactors = FALSE)
    }
    type_info <- merge(type_info, type_source,
                        by.x = "species", by.y = "group", all.x = TRUE)
    type_info$type_code <- ifelse(type_info$type == "feuillu", 1L,
                            ifelse(type_info$type == "r\u00e9sineux", 2L,
                                   NA_integer_))
    rcl_type <- as.matrix(type_info[order(type_info$code),
                                     c("code", "type_code")])
    r_type <- rasters$species
    if (!is.null(terra::levels(r_type)[[1]])) levels(r_type) <- NULL
    r_type <- terra::classify(r_type, rcl_type)

    type_cols <- c("#66c2a5", "#1b7837")
    par(mar = c(3.5, 2, 3.5, 8))
    terra::plot(r_type, type = "classes",
                col = type_cols,
                main = "Types forestiers \u2014 Feuillus / R\u00e9sineux",
                legend = FALSE, axes = TRUE)
    if (!is.null(aoi_proj)) {
      plot(sf::st_geometry(aoi_proj), add = TRUE,
           border = "black", lwd = 2, lty = 2)
    }
    par(xpd = TRUE)
    legend("right", inset = c(-0.08, 0),
           legend = c("Feuillus", "R\u00e9sineux"), fill = type_cols,
           cex = 0.8, bty = "n", title = "Type", title.font = 2)
    par(xpd = FALSE)
    type_stats <- aggregate(surface_ha ~ type, data = statistics, sum, na.rm = TRUE)
    type_stats <- type_stats[type_stats$type %in% c("feuillu", "r\u00e9sineux"), ]
    if (nrow(type_stats) > 0) {
      txt <- paste(
        ifelse(type_stats$type == "feuillu", "Feuillus", "R\u00e9sineux"),
        ":", type_stats$surface_ha, "ha", collapse = "  |  ")
      mtext(txt, side = 1, line = 2.2, cex = 0.8)
    }
  }, error = function(e) {
    plot.new()
    text(0.5, 0.5, paste("Erreur carte types :", e$message), cex = 0.8)
  })

  # ================================================================
  # PAGE 5 : Richesse sp\u00e9cifique (nombre d'essences par pixel)
  # ================================================================
  if (!is.null(rasters$presence)) {
    tryCatch({
      r_richness <- terra::app(rasters$presence, sum, na.rm = TRUE)
      names(r_richness) <- "richesse_specifique"
      max_rich <- max(terra::values(r_richness), na.rm = TRUE)

      par(mar = c(3.5, 2, 3.5, 6))
      terra::plot(r_richness,
                  col = viridis::magma(max(max_rich, 2)),
                  main = "Richesse sp\u00e9cifique \u2014 Nombre d'essences d\u00e9tect\u00e9es par pixel",
                  axes = TRUE,
                  plg = list(title = "Nb essences",
                             title.cex = 0.8, cex = 0.7))
      if (!is.null(aoi_proj)) {
        plot(sf::st_geometry(aoi_proj), add = TRUE,
             border = "white", lwd = 2, lty = 2)
      }
      threshold <- CLASSIF_PARAMS$presence_threshold
      if (is.null(threshold)) threshold <- 0.10
      mtext(paste0("Seuil de pr\u00e9sence : probabilit\u00e9 \u2265 ",
                    threshold * 100, "%"),
            side = 1, line = 2.2, cex = 0.7, font = 3)
    }, error = function(e) {
      plot.new()
      text(0.5, 0.5, paste("Erreur carte richesse :", e$message), cex = 0.8)
    })
  }

  # ================================================================
  # PAGE 6 : Composition foresti\u00e8re (ggplot2)
  # ================================================================
  tryCatch({
    detected <- statistics[statistics$n_pixels > 0, ]
    detected_sp <- detected[!detected$espece %in% c("Coupe/Vide", "Cleared"), ]
    if (nrow(detected_sp) > 0) {
      p <- ggplot(detected_sp,
                   aes(x = reorder(espece, surface_ha),
                       y = surface_ha, fill = espece)) +
        geom_col(alpha = 0.85, show.legend = FALSE) +
        geom_text(aes(label = paste0(pct, "% \u2014 ", surface_ha, " ha")),
                  hjust = -0.05, size = 3) +
        coord_flip(clip = "off") +
        scale_fill_manual(values = species_colors) +
        scale_y_continuous(expand = expansion(mult = c(0, 0.35))) +
        labs(
          title = "Composition foresti\u00e8re \u2014 Surface par essence",
          subtitle = paste0("Surface totale class\u00e9e : ",
                            sum(detected_sp$surface_ha), " ha (",
                            nrow(detected_sp), " essences d\u00e9tect\u00e9es)"),
          x = NULL, y = "Surface (ha)"
        ) +
        theme_minimal(base_size = 11) +
        theme(
          plot.title = element_text(face = "bold", size = 14),
          plot.subtitle = element_text(size = 10, color = "grey40"),
          axis.text.y = element_text(face = "italic", size = 9),
          panel.grid.major.y = element_blank()
        )
      print(p)
    }
  }, error = function(e) {
    plot.new()
    text(0.5, 0.5, paste("Erreur statistiques :", e$message), cex = 0.8)
  })

  # ================================================================
  # PAGE 7 : Confiance moyenne par essence (ggplot2)
  # ================================================================
  tryCatch({
    detected <- statistics[statistics$n_pixels > 0, ]
    detected_sp <- detected[!detected$espece %in% c("Coupe/Vide", "Cleared"), ]
    if (nrow(detected_sp) > 0 && "confiance_moy" %in% names(detected_sp)) {
      p2 <- ggplot(detected_sp,
                    aes(x = reorder(espece, confiance_moy),
                        y = confiance_moy, fill = espece)) +
        geom_col(alpha = 0.85, show.legend = FALSE) +
        geom_text(aes(label = paste0(confiance_moy, "%")),
                  hjust = -0.1, size = 3.5) +
        geom_hline(yintercept = 70, linetype = "dashed",
                   color = "grey50", alpha = 0.7) +
        annotate("text", x = 0.8, y = 72, label = "Seuil 70%",
                 hjust = 0, size = 2.5, color = "grey50",
                 fontface = "italic") +
        coord_flip(clip = "off") +
        scale_fill_manual(values = species_colors) +
        scale_y_continuous(limits = c(0, 105),
                           expand = expansion(mult = c(0, 0.05))) +
        labs(
          title = "Confiance moyenne de classification par essence",
          subtitle = paste0("Probabilit\u00e9 moyenne attribu\u00e9e par le mod\u00e8le ",
                            "aux pixels de chaque essence"),
          x = NULL, y = "Confiance moyenne (%)"
        ) +
        theme_minimal(base_size = 11) +
        theme(
          plot.title = element_text(face = "bold", size = 14),
          plot.subtitle = element_text(size = 10, color = "grey40"),
          axis.text.y = element_text(face = "italic", size = 9),
          panel.grid.major.y = element_blank()
        )
      print(p2)
    }
  }, error = function(e) {
    plot.new()
    text(0.5, 0.5, paste("Erreur confiance :", e$message), cex = 0.8)
  })

  # Fermeture explicite (on.exit aussi en place comme filet de s\u00e9curit\u00e9)
  grDevices::dev.off()
  on.exit(NULL)  # annuler on.exit apr\u00e8s fermeture r\u00e9ussie

  log_msg("Rapport PDF : {pdf_path}", level = "success")
  invisible(pdf_path)
}

