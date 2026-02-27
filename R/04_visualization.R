#!/usr/bin/env Rscript
# ==============================================================================
# TreeSatAI-Time-Series — Visualisation
# Profils phénologiques, cartes de confusion, importance des variables
# ==============================================================================

# Configuration chargée via le package
# ggplot2 et patchwork sont des dépendances du package (Imports)

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

  p <- ggplot(profiles, aes(x = doy, y = mean_val, color = species_name)) +
    geom_ribbon(aes(ymin = mean_val - sd_val, ymax = mean_val + sd_val,
                    fill = species_name), alpha = 0.1, color = NA) +
    geom_line(linewidth = 0.8) +
    scale_color_manual(values = SPECIES_COLORS, name = "Espèce") +
    scale_fill_manual(values = SPECIES_COLORS, guide = "none") +
    scale_x_continuous(
      breaks = c(1, 32, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335),
      labels = c("Jan", "Fév", "Mar", "Avr", "Mai", "Jun",
                 "Jul", "Aoû", "Sep", "Oct", "Nov", "Déc")
    ) +
    labs(
      title = paste("Profils phénologiques —", index_name),
      subtitle = "Moyenne ± écart-type par espèce (TreeSatAI-TS, 20 espèces)",
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
    log_msg("Profils phénologiques sauvegardés : {save_path}", level = "success")
  }

  p
}

#' Profils par type : feuillus vs résineux
#' @param ts_long Dataset long
#' @param save_path Chemin de sauvegarde
#' @return Objet ggplot
plot_profiles_by_type <- function(ts_long, save_path = NULL) {
  # Ajouter le type d'essence
  ts_typed <- ts_long |>
    dplyr::left_join(
      SPECIES |> dplyr::select(french, type, phenologie),
      by = c("species_name" = "french")
    )

  profiles <- ts_typed |>
    dplyr::group_by(species_name, type, phenologie, date) |>
    dplyr::summarise(
      mean_ndvi = mean(NDVI, na.rm = TRUE),
      .groups = "drop"
    ) |>
    dplyr::mutate(doy = as.numeric(format(date, "%j")))

  p <- ggplot(profiles, aes(x = doy, y = mean_ndvi, color = species_name)) +
    geom_line(linewidth = 0.7, alpha = 0.8) +
    facet_grid(type ~ phenologie, scales = "free_y") +
    scale_color_manual(values = SPECIES_COLORS, name = "Espèce") +
    scale_x_continuous(
      breaks = c(1, 91, 182, 274),
      labels = c("Jan", "Avr", "Jul", "Oct")
    ) +
    labs(
      title = "Profils phénologiques par type et régime foliaire",
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
    log_msg("Profils par type sauvegardés : {save_path}", level = "success")
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
    class_names <- rownames(conf_matrix)
  }

  # Normaliser par ligne (recall par classe)
  conf_norm <- sweep(conf_matrix, 1, rowSums(conf_matrix), "/")
  conf_norm[is.nan(conf_norm)] <- 0

  # Conversion en data.frame long
  conf_df <- expand.grid(
    Reference  = class_names,
    Prediction = class_names
  )
  conf_df$Count    <- as.vector(conf_matrix)
  conf_df$Percent  <- as.vector(conf_norm) * 100

  # Labels
  conf_df$label <- ifelse(
    conf_df$Count > 0,
    paste0(conf_df$Count, "\n(", round(conf_df$Percent, 1), "%)"),
    ""
  )

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
      subtitle = paste0("OA = ", round(sum(diag(conf_matrix)) / sum(conf_matrix) * 100, 1), "%"),
      x = "Prédiction",
      y = "Référence"
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
    dplyr::summarise(dplyr::across(dplyr::all_of(ndvi_cols), mean, na.rm = TRUE),
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
    main = "Signatures phénologiques NDVI — 20 espèces (z-score)",
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
    ) |>
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
      subtitle = "20 espèces européennes — Sentinel-2 séries temporelles annuelles",
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

