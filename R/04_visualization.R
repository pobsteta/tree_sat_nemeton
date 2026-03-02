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

#' Trac\u00e9 de l'importance de TOUTES les variables (Random Forest)
#'
#' Affiche l'ensemble des variables class\u00e9es par importance d\u00e9croissante,
#' color\u00e9es par type de feature. Utile pour visualiser la contribution
#' globale de chaque feature, y compris celles de faible importance.
#'
#' @param importance_df data.frame avec colonnes variable, importance
#' @param title Titre
#' @param save_path Chemin de sauvegarde
#' @return Objet ggplot
plot_variable_importance_all <- function(importance_df,
                                          title = "Importance de toutes les variables",
                                          save_path = NULL) {
  imp <- importance_df |>
    dplyr::arrange(dplyr::desc(importance)) |>
    dplyr::mutate(
      feature_type = dplyr::case_when(
        grepl("^pheno_|_pheno_", variable)   ~ "Ph\u00e9nologie",
        grepl("^NDVI|^EVI|^NDWI|^CRI|^NBR", variable) ~ "Indice spectral",
        grepl("^fourier_|_fourier_", variable) ~ "Fourier",
        grepl("^B\\d|^B8A", variable)          ~ "Bande spectrale",
        grepl("^S1_|^s1_", variable)           ~ "Radar S1",
        grepl("^terrain_|^elev|^slope|^aspect|^twi|^tpi", variable) ~ "Terrain",
        TRUE                                    ~ "Autre"
      ),
      rank = dplyr::row_number()
    )

  type_colors <- c(
    "Ph\u00e9nologie"       = "#e41a1c",
    "Indice spectral"  = "#377eb8",
    "Fourier"          = "#4daf4a",
    "Bande spectrale"  = "#984ea3",
    "Radar S1"         = "#ff7f00",
    "Terrain"          = "#a65628",
    "Autre"            = "#999999"
  )

  n_vars <- nrow(imp)
  # Adapter la taille des labels selon le nombre de variables
  label_size <- if (n_vars > 100) 3 else if (n_vars > 50) 4.5 else 6

  p <- ggplot(imp, aes(x = reorder(variable, importance),
                         y = importance, fill = feature_type)) +
    geom_col(alpha = 0.85) +
    coord_flip() +
    scale_fill_manual(values = type_colors, name = "Type de feature") +
    labs(
      title = title,
      subtitle = glue::glue("{n_vars} variables (Random Forest MDA)"),
      x = NULL,
      y = "Mean Decrease Accuracy"
    ) +
    theme_minimal(base_size = VIS_PARAMS$font_size) +
    theme(
      axis.text.y = element_text(size = label_size),
      plot.title = element_text(face = "bold"),
      legend.position = "bottom"
    )

  if (!is.null(save_path)) {
    # Hauteur adaptative selon le nombre de variables
    h <- max(20, n_vars * 0.35)
    ggsave(save_path, p, width = VIS_PARAMS$width_cm, height = h,
           units = "cm", dpi = VIS_PARAMS$dpi, limitsize = FALSE)
    log_msg("Importance (toutes variables) sauvegard\u00e9e : {save_path}", level = "success")
  }

  p
}

#' Trac\u00e9 de l'importance Boruta (d\u00e9cisions Confirmed / Rejected)
#'
#' Affiche l'importance moyenne de chaque variable telle qu'\u00e9valu\u00e9e par
#' Boruta, color\u00e9e par d\u00e9cision (Confirmed = retenue, Rejected = rejet\u00e9e,
#' Tentative = ind\u00e9cise). Permet de comprendre pourquoi certaines variables
#' ont \u00e9t\u00e9 \u00e9limin\u00e9es.
#'
#' @param boruta_importance_df data.frame avec colonnes variable, decision, meanImp
#'   (issu de select_features_boruta()$importance_df ou lu depuis boruta_importance.csv)
#' @param title Titre
#' @param save_path Chemin de sauvegarde
#' @return Objet ggplot
plot_boruta_importance <- function(boruta_importance_df,
                                    title = "S\u00e9lection Boruta \u2014 Importance et d\u00e9cisions",
                                    save_path = NULL) {
  imp <- boruta_importance_df |>
    dplyr::arrange(dplyr::desc(meanImp)) |>
    dplyr::mutate(
      decision_label = dplyr::case_when(
        decision == "Confirmed" ~ "Retenue",
        decision == "Rejected"  ~ "Rejet\u00e9e",
        decision == "Tentative" ~ "Tentative",
        TRUE ~ decision
      ),
      feature_type = dplyr::case_when(
        grepl("^pheno_|_pheno_", variable)   ~ "Ph\u00e9nologie",
        grepl("^NDVI|^EVI|^NDWI|^CRI|^NBR", variable) ~ "Indice spectral",
        grepl("^fourier_|_fourier_", variable) ~ "Fourier",
        grepl("^B\\d|^B8A", variable)          ~ "Bande spectrale",
        grepl("^S1_|^s1_", variable)           ~ "Radar S1",
        grepl("^terrain_|^elev|^slope|^aspect|^twi|^tpi", variable) ~ "Terrain",
        TRUE                                    ~ "Autre"
      )
    )

  decision_colors <- c(
    "Retenue"   = "#1a9850",
    "Rejet\u00e9e"   = "#d73027",
    "Tentative" = "#fee08b"
  )

  n_confirmed <- sum(imp$decision == "Confirmed")
  n_rejected  <- sum(imp$decision == "Rejected")
  n_tentative <- sum(imp$decision == "Tentative")
  n_total     <- nrow(imp)

  label_size <- if (n_total > 100) 3 else if (n_total > 50) 4.5 else 6

  p <- ggplot(imp, aes(x = reorder(variable, meanImp),
                         y = meanImp, fill = decision_label)) +
    geom_col(alpha = 0.85) +
    coord_flip() +
    scale_fill_manual(values = decision_colors, name = "D\u00e9cision Boruta") +
    labs(
      title = title,
      subtitle = glue::glue(
        "{n_total} variables \u00e9valu\u00e9es : ",
        "{n_confirmed} retenues, {n_rejected} rejet\u00e9es",
        if (n_tentative > 0) paste0(", ", n_tentative, " tentatives") else ""
      ),
      x = NULL,
      y = "Importance moyenne (Boruta)"
    ) +
    theme_minimal(base_size = VIS_PARAMS$font_size) +
    theme(
      axis.text.y = element_text(size = label_size),
      plot.title = element_text(face = "bold"),
      legend.position = "bottom"
    )

  if (!is.null(save_path)) {
    h <- max(20, n_total * 0.35)
    ggsave(save_path, p, width = VIS_PARAMS$width_cm, height = h,
           units = "cm", dpi = VIS_PARAMS$dpi, limitsize = FALSE)
    log_msg("Importance Boruta sauvegard\u00e9e : {save_path}", level = "success")
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
#' @param forest_mask SpatRaster binaire du masque forestier (NULL = pas de carte masque)
#' @return Chemin vers le fichier PDF (invisible)
#' @export
generate_prediction_report_pdf <- function(rasters, statistics, output_dir,
                                            s2_rgb = NULL, aoi = NULL,
                                            forest_mask = NULL) {
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
  # PAGE 5 : Masque forestier (OSO + NDVI)
  # ================================================================
  if (!is.null(forest_mask)) {
    tryCatch({
      par(mar = c(3.5, 2, 3.5, 8))
      mask_cols <- c("#f7f7f7", "#1a9850")
      terra::plot(forest_mask, type = "classes",
                  col = mask_cols,
                  main = "Masque forestier (OSO + NDVI)",
                  legend = FALSE, axes = TRUE)
      if (!is.null(aoi_proj)) {
        plot(sf::st_geometry(aoi_proj), add = TRUE,
             border = "black", lwd = 2, lty = 2)
      }
      par(xpd = TRUE)
      legend("right", inset = c(-0.08, 0),
             legend = c("Non-for\u00eat", "For\u00eat"),
             fill = mask_cols,
             cex = 0.8, bty = "n", title = "Masque", title.font = 2)
      par(xpd = FALSE)
      n_forest <- sum(terra::values(forest_mask) == 1L, na.rm = TRUE)
      n_total  <- sum(!is.na(terra::values(forest_mask)))
      pct_forest <- round(n_forest / n_total * 100, 1)
      res_m <- terra::res(forest_mask)[1]
      ha_forest <- round(n_forest * res_m^2 / 10000, 1)
      mtext(paste0("For\u00eat : ", n_forest, " pixels (",
                    pct_forest, "%, ", ha_forest, " ha) | ",
                    "Combinaison : ", FOREST_MASK_PARAMS$combine_method),
            side = 1, line = 2.2, cex = 0.7, font = 3)
    }, error = function(e) {
      plot.new()
      text(0.5, 0.5, paste("Erreur carte masque :", e$message), cex = 0.8)
    })
  }

  # ================================================================
  # PAGE 6 : Richesse sp\u00e9cifique (nombre d'essences par pixel)
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
  # PAGE 7 : Composition foresti\u00e8re (ggplot2)
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
  # PAGE 8 : Confiance moyenne par essence (ggplot2)
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

# --- Helpers internes pour les cartes ggplot2 ---------------------------------

#' Convertir un SpatRaster mono-bande en data.frame pour geom_raster
#' @param r SpatRaster mono-bande
#' @param max_cells Sous-\u00e9chantillonnage si le raster d\u00e9passe ce seuil
#' @return data.frame avec colonnes x, y, value
#' @keywords internal
.rast_to_df <- function(r, max_cells = 500000) {
  if (terra::ncell(r) > max_cells) {
    fact <- ceiling(sqrt(terra::ncell(r) / max_cells))
    r <- terra::aggregate(r, fact = fact, fun = "modal", na.rm = TRUE)
  }
  df <- terra::as.data.frame(r, xy = TRUE, na.rm = TRUE)
  names(df)[3] <- "value"
  df
}

#' Convertir un SpatRaster RGB (3 bandes) en data.frame pour geom_raster
#' @param rgb SpatRaster 3 bandes (R, G, B)
#' @param max_cells Sous-\u00e9chantillonnage si le raster d\u00e9passe ce seuil
#' @return data.frame avec colonnes x, y, hex (couleur hexad\u00e9cimale)
#' @keywords internal
.rgb_to_df <- function(rgb, max_cells = 500000) {
  if (terra::ncell(rgb) > max_cells) {
    fact <- ceiling(sqrt(terra::ncell(rgb) / max_cells))
    rgb <- terra::aggregate(rgb, fact = fact, fun = "mean", na.rm = TRUE)
  }
  df <- terra::as.data.frame(rgb, xy = TRUE, na.rm = TRUE)
  # Stretch lin\u00e9aire 2-98%
  stretch_band <- function(v) {
    q <- quantile(v, c(0.02, 0.98), na.rm = TRUE)
    v <- (v - q[1]) / (q[2] - q[1])
    pmin(pmax(v, 0), 1)
  }
  r_s <- stretch_band(df[[3]])
  g_s <- stretch_band(df[[4]])
  b_s <- stretch_band(df[[5]])
  df$hex <- grDevices::rgb(r_s, g_s, b_s)
  df[, c("x", "y", "hex")]
}

#' Th\u00e8me minimaliste pour les cartes ggplot2
#' @keywords internal
.theme_map <- function(base_size = 9) {
  theme_minimal(base_size = base_size) +
    theme(
      axis.title   = element_blank(),
      axis.text    = element_text(size = 6),
      plot.title   = element_text(face = "bold", size = base_size + 2),
      plot.subtitle = element_text(size = base_size - 1, color = "grey40"),
      legend.key.height = unit(0.8, "cm"),
      legend.key.width  = unit(0.3, "cm"),
      legend.title = element_text(size = base_size - 1, face = "bold"),
      legend.text  = element_text(size = base_size - 2),
      panel.grid   = element_line(color = "grey92", linewidth = 0.2)
    )
}

# --- Rapport cartographique RStudio (ggplot2 + patchwork) ---------------------

#' G\u00e9n\u00e9ration d'un rapport cartographique ggplot2/patchwork
#'
#' Produit un objet patchwork affich\u00e9 dans le plot pane RStudio ET
#' sauvegard\u00e9 en PDF. Toutes les cartes utilisent geom_raster (pas de
#' d\u00e9pendance \u00e0 terra::plot) : compatible RStudio, Quarto, Shiny.
#'
#' @inheritParams generate_prediction_report_pdf
#' @return Liste avec \code{dashboard} (objet patchwork), \code{pdf_path}
#'   (chemin PDF) et \code{plots} (liste individuelle)
#' @export
generate_prediction_report_rstudio <- function(rasters, statistics, output_dir,
                                                s2_rgb = NULL, aoi = NULL,
                                                forest_mask = NULL) {
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  log_msg("G\u00e9n\u00e9ration du rapport cartographique RStudio (ggplot2 + patchwork)...")

  # --- M\u00e9tadonn\u00e9es ---
  legend_df <- rasters$legend
  class_names <- legend_df$species
  n_classes <- length(class_names)
  species_colors <- .match_species_colors(class_names)
  cleared_names <- intersect(class_names, c("Coupe/Vide", "Cleared"))
  cleared_idx   <- which(class_names %in% cleared_names)

  # AOI en sf pour coord_sf
  aoi_geom <- NULL
  if (!is.null(aoi)) {
    aoi_geom <- sf::st_geometry(
      sf::st_transform(aoi, terra::crs(rasters$species))
    )
  }

  plots <- list()

  # ============================================================
  # 1. Carte des essences (fond satellite optionnel)
  # ============================================================
  tryCatch({
    # Raster essences sans Coupe/Vide
    r_sp <- rasters$species
    if (length(cleared_idx) > 0) {
      rcl_na <- cbind(cleared_idx, rep(NA_real_, length(cleared_idx)))
      r_sp <- terra::classify(r_sp, rcl_na)
    }
    if (!is.null(terra::levels(r_sp)[[1]])) levels(r_sp) <- NULL
    df_sp <- .rast_to_df(r_sp)
    df_sp$espece <- class_names[df_sp$value]

    active_colors <- species_colors[!names(species_colors) %in% cleared_names]

    p1 <- ggplot()
    # Fond satellite via annotation_raster (pas de d\u00e9pendance ggnewscale)
    if (!is.null(s2_rgb)) {
      df_rgb <- .rgb_to_df(s2_rgb)
      ext_rgb <- terra::ext(s2_rgb)
      # Construire une matrice de couleurs pour annotation_raster
      rgb_agg <- s2_rgb
      if (terra::ncell(rgb_agg) > 500000) {
        fact <- ceiling(sqrt(terra::ncell(rgb_agg) / 500000))
        rgb_agg <- terra::aggregate(rgb_agg, fact = fact, fun = "mean",
                                     na.rm = TRUE)
      }
      rgb_mat <- terra::as.matrix(rgb_agg, wide = TRUE)
      stretch_v <- function(v) {
        q <- quantile(v, c(0.02, 0.98), na.rm = TRUE)
        pmin(pmax((v - q[1]) / (q[2] - q[1]), 0), 1)
      }
      nr <- terra::nrow(rgb_agg)
      nc <- terra::ncol(rgb_agg)
      r_v <- stretch_v(rgb_mat[, seq_len(nc)])
      g_v <- stretch_v(rgb_mat[, nc + seq_len(nc)])
      b_v <- stretch_v(rgb_mat[, 2 * nc + seq_len(nc)])
      hex_mat <- matrix(grDevices::rgb(r_v, g_v, b_v), nrow = nr, ncol = nc)
      # annotation_raster attend une matrice [nrow, ncol] orient\u00e9e top-to-bottom
      p1 <- p1 +
        annotation_raster(
          hex_mat,
          xmin = ext_rgb[1], xmax = ext_rgb[2],
          ymin = ext_rgb[3], ymax = ext_rgb[4]
        )
    }
    p1 <- p1 +
      geom_raster(data = df_sp,
                   aes(x = x, y = y, fill = espece),
                   alpha = if (!is.null(s2_rgb)) 0.7 else 1)
    if (!is.null(aoi_geom)) {
      p1 <- p1 + geom_sf(data = aoi_geom, fill = NA,
                           color = "white", linewidth = 0.6, linetype = 2,
                           inherit.aes = FALSE)
    }
    p1 <- p1 +
      scale_fill_manual(values = active_colors, name = "Essence",
                         na.translate = FALSE) +
      coord_sf(expand = FALSE) +
      labs(title = "Carte des essences foresti\u00e8res",
           subtitle = "Coupes/Vides = transparents") +
      .theme_map() +
      guides(fill = guide_legend(ncol = 1, override.aes = list(alpha = 1)))
    plots$species <- p1
  }, error = function(e) {
    plots$species <<- ggplot() +
      annotate("text", x = 0.5, y = 0.5,
               label = paste("Erreur carte essences :", e$message)) +
      theme_void()
  })

  # ============================================================
  # 2. Shannon
  # ============================================================
  if (!is.null(rasters$shannon)) {
    tryCatch({
      df_sh <- .rast_to_df(rasters$shannon)
      p2 <- ggplot(df_sh, aes(x = x, y = y, fill = value)) +
        geom_raster()
      if (!is.null(aoi_geom)) {
        p2 <- p2 + geom_sf(data = aoi_geom, fill = NA,
                             color = "white", linewidth = 0.6, linetype = 2,
                             inherit.aes = FALSE)
      }
      p2 <- p2 +
        scale_fill_viridis_c(option = "viridis", name = "Shannon",
                              limits = c(0, 1)) +
        coord_sf(expand = FALSE) +
        labs(title = "Entropie de Shannon",
             subtitle = "0 = pur | 1 = m\u00e9lange maximal") +
        .theme_map()
      plots$shannon <- p2
    }, error = function(e) {
      plots$shannon <<- ggplot() +
        annotate("text", x = 0.5, y = 0.5,
                 label = paste("Erreur Shannon :", e$message)) +
        theme_void()
    })
  }

  # ============================================================
  # 3. Confiance
  # ============================================================
  tryCatch({
    df_cf <- .rast_to_df(rasters$confidence)
    p3 <- ggplot(df_cf, aes(x = x, y = y, fill = value)) +
      geom_raster()
    if (!is.null(aoi_geom)) {
      p3 <- p3 + geom_sf(data = aoi_geom, fill = NA,
                           color = "black", linewidth = 0.6, linetype = 2,
                           inherit.aes = FALSE)
    }
    p3 <- p3 +
      scale_fill_gradientn(
        colors = c("#d73027", "#fc8d59", "#fee08b",
                   "#d9ef8b", "#91cf60", "#1a9850"),
        limits = c(0, 1), name = "Probabilit\u00e9"
      ) +
      coord_sf(expand = FALSE) +
      labs(title = "Carte de confiance",
           subtitle = "Rouge = faible | Vert = forte") +
      .theme_map()
    plots$confidence <- p3
  }, error = function(e) {
    plots$confidence <<- ggplot() +
      annotate("text", x = 0.5, y = 0.5,
               label = paste("Erreur confiance :", e$message)) +
      theme_void()
  })

  # ============================================================
  # 4. Feuillus / R\u00e9sineux
  # ============================================================
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

    df_type <- .rast_to_df(r_type)
    df_type$type_label <- c("Feuillus", "R\u00e9sineux")[df_type$value]

    p4 <- ggplot(df_type, aes(x = x, y = y, fill = type_label)) +
      geom_raster()
    if (!is.null(aoi_geom)) {
      p4 <- p4 + geom_sf(data = aoi_geom, fill = NA,
                           color = "black", linewidth = 0.6, linetype = 2,
                           inherit.aes = FALSE)
    }
    p4 <- p4 +
      scale_fill_manual(values = c("Feuillus" = "#66c2a5",
                                    "R\u00e9sineux" = "#1b7837"),
                         name = "Type") +
      coord_sf(expand = FALSE) +
      labs(title = "Types forestiers",
           subtitle = "Feuillus / R\u00e9sineux") +
      .theme_map()
    plots$type <- p4
  }, error = function(e) {
    plots$type <<- ggplot() +
      annotate("text", x = 0.5, y = 0.5,
               label = paste("Erreur types :", e$message)) +
      theme_void()
  })

  # ============================================================
  # 5. Masque forestier
  # ============================================================
  if (!is.null(forest_mask)) {
    tryCatch({
      df_mask <- .rast_to_df(forest_mask)
      df_mask$label <- ifelse(df_mask$value == 1, "For\u00eat", "Non-for\u00eat")

      p5 <- ggplot(df_mask, aes(x = x, y = y, fill = label)) +
        geom_raster()
      if (!is.null(aoi_geom)) {
        p5 <- p5 + geom_sf(data = aoi_geom, fill = NA,
                             color = "black", linewidth = 0.6, linetype = 2,
                             inherit.aes = FALSE)
      }
      n_forest <- sum(df_mask$value == 1)
      n_total  <- nrow(df_mask)
      pct_f <- round(n_forest / n_total * 100, 1)

      p5 <- p5 +
        scale_fill_manual(values = c("Non-for\u00eat" = "#f7f7f7",
                                      "For\u00eat" = "#1a9850"),
                           name = "Masque") +
        coord_sf(expand = FALSE) +
        labs(title = "Masque forestier (OSO + NDVI)",
             subtitle = paste0("For\u00eat : ", pct_f, "% de la zone")) +
        .theme_map()
      plots$mask <- p5
    }, error = function(e) {
      plots$mask <<- ggplot() +
        annotate("text", x = 0.5, y = 0.5,
                 label = paste("Erreur masque :", e$message)) +
        theme_void()
    })
  }

  # ============================================================
  # 6. Richesse sp\u00e9cifique
  # ============================================================
  if (!is.null(rasters$presence)) {
    tryCatch({
      r_rich <- terra::app(rasters$presence, sum, na.rm = TRUE)
      df_rich <- .rast_to_df(r_rich)
      p6 <- ggplot(df_rich, aes(x = x, y = y, fill = value)) +
        geom_raster()
      if (!is.null(aoi_geom)) {
        p6 <- p6 + geom_sf(data = aoi_geom, fill = NA,
                             color = "white", linewidth = 0.6, linetype = 2,
                             inherit.aes = FALSE)
      }
      p6 <- p6 +
        scale_fill_viridis_c(option = "magma", name = "Nb essences") +
        coord_sf(expand = FALSE) +
        labs(title = "Richesse sp\u00e9cifique",
             subtitle = "Nombre d'essences d\u00e9tect\u00e9es par pixel") +
        .theme_map()
      plots$richness <- p6
    }, error = function(e) {
      plots$richness <<- ggplot() +
        annotate("text", x = 0.5, y = 0.5,
                 label = paste("Erreur richesse :", e$message)) +
        theme_void()
    })
  }

  # ============================================================
  # 7. Composition foresti\u00e8re (barplot)
  # ============================================================
  tryCatch({
    detected <- statistics[statistics$n_pixels > 0, ]
    detected_sp <- detected[!detected$espece %in% c("Coupe/Vide", "Cleared"), ]
    if (nrow(detected_sp) > 0) {
      p7 <- ggplot(detected_sp,
                     aes(x = reorder(espece, surface_ha),
                         y = surface_ha, fill = espece)) +
        geom_col(alpha = 0.85, show.legend = FALSE) +
        geom_text(aes(label = paste0(pct, "%")),
                  hjust = -0.1, size = 2.5) +
        coord_flip(clip = "off") +
        scale_fill_manual(values = species_colors) +
        scale_y_continuous(expand = expansion(mult = c(0, 0.3))) +
        labs(title = "Composition foresti\u00e8re",
             subtitle = paste0(sum(detected_sp$surface_ha), " ha"),
             x = NULL, y = "Surface (ha)") +
        theme_minimal(base_size = 9) +
        theme(
          plot.title = element_text(face = "bold", size = 11),
          plot.subtitle = element_text(size = 8, color = "grey40"),
          axis.text.y = element_text(face = "italic", size = 7),
          panel.grid.major.y = element_blank()
        )
      plots$composition <- p7
    }
  }, error = function(e) {
    plots$composition <<- ggplot() +
      annotate("text", x = 0.5, y = 0.5,
               label = paste("Erreur composition :", e$message)) +
      theme_void()
  })

  # ============================================================
  # 8. Confiance par esp\u00e8ce (barplot)
  # ============================================================
  tryCatch({
    detected <- statistics[statistics$n_pixels > 0, ]
    detected_sp <- detected[!detected$espece %in% c("Coupe/Vide", "Cleared"), ]
    if (nrow(detected_sp) > 0 && "confiance_moy" %in% names(detected_sp)) {
      p8 <- ggplot(detected_sp,
                     aes(x = reorder(espece, confiance_moy),
                         y = confiance_moy, fill = espece)) +
        geom_col(alpha = 0.85, show.legend = FALSE) +
        geom_text(aes(label = paste0(confiance_moy, "%")),
                  hjust = -0.1, size = 2.5) +
        geom_hline(yintercept = 70, linetype = "dashed",
                   color = "grey50", alpha = 0.7) +
        coord_flip(clip = "off") +
        scale_fill_manual(values = species_colors) +
        scale_y_continuous(limits = c(0, 105),
                           expand = expansion(mult = c(0, 0.05))) +
        labs(title = "Confiance par essence",
             subtitle = "Probabilit\u00e9 moyenne (%)",
             x = NULL, y = "Confiance (%)") +
        theme_minimal(base_size = 9) +
        theme(
          plot.title = element_text(face = "bold", size = 11),
          plot.subtitle = element_text(size = 8, color = "grey40"),
          axis.text.y = element_text(face = "italic", size = 7),
          panel.grid.major.y = element_blank()
        )
      plots$conf_species <- p8
    }
  }, error = function(e) {
    plots$conf_species <<- ggplot() +
      annotate("text", x = 0.5, y = 0.5,
               label = paste("Erreur confiance :", e$message)) +
      theme_void()
  })

  # ============================================================
  # Assemblage patchwork
  # ============================================================
  # Ligne 1 : essences + shannon + confiance
  # Ligne 2 : types + masque + richesse
  # Ligne 3 : composition + confiance/esp\u00e8ce
  available <- names(plots)

  # Construire les lignes adaptativement
  row1 <- list()
  if ("species" %in% available) row1 <- c(row1, list(plots$species))
  if ("shannon" %in% available) row1 <- c(row1, list(plots$shannon))
  if ("confidence" %in% available) row1 <- c(row1, list(plots$confidence))

  row2 <- list()
  if ("type" %in% available) row2 <- c(row2, list(plots$type))
  if ("mask" %in% available) row2 <- c(row2, list(plots$mask))
  if ("richness" %in% available) row2 <- c(row2, list(plots$richness))

  row3 <- list()
  if ("composition" %in% available) row3 <- c(row3, list(plots$composition))
  if ("conf_species" %in% available) row3 <- c(row3, list(plots$conf_species))

  # Assembler chaque ligne avec patchwork::wrap_plots
  build_row <- function(plot_list) {
    if (length(plot_list) == 0) return(NULL)
    patchwork::wrap_plots(plot_list, nrow = 1)
  }

  rows <- Filter(Negate(is.null), list(
    build_row(row1), build_row(row2), build_row(row3)
  ))

  dashboard <- patchwork::wrap_plots(rows, ncol = 1) +
    patchwork::plot_annotation(
      title = "TreeSatAI Nemeton \u2014 Rapport cartographique",
      subtitle = glue::glue(
        "{n_classes} classes | {sum(statistics$surface_ha)} ha | ",
        "Sentinel-2 s\u00e9ries temporelles"
      ),
      theme = theme(
        plot.title = element_text(size = 16, face = "bold"),
        plot.subtitle = element_text(size = 11, color = "grey40")
      )
    )

  # --- Afficher dans le plot pane RStudio ---
  print(dashboard)
  log_msg("Rapport cartographique affich\u00e9 dans le viewer RStudio", level = "success")

  # --- Sauvegarde PDF ---
  pdf_path <- file.path(output_dir, "rapport_cartographique_rstudio.pdf")
  ggsave(pdf_path, dashboard,
         width = 42, height = 55, units = "cm",
         dpi = VIS_PARAMS$dpi, limitsize = FALSE)
  log_msg("Rapport PDF sauvegard\u00e9 : {pdf_path}", level = "success")

  invisible(list(dashboard = dashboard, pdf_path = pdf_path, plots = plots))
}

