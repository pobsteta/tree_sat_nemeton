# ==============================================================================
# TreeSatAI-Time-Series — Package hooks
# ==============================================================================

.onLoad <- function(libname, pkgname) {
  # Initialiser les chemins du projet
  root <- getOption("treesatnemeton.project_root", default = getwd())

  assign("PROJECT_ROOT", root, envir = parent.env(environment()))
  assign("DATA_DIR", file.path(root, "data"), envir = parent.env(environment()))
  assign("RAW_DIR", file.path(root, "data", "raw"), envir = parent.env(environment()))
  assign("PROCESSED_DIR", file.path(root, "data", "processed"), envir = parent.env(environment()))
  assign("TS_DIR", file.path(root, "data", "timeseries"), envir = parent.env(environment()))
  assign("OUTPUT_DIR", file.path(root, "output"), envir = parent.env(environment()))
  assign("FIGURES_DIR", file.path(root, "figures"), envir = parent.env(environment()))
  assign("MODELS_DIR", file.path(root, "output", "models"), envir = parent.env(environment()))
}

.onAttach <- function(libname, pkgname) {
  packageStartupMessage(
    cli::format_message(c(
      "{.strong treesatnemeton} {utils::packageVersion(pkgname)}",
      "",
      "Classification d'essences foresti\u00e8res par s\u00e9ries temporelles Sentinel-2",
      "{.emph 20 esp\u00e8ces europ\u00e9ennes \u2014 Signatures ph\u00e9nologiques annuelles}",
      "",
      "i" = "D\u00e9marrage rapide :",
      " " = "  {.code result <- predict_species_map(\"aoi.gpkg\", auto_download = TRUE)}",
      " " = "  {.code download_treesatai_hf()}  # T\u00e9l\u00e9charger le dataset r\u00e9el"
    ))
  )
}
