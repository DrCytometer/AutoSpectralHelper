# app.R — AutoSpectral Control File Builder
# FCS I/O uses AutoSpectral's native readFCSheader / readFCS / writeFCS.
# Table uses DT; constrained columns edited via row-editor panel.
# Greying of mixed fields uses shinyjs.

# ---- Dependencies ----
required.packages <- c("shiny", "shinyFiles", "DT", "dplyr", "shinybusy", "shinyjs")
missing.packages  <- setdiff(required.packages, rownames(installed.packages()))
if (length(missing.packages) > 0)
  invisible(install.packages(missing.packages))
if (!"AutoSpectral" %in% rownames(installed.packages())) {
  if (!"devtools" %in% rownames(installed.packages()))
    invisible(install.packages("devtools"))
  devtools::install_github("DrCytometer/AutoSpectral")
}
invisible(lapply(c(required.packages, "AutoSpectral"), require, character.only = TRUE))

.launch_wd <- getOption("autospectral.launch_wd", default = getwd())

# ---- Cytometer map ----
map <- list(
  "Aurora"          = "aurora",
  "Northern Lights" = "auroraNL",
  "ID7000"          = "id7000",
  "FACSDiscover A8" = "a8",
  "FACSDiscover S8" = "s8",
  "Opteon"          = "opteon",
  "Mosaic"          = "mosaic",
  "Xenith"          = "xenith",
  "Symphony A5SE"   = "a5se"
)

# ---- Filesystem roots ----
.build_roots <- function() {
  base <- c("Working dir" = .launch_wd, "Home" = normalizePath("~"))
  base <- base[!duplicated(base)]
  if (.Platform$OS.type == "windows") {
    ok <- vapply(LETTERS, function(L) {
      tryCatch(isTRUE(file.info(paste0(L, ":/"))$isdir),
               warning = function(w) FALSE, error = function(e) FALSE)
    }, logical(1))
    ready <- LETTERS[ok]
    if (length(ready) > 0) {
      drive_paths <- setNames(paste0(ready, ":/"), paste0(ready, ":"))
      base <- c(base, drive_paths[!drive_paths %in% base])
    }
  } else {
    base <- c(base, c("/" = "/"))
  }
  base
}
.roots <- .build_roots()

# ---- Helper: build control data.frame ----
create.control.df <- function(control.dir, asp) {
  data.path <- system.file("extdata", "fluorophore_database.csv", package = "AutoSpectral")
  fluorophore.database <- read.csv(data.path, stringsAsFactors = FALSE)
  fluorophore.database[fluorophore.database == ""] <- NA

  control.files <- list.files(control.dir, pattern = "\\.fcs$", ignore.case = TRUE)
  if (length(control.files) < 1) stop("No FCS files found in the selected directory.")

  control.colnames <- c("filename", "fluorophore", "marker", "channel",
                        "control.type", "universal.negative",
                        "large.gate", "gate.name", "gate.define")

  df <- data.frame(matrix(ncol = length(control.colnames), nrow = length(control.files)),
                   stringsAsFactors = FALSE)
  colnames(df) <- control.colnames
  df$filename  <- control.files

  fluor.matches <- suppressMessages(
    AutoSpectral::match.fluorophores(control.files, fluorophore.database))
  df$fluorophore <- fluor.matches[df$filename]

  detectors    <- .get.detectors(asp, fluorophore.database)
  df$channel   <- detectors[match(df$fluorophore, names(detectors))]

  df$control.type <- sapply(df$filename, function(f) {
    if (grepl("cells", f, ignore.case = TRUE)) "cells"
    else if (grepl("beads", f, ignore.case = TRUE)) "beads"
    else ""
  })

  marker.data.path <- system.file("extdata", "marker_database.csv", package = "AutoSpectral")
  marker.database  <- read.csv(marker.data.path, stringsAsFactors = FALSE)
  marker.database[marker.database == ""] <- NA
  marker.matches <- suppressMessages(
    AutoSpectral::match.markers(control.files, marker.database))
  df$marker <- marker.matches[df$filename]

  merged <- merge(df, fluorophore.database, by = "fluorophore", all.x = TRUE)
  laser.order <- c("DeepUV", "UV", "Violet", "Blue", "YellowGreen", "Red", "IR")
  merged$excitation.laser <- factor(merged$excitation.laser, levels = laser.order)
  merged <- merged[order(merged$excitation.laser, merged$nominal.wavelength), ]
  df <- merged[, c(control.colnames, "is.viability")]

  idx_uns <- grepl("Unstained", df$filename, ignore.case = TRUE)
  df$fluorophore[idx_uns] <- ifelse(df$control.type[idx_uns] == "cells", "AF", "Negative")
  df$fluorophore[grepl("Negative", df$filename, ignore.case = TRUE)] <- "Negative"

  df$gate.name   <- ""
  df$gate.define <- TRUE

  text.cols <- setdiff(colnames(df), c("gate.define", "is.viability", "large.gate"))
  for (cn in text.cols) {
    df[[cn]][is.na(df[[cn]])] <- ""
    df[[cn]][df[[cn]] == "NA"] <- ""
  }
  df$is.viability <- suppressWarnings(as.logical(df$is.viability))
  df$large.gate   <- suppressWarnings(as.logical(df$large.gate))
  df$gate.define  <- suppressWarnings(as.logical(df$gate.define))
  df$is_unstained <- FALSE

  database.path <- system.file("extdata", "cytometer_database.csv", package = "AutoSpectral")
  cytometers    <- read.csv(database.path, stringsAsFactors = FALSE)
  all.detectors <- .get.cytometer.channels(asp, cytometers)

  return(list(df = df, detectors = all.detectors))
}

.get.detectors <- function(asp, fluorophore.database) {
  col <- switch(asp$cytometer,
                "Aurora"          = if (identical(asp$cytometer.version, "NL")) "channel.NL" else "channel.Aurora",
                "ID7000"          = "channel.ID7000",
                "FACSDiscover A8" = "channel.s8",
                "FACSDiscover S8" = "channel.s8",
                "Opteon"          = "channel.opteon",
                "Mosaic"          = "channel.mosaic",
                "Xenith"          = "channel.xenith",
                "Symphony"        = "channel.A5SE",
                stop("Unsupported cytometer"))
  stats::setNames(fluorophore.database[[col]], fluorophore.database$fluorophore)
}

.get.cytometer.channels <- function(asp, cytometers) {
  col <- switch(asp$cytometer,
                "Aurora"          = if (identical(asp$cytometer.version, "NL")) "NorthernLights" else "Aurora",
                "ID7000"          = "ID7000",
                "FACSDiscover A8" = ,
                "FACSDiscover S8" = "Discover",
                "Opteon"          = "Opteon",
                "Mosaic"          = "Mosaic",
                "Xenith"          = "Xenith",
                "Symphony"        = "A5SE",
                stop("Unsupported cytometer"))
  ch <- cytometers[[col]]
  ch[!is.na(ch) & ch != ""]
}

normalize_check_result <- function(x) {
  empty <- data.frame(severity = character(0), rule = character(0),
                      message  = character(0), stringsAsFactors = FALSE)
  if (is.null(x) || (is.logical(x) && isTRUE(x))) return(empty)
  if (is.data.frame(x)) {
    if (nrow(x) == 0) return(empty)
    required <- c("severity", "rule", "message")
    for (m in setdiff(required, colnames(x))) x[[m]] <- NA_character_
    return(x[, required, drop = FALSE])
  }
  data.frame(severity = "error", rule = "legacy_return",
             message  = "check.control.file() returned an unexpected format",
             stringsAsFactors = FALSE)
}

prepare.for.export <- function(df) {
  df <- df[, setdiff(colnames(df), "is_unstained"), drop = FALSE]
  for (cn in c("large.gate", "is.viability", "gate.define"))
    if (cn %in% colnames(df))
      df[[cn]] <- ifelse(is.na(df[[cn]]), NA, as.character(df[[cn]]))
  df
}

.empty.row <- function(df) {
  r <- df[NA_integer_, , drop = FALSE]
  for (cn in colnames(r)) {
    if      (is.logical(r[[cn]])) r[[cn]] <- FALSE
    else if (is.numeric(r[[cn]])) r[[cn]] <- NA_real_
    else                          r[[cn]] <- ""
  }
  r
}

# Field spec table — maps UI div id, data column, overwrite checkbox id, and label
.field_specs <- list(
  list(div = "fld_control_type",  col = "control.type",       ow = "ow_control_type",  label = "control.type"),
  list(div = "fld_channel",       col = "channel",            ow = "ow_channel",       label = "channel"),
  list(div = "fld_universal_neg", col = "universal.negative", ow = "ow_universal_neg", label = "universal.negative"),
  list(div = "fld_gate_name",     col = "gate.name",          ow = "ow_gate_name",     label = "gate.name"),
  list(div = "fld_is_unstained",  col = "is_unstained",       ow = "ow_is_unstained",  label = "Is unstained control"),
  list(div = "fld_large_gate",    col = "large.gate",         ow = "ow_large_gate",    label = "large.gate"),
  list(div = "fld_is_viability",  col = "is.viability",       ow = "ow_is_viability",  label = "is.viability"),
  list(div = "fld_gate_define",   col = "gate.define",        ow = "ow_gate_define",   label = "gate.define")
)

# ============================================================
# ---- UI ----
# ============================================================
ui <- fluidPage(
  tags$head(
    tags$link(rel = "stylesheet",
              href = "https://fonts.googleapis.com/css2?family=DM+Mono:ital,wght@0,300;0,400;0,500;1,400&family=DM+Sans:opsz,wght@9..40,300;9..40,400;9..40,500;9..40,600&display=swap"),
    tags$style(HTML("
      *, *::before, *::after { box-sizing: border-box; }
      body { font-family: 'DM Sans', sans-serif; font-size: 13.5px;
             background: #f4f5f7; color: #1a1d23; margin: 0; }

      .as-banner { background: #0f1923; color: #e8edf2; padding: 18px 32px 14px;
        display: flex; align-items: baseline; gap: 16px; border-bottom: 3px solid #00b4d8; }
      .as-banner h1 { font-family: 'DM Mono', monospace; font-size: 19px;
        font-weight: 500; letter-spacing: 0.03em; margin: 0; color: #fff; }
      .as-banner .as-subtitle { font-size: 11.5px; color: #7a92a8;
        font-family: 'DM Mono', monospace; letter-spacing: 0.06em; text-transform: uppercase; }

      .as-outer { display: flex; min-height: calc(100vh - 57px); }

      .as-sidebar { width: 270px; min-width: 270px; background: #0f1923;
        color: #c8d8e4; padding: 24px 20px; display: flex; flex-direction: column;
        gap: 0; border-right: 1px solid #1e2d3d; }
      .as-sidebar .as-section-label { font-family: 'DM Mono', monospace; font-size: 9.5px;
        letter-spacing: 0.14em; text-transform: uppercase; color: #4a6278; margin: 22px 0 8px; }
      .as-sidebar .as-section-label:first-child { margin-top: 0; }
      .as-sidebar label { color: #8aa8bf !important; font-size: 12px !important;
        font-weight: 500 !important; margin-bottom: 4px !important; }
      .as-sidebar .form-control, .as-sidebar .selectize-input {
        background: #162130 !important; border: 1px solid #253d52 !important;
        color: #d4e6f1 !important; border-radius: 6px !important;
        font-size: 12.5px !important; font-family: 'DM Sans', sans-serif !important; }
      .as-sidebar .selectize-dropdown { background: #162130 !important;
        border: 1px solid #253d52 !important; color: #d4e6f1 !important; }
      .as-sidebar .selectize-dropdown .option { color: #d4e6f1 !important; }
      .as-sidebar .selectize-dropdown .option.active { background: #00b4d8 !important; color: #fff !important; }

      .dir-display { font-family: 'DM Mono', monospace; font-size: 10.5px; color: #4a8fa8;
        background: #0a1520; border: 1px solid #1e3348; border-radius: 5px; padding: 6px 10px;
        margin-top: 4px; word-break: break-all; min-height: 30px; line-height: 1.5; }

      .btn-as-primary { background: #00b4d8 !important; color: #0a1520 !important;
        border: none !important; border-radius: 6px !important; font-weight: 600 !important;
        font-size: 12.5px !important; padding: 7px 14px !important; width: 100%;
        margin-top: 6px; transition: background 0.15s; }
      .btn-as-primary:hover { background: #48cae4 !important; }
      .btn-as-save { background: #06d6a0 !important; color: #0a1520 !important;
        border: none !important; border-radius: 6px !important; font-weight: 700 !important;
        font-size: 13px !important; padding: 9px 14px !important; width: 100%;
        margin-top: 10px; transition: background 0.15s; }
      .btn-as-save:hover { background: #38e8b8 !important; }

      a.btn-dir { display: block; background: #162130 !important; color: #7ab8d0 !important;
        border: 1px dashed #253d52 !important; border-radius: 6px !important;
        font-size: 12px !important; padding: 7px 12px !important; width: 100%; cursor: pointer;
        text-align: left; text-decoration: none !important;
        transition: border-color 0.15s, color 0.15s; }
      a.btn-dir:hover { border-color: #00b4d8 !important; color: #00b4d8 !important; }

      .as-instructions { list-style: none; padding: 0; margin: 0; }
      .as-instructions li { padding: 5px 0 5px 20px; position: relative;
        font-size: 12px; color: #8aa8bf; line-height: 1.5; }
      .as-instructions li::before { content: '>'; position: absolute; left: 4px;
        color: #00b4d8; font-weight: 700; }
      hr.as-hr { border: none; border-top: 1px solid #1e2d3d; margin: 16px 0; }

      .as-main { flex: 1; padding: 28px 32px; overflow-x: auto; }
      .as-card { background: #fff; border-radius: 10px; border: 1px solid #dde3ea;
        padding: 20px 24px; margin-bottom: 22px; box-shadow: 0 1px 4px rgba(0,0,0,0.05); }
      .as-card-header { font-family: 'DM Mono', monospace; font-size: 10.5px;
        letter-spacing: 0.12em; text-transform: uppercase; color: #7a92a8;
        margin-bottom: 14px; padding-bottom: 10px; border-bottom: 1px solid #edf0f4;
        display: flex; align-items: center; gap: 8px; }
      .as-card-header .dot { width: 6px; height: 6px; border-radius: 50%;
        background: #00b4d8; display: inline-block; }

      .tbl-toolbar { display: flex; gap: 8px; flex-wrap: wrap; margin-bottom: 12px; align-items: center; }
      .btn-tool { background: #f4f5f7 !important; color: #3a5068 !important;
        border: 1px solid #dde3ea !important; border-radius: 6px !important;
        font-size: 12px !important; padding: 5px 12px !important; cursor: pointer;
        transition: background 0.12s, border-color 0.12s; white-space: nowrap; }
      .btn-tool:hover { background: #e8f4f8 !important; border-color: #00b4d8 !important;
        color: #006d8a !important; }
      .btn-tool-danger { background: #fff5f5 !important; color: #c0392b !important;
        border: 1px solid #f5c6c2 !important; border-radius: 6px !important;
        font-size: 12px !important; padding: 5px 12px !important; cursor: pointer; white-space: nowrap; }
      .btn-tool-danger:hover { background: #fde8e6 !important; }

      /* ── Row editor ── */
      .row-editor { background: #f0f8fb; border: 1px solid #c8e6f0; border-radius: 8px;
        padding: 16px 18px 180px; margin-top: 14px; }

      /* Grid of editor fields: fixed columns, no wrapping chaos */
      .editor-grid {
        display: grid;
        grid-template-columns: 140px 170px 230px 170px 160px 110px 110px 110px;
        gap: 10px 14px;
        align-items: end;
        overflow-x: auto;
      }

      .editor-field { display: flex; flex-direction: column; }
      .editor-field > label {
        font-family: 'DM Mono', monospace;
        font-size: 10.5px !important;
        font-weight: 600 !important;
        color: #4a7a8a !important;
        letter-spacing: 0.04em;
        text-transform: uppercase;
        margin-bottom: 4px !important;
        white-space: nowrap;
        transition: color 0.2s;
      }
      .editor-field .form-control, .editor-field .selectize-input {
        font-size: 12px !important; border-radius: 5px !important;
        border: 1px solid #b0d4e0 !important; background: #fff !important;
        transition: opacity 0.2s, border-color 0.2s; }
      .editor-field .checkbox label { font-size: 12px !important; color: #2d4356 !important; }

      /* Greyed-out mixed field */
      .field-mixed > label { color: #b0c8d4 !important; }
      .field-mixed .form-control,
      .field-mixed .selectize-input { opacity: 0.4 !important; border-style: dashed !important; }
      .field-mixed .checkbox label { opacity: 0.4; }
      .field-mixed .selectize-control { opacity: 0.4; }

      /* Overwrite row: single line, no wrapping */
      .overwrite-row {
        display: flex;
        flex-direction: row;
        flex-wrap: nowrap;
        gap: 0;
        margin-top: 12px;
        padding-top: 10px;
        border-top: 1px dashed #c8e6f0;
        overflow-x: auto;
        align-items: center;
      }
      .overwrite-label-hdr {
        font-family: 'DM Mono', monospace;
        font-size: 10px;
        color: #7a92a8;
        white-space: nowrap;
        margin-right: 10px;
        flex-shrink: 0;
      }
      /* Each overwrite item aligns under its editor-field column */
      .overwrite-item {
        display: flex;
        flex-direction: column;
        align-items: center;
        font-size: 10px;
        color: #7a92a8;
        line-height: 1.3;
        text-align: center;
        /* match editor-grid column widths so items line up */
      }
      .overwrite-item:nth-child(2)  { min-width: 140px; }
      .overwrite-item:nth-child(3)  { min-width: 170px; }
      .overwrite-item:nth-child(4)  { min-width: 230px; }
      .overwrite-item:nth-child(5)  { min-width: 170px; }
      .overwrite-item:nth-child(6)  { min-width: 160px; }
      .overwrite-item:nth-child(7)  { min-width: 110px; }
      .overwrite-item:nth-child(8)  { min-width: 110px; }
      .overwrite-item:nth-child(9)  { min-width: 110px; }
      .overwrite-item input[type=checkbox] { cursor: pointer; margin-bottom: 2px; }

      .editor-apply-row { margin-top: 12px; display: flex; align-items: center; gap: 14px; }
      .row-editor-hint { font-size: 11px; color: #7a92a8; font-style: italic; margin: 0; }

      /* DT */
      .as-dt-wrapper .dataTables_wrapper { font-family: 'DM Sans', sans-serif; }
      .as-dt-wrapper table.dataTable thead th { font-family: 'DM Mono', monospace;
        font-size: 11px; font-weight: 500; letter-spacing: 0.04em; background: #f0f4f7;
        color: #3a5068; border-bottom: 2px solid #dde3ea !important; }
      .as-dt-wrapper table.dataTable tbody td { font-family: 'DM Mono', monospace;
        font-size: 11.5px; vertical-align: middle; padding: 5px 10px; }
      .as-dt-wrapper table.dataTable tbody tr.selected td { background: #e0f4fa !important; color: #003d52 !important; }
      .as-dt-wrapper table.dataTable tbody tr:hover td { background: #f5fbfd; }
      .as-dt-wrapper .dataTables_info { font-size: 11.5px; color: #7a92a8; margin-top: 8px; }
      .as-dt-wrapper table.dataTable td.dt-center { text-align: center; }

      #check_output { font-family: 'DM Mono', monospace; font-size: 12px;
        background: #f8f9fb; border: 1px solid #dde3ea; border-radius: 7px;
        padding: 14px 18px; min-height: 60px; color: #2d4356; white-space: pre-wrap;
        line-height: 1.7; display: block; width: 100%; }
      #check_output.has-errors   { border-color: #e74c3c; background: #fff9f9; }
      #check_output.has-warnings { border-color: #f39c12; background: #fffdf5; }
      #check_output.all-ok       { border-color: #27ae60; background: #f4fef8; }

      .modal-content { border-radius: 10px; }
      .modal-header  { background: #0f1923; color: #fff; border-radius: 10px 10px 0 0; }
      .modal-title   { font-family: 'DM Mono', monospace; font-size: 15px; }
      .shiny-notification { font-family: 'DM Sans', sans-serif !important;
        border-radius: 8px !important; font-size: 13px !important; }
    "))
  ),

  useShinyjs(),
  add_busy_spinner(spin = "fading-circle", color = "#00b4d8", position = "bottom-right"),

  div(class = "as-banner",
      h1("AutoSpectral"),
      span(class = "as-subtitle", "Control File Builder")
  ),

  div(class = "as-outer",

      div(class = "as-sidebar",
          div(class = "as-section-label", "Setup"),
          selectInput("cytometer", "Cytometer", choices = names(map), selected = names(map)[1]),
          div(class = "as-section-label", "Control directory"),
          shinyDirButton("control_dir", label = "Browse...",
                         title = "Select control directory", class = "btn-dir"),
          div(class = "dir-display", textOutput("selected_dir")),
          actionButton("load_controls", "Load controls", class = "btn-as-primary"),
          hr(class = "as-hr"),
          div(class = "as-section-label", "Export"),
          textInput("output_filename", label = "Output filename", value = "fcs_control_file.csv"),
          selectInput("save_location", "Save to",
                      choices = c("Working directory" = "wd", "Control directory" = "fcs"),
                      selected = "wd"),
          actionButton("save_controls", "Save control file", class = "btn-as-save"),
          hr(class = "as-hr"),
          div(class = "as-section-label", "Instructions"),
          tags$ul(class = "as-instructions",
                  tags$li("Select cytometer and browse to FCS directory."),
                  tags$li("Click 'Load controls' to auto-fill the table."),
                  tags$li("Select one or more rows to open the editor."),
                  tags$li("When multiple rows are selected, fields that differ are greyed and inactive. Tick the corresponding overwrite checkbox to activate them."),
                  tags$li("Free-text cells (fluorophore, marker, gate.name) can also be double-clicked to edit inline."),
                  tags$li("Save — validation runs before writing.")
          )
      ),

      div(class = "as-main",

          div(class = "as-card",
              div(class = "as-card-header", span(class = "dot"), "Control table"),

              div(class = "tbl-toolbar",
                  actionButton("row_up",                  "Move up",             class = "btn-tool"),
                  actionButton("row_down",                "Move down",           class = "btn-tool"),
                  actionButton("add_row",                 "Add row",             class = "btn-tool"),
                  actionButton("delete_row",              "Delete selected",     class = "btn-tool-danger"),
                  actionButton("fill_control_type",       "Fill control.type",   class = "btn-tool"),
                  actionButton("fill_universal_negative", "Fill universal.neg.", class = "btn-tool"),
                  actionButton("fill_gate_names",         "Auto gate names",     class = "btn-tool"),
                  actionButton("clear_table",             "Clear all",           class = "btn-tool-danger")
              ),

              div(class = "as-dt-wrapper", DTOutput("controls_table")),

              # ---- Row editor ----
              conditionalPanel(
                condition = "output.row_selected",
                div(class = "row-editor",

                    # Editor fields in a fixed-column grid
                    div(class = "editor-grid",

                        div(class = "editor-field", id = "fld_control_type",
                            tags$label("control.type"),
                            selectInput("edit_control_type", label = NULL,
                                        choices = c("", "cells", "beads"), width = "100%")
                        ),
                        div(class = "editor-field", id = "fld_channel",
                            tags$label("channel"),
                            selectInput("edit_channel", label = NULL, choices = c(""), width = "100%")
                        ),
                        div(class = "editor-field", id = "fld_universal_neg",
                            tags$label("universal.negative"),
                            selectInput("edit_universal_neg", label = NULL, choices = c(""), width = "100%")
                        ),
                        div(class = "editor-field", id = "fld_gate_name",
                            tags$label("gate.name"),
                            selectizeInput("edit_gate_name", label = NULL, choices = c(""), width = "100%",
                                           options = list(create = TRUE, placeholder = "Select or type..."))
                        ),
                        div(class = "editor-field", id = "fld_is_unstained",
                            tags$label("Is unstained control"),
                            checkboxInput("edit_is_unstained", label = NULL, value = FALSE)
                        ),
                        div(class = "editor-field", id = "fld_large_gate",
                            tags$label("large.gate"),
                            checkboxInput("edit_large_gate", label = NULL, value = FALSE)
                        ),
                        div(class = "editor-field", id = "fld_is_viability",
                            tags$label("is.viability"),
                            checkboxInput("edit_is_viability", label = NULL, value = FALSE)
                        ),
                        div(class = "editor-field", id = "fld_gate_define",
                            tags$label("gate.define"),
                            checkboxInput("edit_gate_define", label = NULL, value = TRUE)
                        )
                    ),

                    # Overwrite row — single horizontal line, columns aligned to grid above
                    conditionalPanel(
                      condition = "output.multi_row_selected",
                      div(class = "overwrite-row",
                          span(class = "overwrite-label-hdr", "Overwrite mixed:"),
                          div(class = "overwrite-item",
                              checkboxInput("ow_control_type",  NULL, value = FALSE),
                              "control.type"),
                          div(class = "overwrite-item",
                              checkboxInput("ow_channel",       NULL, value = FALSE),
                              "channel"),
                          div(class = "overwrite-item",
                              checkboxInput("ow_universal_neg", NULL, value = FALSE),
                              "universal.neg"),
                          div(class = "overwrite-item",
                              checkboxInput("ow_gate_name",     NULL, value = FALSE),
                              "gate.name"),
                          div(class = "overwrite-item",
                              checkboxInput("ow_is_unstained",  NULL, value = FALSE),
                              "unstained"),
                          div(class = "overwrite-item",
                              checkboxInput("ow_large_gate",    NULL, value = FALSE),
                              "large.gate"),
                          div(class = "overwrite-item",
                              checkboxInput("ow_is_viability",  NULL, value = FALSE),
                              "viability"),
                          div(class = "overwrite-item",
                              checkboxInput("ow_gate_define",   NULL, value = FALSE),
                              "gate.define")
                      )
                    ),

                    div(class = "editor-apply-row",
                        actionButton("apply_edit", "Apply to selected row(s)",
                                     class = "btn-as-primary",
                                     style = "width:auto; padding: 6px 16px;"),
                        p(class = "row-editor-hint",
                          "Greyed fields differ across selected rows and won't be applied unless the overwrite checkbox is ticked.")
                    )
                )
              )
          ),

          div(class = "as-card",
              div(class = "as-card-header",
                  span(class = "dot", style = "background:#06d6a0;"), "Validation results"),
              verbatimTextOutput("check_output")
          )
      )
  ),

  tags$script(HTML("
    Shiny.addCustomMessageHandler('setCheckClass', function(cls) {
      var el = document.getElementById('check_output');
      if (el) el.className = cls;
    });
  "))
)

# ============================================================
# ---- Server ----
# ============================================================
server <- function(input, output, session) {

  roots <- .roots
  shinyDirChoose(input, "control_dir", roots = roots, session = session,
                 allowDirCreate = FALSE)

  fcs_dir <- reactiveVal(NULL)
  observeEvent(input$control_dir, {
    p <- parseDirPath(roots, input$control_dir)
    if (length(p) > 0) fcs_dir(as.character(p))
  })

  output$selected_dir <- renderText({
    d <- fcs_dir(); if (is.null(d)) "No directory selected" else d
  })

  asp <- reactiveVal(
    tryCatch(AutoSpectral::get.autospectral.param(map[[names(map)[1]]]),
             error = function(e) NULL))
  observeEvent(input$cytometer, {
    p <- tryCatch(AutoSpectral::get.autospectral.param(map[[input$cytometer]]),
                  error = function(e) NULL)
    asp(p)
  })

  rv <- reactiveValues(tbl = NULL, detectors = character(0))

  # ---- Load controls ----
  observeEvent(input$load_controls, {
    req(fcs_dir(), asp())
    show_modal_spinner(spin = "fading-circle", color = "#00b4d8",
                       text = "Matching fluorophores and markers...")
    on.exit(remove_modal_spinner())
    result <- tryCatch(
      create.control.df(fcs_dir(), asp()),
      error = function(e) { showNotification(e$message, type = "error"); NULL })
    if (is.null(result)) return()
    rv$tbl       <- result$df
    rv$detectors <- result$detectors
    updateSelectInput(session, "edit_channel", choices = c("", result$detectors))
    showNotification(paste0(nrow(result$df), " control files loaded."), type = "message")
  })

  # ---- Render DT ----
  output$controls_table <- renderDT({
    req(rv$tbl)
    df        <- rv$tbl
    bool_cols <- c("large.gate", "is.viability", "gate.define", "is_unstained")
    bool_cols <- bool_cols[bool_cols %in% colnames(df)]

    display <- df
    for (cn in bool_cols)
      display[[cn]] <- ifelse(is.na(display[[cn]]), FALSE, display[[cn]])

    editor_cols <- c("channel", "control.type", "universal.negative",
                     "large.gate", "is.viability", "gate.define", "is_unstained")
    disable_idx <- which(colnames(df) %in% editor_cols) - 1L
    bool_idx    <- which(colnames(df) %in% bool_cols) - 1L

    bool_defs <- lapply(as.list(bool_idx), function(i)
      list(targets = i, width = "55px", className = "dt-center",
           render = JS("function(data,type,row,meta){",
                       "  if(type!=='display') return data;",
                       "  return data ? '&#10003;' : '&#10007;';",
                       "}")))

    col_defs <- c(
      list(list(targets = 0, width = "190px")),
      list(list(targets = disable_idx, className = "text-muted")),
      bool_defs)

    datatable(display,
              selection = list(mode = "multiple", target = "row"),
              editable  = list(target = "cell", disable = list(columns = disable_idx)),
              rownames  = FALSE, escape = TRUE,
              options   = list(dom = "tip", scrollX = TRUE, scrollY = "240px",
                               scrollCollapse = TRUE, paging = FALSE, ordering = FALSE,
                               columnDefs = col_defs),
              class = "stripe hover compact")
  }, server = TRUE)

  # ---- Free-text cell edits ----
  observeEvent(input$controls_table_cell_edit, {
    info <- input$controls_table_cell_edit
    df   <- rv$tbl
    row  <- info$row
    col  <- info$col + 1L
    if (col < 1 || col > ncol(df)) return()
    df[row, col] <- info$value
    rv$tbl <- df
  })

  # ---- Selection state ----
  output$row_selected <- reactive({
    !is.null(input$controls_table_rows_selected) &&
      length(input$controls_table_rows_selected) > 0
  })
  outputOptions(output, "row_selected", suspendWhenHidden = FALSE)

  output$multi_row_selected <- reactive({
    !is.null(input$controls_table_rows_selected) &&
      length(input$controls_table_rows_selected) > 1
  })
  outputOptions(output, "multi_row_selected", suspendWhenHidden = FALSE)

  # ---- Helper: uniform check ----
  # For logicals, NA and FALSE are both treated as "not set".
  .uniform <- function(df, sel, col) {
    vals <- df[[col]][sel]
    if (is.logical(vals)) vals <- ifelse(is.na(vals), FALSE, vals)
    length(unique(vals)) <= 1
  }

  # ---- Populate editor when selection changes ----
  observeEvent(input$controls_table_rows_selected, {
    sel <- input$controls_table_rows_selected
    if (is.null(sel) || length(sel) == 0) return()
    df  <- rv$tbl
    r   <- sel[1]

    for (id in vapply(.field_specs, `[[`, character(1), "ow"))
      updateCheckboxInput(session, id, value = FALSE)

    updateSelectInput(session, "edit_control_type",
                      selected = if (is.na(df$control.type[r])) "" else df$control.type[r])
    updateSelectInput(session, "edit_channel",
                      selected = if (is.na(df$channel[r])) "" else df$channel[r])

    ct <- df$control.type[r]
    un_choices <- c("",
                    df$filename[!is.na(df$is_unstained) & df$is_unstained &
                                  (df$control.type == ct | is.na(df$control.type))])
    updateSelectInput(session, "edit_universal_neg", choices = un_choices,
                      selected = if (is.na(df$universal.negative[r])) ""
                      else df$universal.negative[r])

    existing_gates <- sort(unique(df$gate.name[!is.na(df$gate.name) & df$gate.name != ""]))
    current_gate   <- if (is.na(df$gate.name[r]) || df$gate.name[r] == "") "" else df$gate.name[r]
    updateSelectizeInput(session, "edit_gate_name",
                         choices = unique(c("", existing_gates, current_gate)),
                         selected = current_gate, server = FALSE)

    updateCheckboxInput(session, "edit_is_unstained", value = isTRUE(df$is_unstained[r]))
    updateCheckboxInput(session, "edit_large_gate",   value = isTRUE(df$large.gate[r]))
    updateCheckboxInput(session, "edit_is_viability", value = isTRUE(df$is.viability[r]))
    updateCheckboxInput(session, "edit_gate_define",  value = isTRUE(df$gate.define[r]))
  })

  # ---- Grey mixed fields via shinyjs ----
  # Runs whenever selection or any overwrite checkbox changes.
  observe({
    sel <- input$controls_table_rows_selected
    df  <- rv$tbl
    # read all overwrite inputs to register dependencies
    ow_vals <- setNames(
      vapply(.field_specs, function(s) isTRUE(input[[s$ow]]), logical(1)),
      vapply(.field_specs, `[[`, character(1), "div"))

    if (is.null(sel) || length(sel) <= 1 || is.null(df)) {
      for (spec in .field_specs) shinyjs::removeClass(spec$div, "field-mixed")
      return()
    }
    for (spec in .field_specs) {
      uni   <- .uniform(df, sel, spec$col)
      mixed <- !uni && !ow_vals[[spec$div]]
      if (mixed) shinyjs::addClass(spec$div, "field-mixed")
      else       shinyjs::removeClass(spec$div, "field-mixed")
    }
  })

  # ---- Apply editor ----
  # Evaluate all field decisions ONCE on the unmodified df BEFORE the row loop,
  # so mutations to df inside the loop don't affect uniformity checks.
  observeEvent(input$apply_edit, {
    req(rv$tbl)
    sel <- input$controls_table_rows_selected
    if (is.null(sel) || length(sel) == 0) {
      showNotification("No rows selected.", type = "warning"); return()
    }
    df    <- rv$tbl
    multi <- length(sel) > 1

    decide <- function(col, ow_id) {
      if (!multi) return(TRUE)
      .uniform(df, sel, col) || isTRUE(input[[ow_id]])
    }

    do_control_type  <- decide("control.type",       "ow_control_type")
    do_channel       <- decide("channel",            "ow_channel")
    do_universal_neg <- decide("universal.negative", "ow_universal_neg")
    do_gate_name     <- decide("gate.name",          "ow_gate_name")
    do_is_unstained  <- decide("is_unstained",       "ow_is_unstained")
    do_large_gate    <- decide("large.gate",         "ow_large_gate")
    do_is_viability  <- decide("is.viability",       "ow_is_viability")
    do_gate_define   <- decide("gate.define",        "ow_gate_define")

    for (r in sel) {
      if (do_control_type)  df$control.type[r]       <- input$edit_control_type
      if (do_channel)       df$channel[r]            <- input$edit_channel
      if (do_universal_neg) df$universal.negative[r] <- input$edit_universal_neg
      if (do_gate_name)     df$gate.name[r]          <- if (is.null(input$edit_gate_name)) "" else input$edit_gate_name
      if (do_is_unstained)  df$is_unstained[r]       <- isTRUE(input$edit_is_unstained)
      if (do_large_gate)    df$large.gate[r]         <- isTRUE(input$edit_large_gate)
      if (do_is_viability)  df$is.viability[r]       <- isTRUE(input$edit_is_viability)
      if (do_gate_define)   df$gate.define[r]        <- isTRUE(input$edit_gate_define)
    }

    rv$tbl <- df

    ct <- df$control.type[sel[1]]
    un_choices <- c("",
                    df$filename[!is.na(df$is_unstained) & df$is_unstained &
                                  (df$control.type == ct | is.na(df$control.type))])
    updateSelectInput(session, "edit_universal_neg", choices = un_choices)

    existing_gates <- sort(unique(df$gate.name[!is.na(df$gate.name) & df$gate.name != ""]))
    updateSelectizeInput(session, "edit_gate_name",
                         choices = unique(c("", existing_gates)), server = FALSE)

    showNotification(paste0("Applied to ", length(sel), " row(s)."), type = "message")
  })

  # ---- Row move up ----
  observeEvent(input$row_up, {
    req(rv$tbl)
    sel <- input$controls_table_rows_selected
    if (is.null(sel) || length(sel) == 0) {
      showNotification("Select a row first.", type = "warning"); return()
    }
    r <- sel[1]; if (r <= 1) return()
    df <- rv$tbl; tmp <- df[r-1,]; df[r-1,] <- df[r,]; df[r,] <- tmp; rv$tbl <- df
  })

  # ---- Row move down ----
  observeEvent(input$row_down, {
    req(rv$tbl)
    sel <- input$controls_table_rows_selected
    if (is.null(sel) || length(sel) == 0) {
      showNotification("Select a row first.", type = "warning"); return()
    }
    r <- sel[1]; if (r >= nrow(rv$tbl)) return()
    df <- rv$tbl; tmp <- df[r+1,]; df[r+1,] <- df[r,]; df[r,] <- tmp; rv$tbl <- df
  })

  # ---- Add row ----
  observeEvent(input$add_row, {
    req(rv$tbl); rv$tbl <- rbind(rv$tbl, .empty.row(rv$tbl))
  })

  # ---- Delete selected rows ----
  observeEvent(input$delete_row, {
    req(rv$tbl)
    sel <- input$controls_table_rows_selected
    if (is.null(sel) || length(sel) == 0) {
      showNotification("Select one or more rows first.", type = "warning"); return()
    }
    showModal(modalDialog(title = "Delete rows",
                          paste0("Delete ", length(sel), " selected row(s)?"),
                          footer = tagList(modalButton("Cancel"),
                                           actionButton("confirm_delete", "Delete", class = "btn-danger"))))
  })
  observeEvent(input$confirm_delete, {
    req(rv$tbl)
    rv$tbl <- rv$tbl[-input$controls_table_rows_selected, , drop = FALSE]
    removeModal()
  })

  # ---- Fill control.type ----
  observeEvent(input$fill_control_type, {
    req(rv$tbl)
    showModal(modalDialog(title = "Fill empty control.type",
                          radioButtons("fill_ct_choice", "Fill empty cells with:", choices = c("cells", "beads")),
                          footer = tagList(modalButton("Cancel"),
                                           actionButton("confirm_fill_ct", "Fill", class = "btn-primary"))))
  })
  observeEvent(input$confirm_fill_ct, {
    req(rv$tbl)
    df  <- rv$tbl
    idx <- which(is.na(df$control.type) | df$control.type == "")
    if (length(idx) > 0) {
      df$control.type[idx] <- input$fill_ct_choice; rv$tbl <- df
      showNotification(paste("Filled", length(idx), "rows."), type = "message")
    } else {
      showNotification("No empty control.type cells.", type = "warning")
    }
    removeModal()
  })

  # ---- Fill universal.negative ----
  observeEvent(input$fill_universal_negative, {
    req(rv$tbl)
    df  <- rv$tbl
    idx <- which(!is.na(df$is_unstained) & df$is_unstained)
    if (length(idx) == 0) {
      showNotification(
        "No rows marked as unstained. Select rows and tick 'Is unstained control' in the editor.",
        type = "warning"); return()
    }
    cell_negs <- df$filename[idx[df$control.type[idx] == "cells"]]
    bead_negs <- df$filename[idx[df$control.type[idx] == "beads"]]
    showModal(modalDialog(title = "Fill universal.negative",
                          if (length(cell_negs) > 0)
                            radioButtons("fill_un_cells", "Negative for CELLS:", choices = c("(leave blank)" = "", cell_negs)),
                          if (length(bead_negs) > 0)
                            radioButtons("fill_un_beads", "Negative for BEADS:", choices = c("(leave blank)" = "", bead_negs)),
                          footer = tagList(modalButton("Cancel"),
                                           actionButton("confirm_fill_un", "Fill", class = "btn-primary"))))
  })
  observeEvent(input$confirm_fill_un, {
    req(rv$tbl)
    df      <- rv$tbl
    stained <- which(is.na(df$is_unstained) | !df$is_unstained)
    cells_c <- if (!is.null(input$fill_un_cells)) input$fill_un_cells else ""
    beads_c <- if (!is.null(input$fill_un_beads)) input$fill_un_beads else ""
    if (cells_c != "")
      df$universal.negative[stained[df$control.type[stained] == "cells"]] <- cells_c
    if (beads_c != "")
      df$universal.negative[stained[df$control.type[stained] == "beads"]] <- beads_c
    rv$tbl <- df; removeModal()
    showNotification("universal.negative filled.", type = "message")
  })

  # ---- Auto-fill gate names ----
  # Key design decisions from define.flow.control():
  #   1. control.table$sample is set to fluorophore BEFORE assign.gates() runs.
  #      assign.gates() then modifies $sample for replicated neg rows but never
  #      touches $fluorophore, $marker, or $universal.negative in original rows.
  #   2. assign.gates() skips rows where gate.name is already set (not NA).
  #      So we must clear all gate names first to allow full reassignment when
  #      large.gate / is.viability etc. change.
  #   3. We write back ONLY gate.name and gate.define from result.
  #      For extra rows appended by assign.gates() (replicated negatives), we
  #      copy only structural columns -- NOT fluorophore/marker/universal.negative,
  #      which assign.gates() may have renamed with descriptive gate suffixes.
  observeEvent(input$fill_gate_names, {
    req(rv$tbl)
    df     <- rv$tbl
    n_orig <- nrow(df)

    # Build clean input for assign.gates, mirroring define.flow.control()
    df_ag <- df[, setdiff(colnames(df), "is_unstained"), drop = FALSE]

    # "" -> NA for assign.gates to handle correctly
    df_ag$gate.name[!is.na(df_ag$gate.name) & df_ag$gate.name == ""] <- NA
    df_ag$universal.negative[!is.na(df_ag$universal.negative) &
                               df_ag$universal.negative == ""] <- NA

    # Clear ALL existing gate names so assign.gates() does a full fresh reassignment.
    # If we leave existing names in, assign.gates() skips those rows and won't
    # update them when large.gate / is.viability has changed.
    df_ag$gate.name <- NA

    # Initialise $sample from $fluorophore, exactly as define.flow.control() does
    df_ag$sample <- df_ag$fluorophore

    result <- tryCatch(
      withCallingHandlers(
        AutoSpectral::assign.gates(df_ag, gating.system = "density",
                                   gate = TRUE, verbose = FALSE),
        warning = function(w) {
          message("fill_gate_names WARNING: ", conditionMessage(w))
          invokeRestart("muffleWarning")
        }
      ),
      error = function(e) {
        message("fill_gate_names ERROR: ", conditionMessage(e))
        showNotification(paste("assign.gates:", conditionMessage(e)), type = "error")
        NULL
      })

    if (is.null(result)) return()

    # Write back gate.name and gate.define for the original rows by position.
    # assign.gates() preserves input row order for original rows.
    n_use <- min(n_orig, nrow(result))
    if ("gate.name" %in% colnames(result)) {
      vals <- result$gate.name[seq_len(n_use)]
      vals[is.na(vals)] <- ""
      df$gate.name[seq_len(n_use)] <- vals
    }
    if ("gate.define" %in% colnames(result))
      df$gate.define[seq_len(n_use)] <- result$gate.define[seq_len(n_use)]


    # Append truly new extra rows (replicated negatives from assign.gates()).
    # Check filename+gate.name to avoid duplicates when re-running.
    # Name them "Negative N" -- unique, contains "negative" so clean.controls()
    # correctly excludes them from stained-sample processing.
    if (nrow(result) > n_orig) {
      extra      <- result[(n_orig + 1):nrow(result), , drop = FALSE]
      df_keys    <- paste(df$filename, df$gate.name, sep = "\t")
      extra_keys <- paste(extra$filename,
                          ifelse(is.na(extra$gate.name), "", extra$gate.name),
                          sep = "\t")
      truly_new  <- !extra_keys %in% df_keys
      if (any(truly_new)) {
        en       <- extra[truly_new, , drop = FALSE]
        new_rows <- .empty.row(df)[rep(1L, nrow(en)), ]

        safe_cols <- c("filename", "control.type", "channel",
                       "large.gate", "is.viability", "gate.name", "gate.define")
        for (cn in intersect(safe_cols, colnames(en))) new_rows[[cn]] <- en[[cn]]
        new_rows$gate.name[is.na(new_rows$gate.name)] <- ""

        # Unique "Negative N" names -- contains "negative" for clean.controls(),
        # unique so fluorophore uniqueness validation still passes.
        existing_neg_n <- sum(grepl("^Negative \\d+$", df$fluorophore))
        new_rows$fluorophore <- paste0("Negative ",
                                       seq_len(nrow(new_rows)) + existing_neg_n)
        # Mark as unstained so they appear in the universal.negative dropdown
        new_rows$is_unstained <- TRUE
        df <- rbind(df, new_rows)
      }
    }

    rv$tbl <- df
    n_named <- sum(!is.na(df$gate.name) & df$gate.name != "", na.rm = TRUE)
    showNotification(paste0("Gate names assigned to ", n_named, " rows."),
                     type = "message")
  })

  # ---- Clear table ----
  observeEvent(input$clear_table, {
    showModal(modalDialog(title = "Clear table", "Remove all rows?",
                          footer = tagList(modalButton("Cancel"),
                                           actionButton("confirm_clear", "Clear", class = "btn-danger"))))
  })
  observeEvent(input$confirm_clear, {
    rv$tbl <- NULL; removeModal()
    showNotification("Table cleared.", type = "warning")
  })

  # ---- Validation display ----
  check_result_text <- reactiveVal("No checks performed yet.")
  check_has_errors  <- reactiveVal(FALSE)
  check_has_warn    <- reactiveVal(FALSE)
  check_all_ok      <- reactiveVal(FALSE)

  output$check_output <- renderText({ check_result_text() })
  outputOptions(output, "check_output", suspendWhenHidden = FALSE)

  observe({
    cls <- if      (check_has_errors()) "has-errors"
    else if (check_has_warn())   "has-warnings"
    else if (check_all_ok())     "all-ok"
    else                         ""
    session$sendCustomMessage("setCheckClass", cls)
  })

  # ---- Save ----
  observeEvent(input$save_controls, {
    req(rv$tbl, fcs_dir())
    df      <- prepare.for.export(rv$tbl)
    tmpfile <- tempfile(fileext = ".csv")
    write.csv(df, tmpfile, row.names = FALSE, na = "")

    raw_res <- tryCatch(
      suppressMessages(suppressWarnings(
        AutoSpectral::check.control.file(
          control.def.file = tmpfile, control.dir = fcs_dir(),
          asp = asp(), strict = FALSE))),
      error = function(e)
        data.frame(severity = "fatal", rule = "execution_error",
                   message = e$message, stringsAsFactors = FALSE))

    check_res <- normalize_check_result(raw_res)
    errors    <- subset(check_res, severity == "error")
    warnings  <- subset(check_res, severity == "warning")
    fatals    <- subset(check_res, severity == "fatal")

    txt <- if (nrow(fatals) > 0) {
      paste0("FATAL ERROR:\n", paste(fatals$message, collapse = "\n"))
    } else if (nrow(errors) == 0 && nrow(warnings) == 0) {
      "No validation issues found."
    } else {
      parts <- character(0)
      if (nrow(errors) > 0)
        parts <- c(parts, paste0("ERRORS:\n",
                                 paste(paste0("  - ", unique(errors$message)), collapse = "\n")))
      if (nrow(warnings) > 0)
        parts <- c(parts, paste0("WARNINGS:\n",
                                 paste(paste0("  - ", unique(warnings$message)), collapse = "\n")))
      paste(parts, collapse = "\n\n")
    }

    check_result_text(txt)
    check_has_errors(nrow(fatals) > 0 || nrow(errors) > 0)
    check_has_warn(nrow(warnings) > 0 && nrow(errors) == 0 && nrow(fatals) == 0)
    check_all_ok(nrow(check_res) == 0)

    if (nrow(fatals) > 0 || nrow(errors) > 0) {
      showNotification("Validation errors — file not saved.", type = "error"); return()
    }

    save.dir  <- if (input$save_location == "fcs") fcs_dir() else .launch_wd
    outname   <- trimws(input$output_filename)
    if (outname == "") outname <- "fcs_control_file.csv"
    if (!grepl("\\.csv$", outname, ignore.case = TRUE)) outname <- paste0(outname, ".csv")
    base      <- tools::file_path_sans_ext(file.path(save.dir, outname))
    final_out <- file.path(save.dir, outname)
    counter   <- 1L
    while (file.exists(final_out)) {
      final_out <- paste0(base, "_", counter, ".csv"); counter <- counter + 1L
    }
    write.csv(df, final_out, row.names = FALSE, na = "")
    showNotification(paste0("Saved: ", final_out), type = "message", duration = 8)
  })

}

# ============================================================
shinyApp(ui, server)
