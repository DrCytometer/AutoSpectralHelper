# app.R

# install any missing packages and start-up
required.packages <- c("shiny", "shinyFiles", "rhandsontable", "dplyr", "flowCore")
missing.packages <- setdiff(required.packages, installed.packages())
if(length(missing.packages) != 0)
  invisible(BiocManager::install(missing.packages))
if(!"AutoSpectral" %in% installed.packages()) {
  if(!"devtools" %in% installed.packages()) invisible(install.packages("devtools"))
  devtools::install_github('DrCytometer/AutoSpectral')
}
invisible(lapply(c(required.packages, "AutoSpectral"), require, character.only = TRUE))

# cytometer mapping
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

# helper function returning dataframe for user modification
create.control.df <- function(control.dir, asp) {
  # read fluorophore database
  data.path <- system.file("extdata", "fluorophore_database.csv", package = "AutoSpectral")
  fluorophore.database <- read.csv(data.path, stringsAsFactors = FALSE)
  fluorophore.database[fluorophore.database == ""] <- NA

  control.files <- list.files(control.dir, pattern = ".fcs", ignore.case = TRUE)
  if (length(control.files) < 1) stop("Single-stained control files not found.")

  control.colnames <- c("filename", "fluorophore", "marker", "channel",
                        "control.type", "universal.negative", "large.gate")
  control.def.file <- data.frame(matrix(ncol = length(control.colnames),
                                        nrow = length(control.files)),
                                 stringsAsFactors = FALSE)
  colnames(control.def.file) <- control.colnames
  control.def.file$filename <- control.files

  # match fluorophores to database
  fluorophore.matches <- suppressMessages(
    AutoSpectral::match.fluorophores(control.files, fluorophore.database)
  )
  control.def.file$fluorophore <- fluorophore.matches[control.def.file$filename]

  # channel mapping based on cytometer
  if (asp$cytometer == "Aurora") {
    if (asp$cytometer.version == "NL") {
      detectors <- setNames(fluorophore.database$channel.NL, fluorophore.database$fluorophore)
    } else {
      detectors <- setNames(fluorophore.database$channel.Aurora, fluorophore.database$fluorophore)
    }
  } else if (asp$cytometer == "ID7000") {
    detectors <- setNames(fluorophore.database$channel.ID7000, fluorophore.database$fluorophore)
  } else if (asp$cytometer == "FACSDiscover A8" || asp$cytometer == "FACSDiscover S8") {
    detectors <- setNames(fluorophore.database$channel.s8, fluorophore.database$fluorophore)
  } else if (asp$cytometer == "Opteon") {
    detectors <- setNames(fluorophore.database$channel.opteon, fluorophore.database$fluorophore)
  } else if (asp$cytometer == "Mosaic") {
    detectors <- setNames(fluorophore.database$channel.mosaic, fluorophore.database$fluorophore)
  } else if (asp$cytometer == "Xenith") {
    detectors <- setNames(fluorophore.database$channel.xenith, fluorophore.database$fluorophore)
  } else if (asp$cytometer == "Symphony" || asp$cytometer == "A5SE") {
    detectors <- setNames(fluorophore.database$channel.A5SE, fluorophore.database$fluorophore)
  } else stop("Unsupported cytometer")

  detector.idx <- match(control.def.file$fluorophore, names(detectors))
  control.def.file$channel <- detectors[detector.idx]

  # control.type
  control.def.file$control.type <- sapply(control.def.file$filename, function(filename) {
    if (grepl("cells", filename, ignore.case = TRUE)) "cells"
    else if (grepl("beads", filename, ignore.case = TRUE)) "beads"
    else ""
  })

  # markers
  marker.data.path <- system.file("extdata", "marker_database.csv", package = "AutoSpectral")
  marker.database <- read.csv(marker.data.path)
  marker.database[marker.database == ""] <- NA
  marker.matches <- suppressMessages(match.markers(control.files, marker.database))
  control.def.file$marker <- marker.matches[control.def.file$filename]

  # merge and reorder
  control.def.file.merged <- merge(control.def.file, fluorophore.database,
                                   by = "fluorophore", all.x = TRUE)
  laser.order <- c("DeepUV", "UV", "Violet", "Blue", "YellowGreen", "Red", "IR")
  control.def.file.merged$excitation.laser <- factor(control.def.file.merged$excitation.laser,
                                                     levels = laser.order)
  control.def.file.merged <- control.def.file.merged[order(control.def.file.merged$excitation.laser,
                                                           control.def.file.merged$nominal.wavelength), ]
  desired.col <- c(control.colnames, "is.viability")
  control.def.file <- control.def.file.merged[, desired.col]

  # fill AF / Negative
  idx_unstained <- grepl("Unstained", control.def.file$filename, ignore.case = TRUE)
  control.def.file$fluorophore[idx_unstained] <- ifelse(control.def.file$control.type[idx_unstained] == "cells",
                                                        "AF", "Negative")
  control.def.file$fluorophore[grepl("Negative", control.def.file$filename, ignore.case = TRUE)] <- "Negative"

  control.def.file[is.na(control.def.file)] <- ""
  control.def.file[control.def.file == "NA"] <- ""

  # helper boolean column
  control.def.file$is_unstained <- FALSE

  # detectors for selection
  database.path <- system.file("extdata", "cytometer_database.csv", package = "AutoSpectral")
  cytometers <- read.csv(database.path)
  if (asp$cytometer == "Aurora") {
    detectors <- if (asp$cytometer.version == "NL") cytometers$NorthernLights else cytometers$Aurora
  } else if (asp$cytometer == "ID7000") {
    detectors <- cytometers$ID7000
  } else if (grepl("Discover", asp$cytometer)) {
    detectors <- cytometers$Discover
  } else if (asp$cytometer == "Opteon") {
    detectors <- cytometers$Opteon
  } else if (asp$cytometer == "Mosaic") {
    detectors <- cytometers$Mosaic
  } else if (asp$cytometer == "Xenith") {
    detectors <- cytometers$Xenith
  } else if (asp$cytometer == "Symphony") {
    detectors <- cytometers$A5SE
  } else stop("Unsupported cytometer")

  return(list(df = control.def.file, detectors = detectors))
}

# ---- UI ----
ui <- fluidPage(
  titlePanel("AutoSpectral: Interactive Control File Builder"),
  sidebarLayout(
    sidebarPanel(
      h4("Instructions"),
      tags$ul(
        tags$li("Select your cytometer."),
        tags$li("Choose the control directory."),
        tags$li("Click 'Load controls' to import files."),
        tags$li("Edit the control table as needed."),
        tags$li("Save the control file when you're done."),
        tags$li("Saving will run checks."),
        tags$li("Errors must be fixed if there are any."),
      ),
      selectInput("cytometer", "Cytometer", choices = names(map), selected = names(map)[1]),
      shinyDirButton("control_dir", "Select control directory", "Please select a folder"),
      verbatimTextOutput("selected_dir"),
      actionButton("load_controls", "Load controls from directory"),
      hr(),
      textInput("output_filename", "Output filename", value = "fcs_control_file.csv"),
      actionButton("save_controls", "Save control file"),
      width = 3
    ),
    mainPanel(
      h4("Control table"),
      helpText("Edit cells directly. Right-click to remove or add rows."),
      rHandsontableOutput("controls_table", height = "500px"),
      fluidRow(
        column(3, actionButton("fill_control_type", "Fill empty control.type")),
        column(6, actionButton("fill_universal_negative", "Fill all universal.negative"))
      ),
      hr(),
      h4("Check results"),
      verbatimTextOutput("check_output"),
      width = 9
    )
  )
)

# ---- Server ----
server <- function(input, output, session) {
  # validation helper function: allow compatibility with AutoSpectral v0.8.7 and earlier
  normalize_check_result <- function(x) {
    
    # Case 1: NULL → empty validation table
    if (is.null(x)) {
      return(data.frame(
        severity = character(0),
        rule     = character(0),
        message  = character(0),
        stringsAsFactors = FALSE
      ))
    }
    
    # Case 2: already a data.frame
    if (is.data.frame(x)) {
      required <- c("severity", "rule", "message")
      missing  <- setdiff(required, colnames(x))
      
      if (length(missing) > 0) {
        for (m in missing) {
          x[[m]] <- NA_character_
        }
      }
      
      return(x[, required, drop = FALSE])
    }
    
    # Case 3: legacy list return (pre-refactor)
    if (is.list(x)) {
      return(data.frame(
        severity = "error",
        rule     = "legacy_return",
        message  = "check.control.file() returned legacy list format",
        stringsAsFactors = FALSE
      ))
    }
    
    # Case 4: anything else (should never happen)
    data.frame(
      severity = "fatal",
      rule     = "unknown_return_type",
      message  = paste("Unexpected return type:", class(x)[1]),
      stringsAsFactors = FALSE
    )
  }
  
  # reactive
  roots <- c(WorkingDir = getwd(), Home = normalizePath("~"))
  shinyDirChoose(input, "control_dir", roots = roots, session = session)
  fcs_dir  <- reactiveVal(NULL)
  save_dir <- reactiveVal(NULL)
  observeEvent(input$control_dir, {
    dirinfo <- parseDirPath(roots, input$control_dir)
    if (length(dirinfo) > 0) {
      fcs_dir(dirinfo)
      save_dir(getwd())  # preserve legacy behavior: default save = FCS dir
    }
  })
  output$selected_dir <- renderText({
    d <- fcs_dir();
    if (is.null(d)) "No directory selected" else d
    })

  asp <- reactiveVal(AutoSpectral::get.autospectral.param(map[[names(map)[1]]]))
  observeEvent(input$cytometer, {
    selected <- input$cytometer
    aspp <- tryCatch(AutoSpectral::get.autospectral.param(map[[selected]]), error = function(e) NULL)
    asp(aspp)
  })

  rv <- reactiveValues(tbl = NULL, detectors = character(0), unstained_choices = character(0))
  
  # load controls
  observeEvent(input$load_controls, {
    req(fcs_dir(), asp())
    res <- create.control.df(fcs_dir(), asp())
    df <- res$df
    detectors <- res$detectors
    if (!"is.viability" %in% colnames(df)) df$is.viability <- ""
    if (!"universal.negative" %in% colnames(df)) df$universal.negative <- ""
    df$large.gate <- ifelse(df$large.gate == "", NA, df$large.gate)
    df$is.viability <- ifelse(df$is.viability == "", NA, df$is.viability)
    rv$tbl <- df
    rv$detectors <- detectors
    rv$unstained_choices <- character(0)
    showNotification("Control table created from directory", type = "message")
  })

  # update rv$tbl when table edited
  observeEvent(input$controls_table, {
    req(input$controls_table)
    tbl <- hot_to_r(input$controls_table)
    if (!"is_unstained" %in% colnames(tbl)) tbl$is_unstained <- FALSE
    if (!"universal.negative" %in% colnames(tbl)) tbl$universal.negative <- ""

    # auto-fill ALL rows with universal.negative if first row is changed to a non-empty value
    if (!is.null(input$controls_table_changes)) {
      changes <- input$controls_table_changes
      for (ch in changes) {
        col_name <- names(tbl)[ch$col + 1]
        if (col_name == "universal.negative" && ch$row == 0 && ch$newValue != "") {
          # fill all subsequent rows (not just empty ones)
          if (nrow(tbl) > 1) {
            tbl$universal.negative[2:nrow(tbl)] <- ch$newValue
          }
        }
      }
    }

    rv$tbl <- tbl
  })

  # render rhandsontables
  output$controls_table <- renderRHandsontable({
    req(rv$tbl)
    df <- rv$tbl

    channel_choices <- if(length(rv$detectors)>0) rv$detectors else sort(unique(df$channel[df$channel!=""]))
    control_type_choices <- c("cells","beads")
    bool_choices <- c("TRUE","FALSE")

    # get all unstained options - will show all regardless of control.type
    unstained_idx <- which(df$is_unstained == TRUE)
    if(length(unstained_idx) > 0) {
      universal_negative_choices <- c("", unique(df$filename[unstained_idx]))
    } else {
      universal_negative_choices <- c("")
    }

    rhandsontable(df, useTypes=TRUE, stretchH="all") %>%
      hot_table(highlightCol=TRUE, highlightRow=TRUE, allowRowEdit=TRUE) %>%
      hot_col("filename", readOnly=TRUE) %>%
      hot_col("fluorophore", type="text") %>%
      hot_col("marker", type="text") %>%
      hot_col("channel", type="dropdown", source=channel_choices) %>%
      hot_col("control.type", type="dropdown", source=control_type_choices) %>%
      hot_col("universal.negative", type = "dropdown", source = universal_negative_choices) %>%
      hot_col("large.gate", type="dropdown", source=bool_choices) %>%
      hot_col("is.viability", type="dropdown", source=bool_choices) %>%
      hot_col("is_unstained", type="checkbox")
  })

  # fill empty control.type
  observeEvent(input$fill_control_type, {
    showModal(modalDialog(
      title="Fill empty control.type fields",
      radioButtons("fill_choice","Choose value to fill empty control.type with:", choices=c("cells","beads")),
      footer=tagList(modalButton("Cancel"), actionButton("confirm_fill","Fill"))
    ))
  })
  observeEvent(input$confirm_fill, {
    req(rv$tbl)
    choice <- input$fill_choice
    df <- rv$tbl
    empty_idx <- which(is.na(df$control.type)|df$control.type=="")
    if(length(empty_idx)>0){df$control.type[empty_idx]<-choice; rv$tbl<-df; removeModal(); showNotification(paste("Filled", length(empty_idx), "rows"), type="message")}
    else{removeModal(); showNotification("No empty control.type fields", type="message")}
  })

  # fill all universal.negative
  observeEvent(input$fill_universal_negative, {
    req(rv$tbl)
    df <- rv$tbl

    # get all unstained controls
    unstained_idx <- which(df$is_unstained == TRUE)

    if(length(unstained_idx) == 0) {
      showNotification("No controls are marked as unstained. Please check the is_unstained boxes first.", type="warning")
      return()
    }

    unstained_choices <- unique(df$filename[unstained_idx])

    showModal(modalDialog(
      title="Fill all universal.negative fields",
      radioButtons("fill_universal_choice","Choose unstained control to fill all rows with:",
                   choices=unstained_choices),
      footer=tagList(modalButton("Cancel"), actionButton("confirm_fill_universal","Fill All"))
    ))
  })

  observeEvent(input$confirm_fill_universal, {
    req(rv$tbl)
    choice <- input$fill_universal_choice
    df <- rv$tbl
    df$universal.negative <- choice
    rv$tbl <- df
    removeModal()
    showNotification(paste("Filled all rows with", choice), type="message")
  })

  # placeholder check output
  output$check_output <- renderText({"No checks performed yet."})
  observeEvent(input$save_controls, {
    req(rv$tbl)
    req(save_dir())
    req(fcs_dir())
    df <- rv$tbl

    # write temp csv
    tmpfile <- tempfile(fileext = ".csv")
    write.csv(df[, setdiff(colnames(df), "is_unstained")], tmpfile, row.names = FALSE, na = "")

    # Call check.control.file
    raw_check_res <- tryCatch({
      suppressMessages(
        AutoSpectral::check.control.file(
          control.def.file = tmpfile,
          control.dir = fcs_dir(),
          asp = asp(),
          strict = FALSE
        )
      )
    }, error = function(e) {
      data.frame(
        severity = "fatal",
        rule     = "execution_error",
        message  = e$message,
        stringsAsFactors = FALSE
      )
    })
    
    check_res <- normalize_check_result(raw_check_res)
    
    # split errors by severity
    errors   <- subset(check_res, severity == "error")
    warnings <- subset(check_res, severity == "warning")
    fatals   <- subset(check_res, severity == "fatal")
    
    # format and show check results
    formatted <- capture.output({
      
      if (nrow(fatals) > 0) {
        cat("Fatal error:\n")
        print(fatals$message)
        return()
      }
      
      if (nrow(errors) == 0 && nrow(warnings) == 0) {
        cat("No validation issues found.")
        return()
      }
      
      if (nrow(errors) > 0) {
        cat("Errors:\n")
        print(unique(errors$message))
        cat("\n")
      }
      
      if (nrow(warnings) > 0) {
        cat("Warnings:\n")
        print(unique(warnings$message))
      }
    })
    output$check_output <- renderText(paste(formatted, collapse = "\n"))

    # if there are critical errors then do not save
    has_fatal  <- any(check_res$severity == "fatal")
    has_errors <- any(check_res$severity == "error")
    
    if (has_fatal || has_errors) {
      showNotification(
        "Critical errors detected — file not saved. See check results.",
        type = "error"
      )
      return()
    }

    # otherwise save to the chosen folder
    outname <- input$output_filename
    if (is.null(outname) || outname == "") outname <- "fcs_control_file.csv"
    req(save_dir())
    outpath <- file.path(save_dir(), outname)

    # if exists, generate unique filename
    base <- tools::file_path_sans_ext(outpath)
    ext <- tools::file_ext(outpath)
    if (ext == "") ext <- "csv"
    final_out <- outpath
    idx <- 1
    while (file.exists(final_out)) {
      final_out <- paste0(base, "_", idx, ".", ext)
      idx <- idx + 1
    }
    write.csv(df[, setdiff(colnames(df), "is_unstained")], final_out, row.names = FALSE, na = "")
    showNotification(paste("Control file saved to", final_out), type = "message")
  })
}

shinyApp(ui, server)
