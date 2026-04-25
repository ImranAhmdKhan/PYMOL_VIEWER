# =============================================================================
#  Molecular Viewer & Analytics — R / Shiny
#  Protein–Ligand Docking Analysis Dashboard
# =============================================================================

# ── Package bootstrap ─────────────────────────────────────────────────────────
required_pkgs <- c(
  "shiny", "shinydashboard", "shinyWidgets", "colourpicker",
  "bio3d",
  "ggplot2", "ggrepel",
  "plotly",
  "DT",
  "dplyr", "tidyr",
  "scales",
  "fmsb"          # radar / spider charts
)

missing_pkgs <- required_pkgs[!sapply(required_pkgs, requireNamespace, quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  message("Installing missing packages: ", paste(missing_pkgs, collapse = ", "))
  install.packages(missing_pkgs, repos = "https://cloud.r-project.org", quiet = TRUE)
}

suppressPackageStartupMessages({
  library(shiny)
  library(shinydashboard)
  library(shinyWidgets)
  library(colourpicker)
  library(bio3d)
  library(ggplot2)
  library(ggrepel)
  library(plotly)
  library(DT)
  library(dplyr)
  library(tidyr)
  library(scales)
  library(fmsb)
})

# ── Molecular analysis helpers ────────────────────────────────────────────────

#' Read a PDB file and return a bio3d pdb object (or NULL on error)
read_structure <- function(path) {
  tryCatch(read.pdb(path, verbose = FALSE), error = function(e) NULL)
}

#' Extract HETATM / ligand records from a pdb object
get_ligand_atoms <- function(pdb) {
  idx <- pdb$atom$type == "HETATM" & pdb$atom$resid != "HOH"
  pdb$atom[idx, ]
}

#' Extract protein CA atoms
get_protein_ca <- function(pdb) {
  idx <- pdb$atom$type == "ATOM" & pdb$atom$elety == "CA"
  pdb$atom[idx, ]
}

#' Euclidean distance between two xyz vectors
atom_dist <- function(a, b) sqrt(sum((a - b)^2))

#' Compute all pairwise distances between two sets of atoms (data.frames with x,y,z)
pairwise_dist <- function(set_a, set_b) {
  coords_a <- as.matrix(set_a[, c("x", "y", "z")])
  coords_b <- as.matrix(set_b[, c("x", "y", "z")])
  # outer product via sweep
  sapply(seq_len(nrow(coords_b)), function(j) {
    sqrt(rowSums((coords_a - matrix(coords_b[j, ], nrow = nrow(coords_a),
                                    ncol = 3, byrow = TRUE))^2))
  })
}

#' Analyse one PDB file and return a named list of metrics
analyse_complex <- function(path, label = NULL,
                             contact_cutoff = 4.0,
                             hbond_cutoff   = 3.5,
                             site_radius    = 5.0) {

  pdb <- read_structure(path)
  if (is.null(pdb)) return(NULL)

  name <- if (!is.null(label) && nzchar(label)) label else
    tools::file_path_sans_ext(basename(path))

  lig  <- get_ligand_atoms(pdb)
  prot <- get_protein_ca(pdb)

  if (nrow(lig) == 0 || nrow(prot) == 0) {
    return(list(
      name = name, filename = basename(path), datapath = path,
      hbond_count = NA, contact_count = NA,
      binding_site_size = NA, binding_site_residues = "",
      ligand_atoms = nrow(lig),
      contact_cutoff = contact_cutoff,
      site_radius = site_radius,
      bfactor_mean = NA, bfactor_sd = NA,
      n_residues = nrow(prot),
      n_chains = length(unique(pdb$atom$chain[pdb$atom$type == "ATOM"]))
    ))
  }

  # Distances between every ligand atom and every protein CA
  dist_mat <- pairwise_dist(lig[, c("x","y","z")], prot[, c("x","y","z")])

  # --- Close contacts ---
  contact_count <- sum(dist_mat <= contact_cutoff)

  # --- Binding-site residues (CA within site_radius of any ligand atom) ---
  min_lig_dist <- apply(dist_mat, 2, min)   # min dist from each CA to any lig atom
  site_idx     <- which(min_lig_dist <= site_radius)
  site_res     <- prot[site_idx, ]
  bs_labels    <- paste0(trimws(site_res$resid), site_res$resno)
  bs_string    <- paste(unique(bs_labels), collapse = ";")

  # --- Approximate H-bond count (N/O atoms within hbond_cutoff) ---
  donor_accept <- c("N", "O", "F")
  lig_no  <- lig[lig$elesy  %in% donor_accept, ]
  prot_no <- pdb$atom[pdb$atom$type == "ATOM" & pdb$atom$elesy %in% donor_accept, ]

  hbond_count <- 0L
  if (nrow(lig_no) > 0 && nrow(prot_no) > 0) {
    hb_mat <- pairwise_dist(lig_no[, c("x","y","z")], prot_no[, c("x","y","z")])
    # Lower bound of 0.5 Å excludes bonded atoms (e.g. covalently bound N-H)
    # that appear in the same residue and are trivially close.
    hbond_count <- sum(hb_mat > 0.5 & hb_mat <= hbond_cutoff)
  }

  # --- B-factors (temperature factors) ---
  bf <- pdb$atom$b[pdb$atom$type == "ATOM"]
  bf <- suppressWarnings(as.numeric(bf))
  bf <- bf[!is.na(bf)]

  list(
    name                 = name,
    filename             = basename(path),
    datapath             = path,          # full temp path; used for atom-table lookup
    hbond_count          = hbond_count,
    contact_count        = contact_count,
    binding_site_size    = length(site_idx),
    binding_site_residues = bs_string,
    ligand_atoms         = nrow(lig),
    contact_cutoff       = contact_cutoff,
    site_radius          = site_radius,
    bfactor_mean         = if (length(bf) > 0) round(mean(bf), 2) else NA,
    bfactor_sd           = if (length(bf) > 0) round(sd(bf), 2) else NA,
    n_residues           = nrow(prot),
    n_chains             = length(unique(pdb$atom$chain[pdb$atom$type == "ATOM"]))
  )
}

# ── Plot templates ────────────────────────────────────────────────────────────

theme_mol <- function() {
  theme_minimal(base_size = 13) +
    theme(
      plot.title    = element_text(face = "bold", hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5, colour = "grey50"),
      legend.position = "bottom"
    )
}

plot_hbond_bar <- function(df) {
  ggplot(df, aes(x = reorder(name, -hbond_count), y = hbond_count,
                 fill = hbond_count)) +
    geom_col(width = 0.65, colour = "white") +
    geom_text(aes(label = hbond_count), vjust = -0.4, size = 3.5, fontface = "bold") +
    scale_fill_gradient(low = "#AED6F1", high = "#1A5276", guide = "none") +
    labs(title = "H-Bond Count per Complex",
         x = "Complex", y = "H-Bond Count") +
    theme_mol()
}

plot_scatter_hb_contacts <- function(df) {
  ggplot(df, aes(x = hbond_count, y = contact_count,
                 size = binding_site_size, colour = name, label = name)) +
    geom_point(alpha = 0.8) +
    geom_text_repel(size = 3, show.legend = FALSE) +
    scale_size_continuous(name = "Binding-site\nresidues", range = c(3, 12)) +
    labs(title = "H-Bonds vs. Close Contacts",
         x = "H-Bond Count", y = "Contact Count", colour = "Complex") +
    theme_mol()
}

plot_bsite_freq <- function(df) {
  res_vec <- unlist(lapply(seq_len(nrow(df)), function(i) {
    r <- trimws(strsplit(df$binding_site_residues[i], ";")[[1]])
    r[nzchar(r)]
  }))
  if (length(res_vec) == 0) return(NULL)
  freq_df <- as.data.frame(sort(table(res_vec), decreasing = TRUE))
  names(freq_df) <- c("residue", "count")
  top_df <- head(freq_df, 20)
  ggplot(top_df, aes(x = reorder(residue, count), y = count)) +
    geom_col(fill = "#1F618D", width = 0.7) +
    coord_flip() +
    labs(title = "Top Binding-Site Residue Frequency",
         x = "Residue", y = "Frequency") +
    theme_mol()
}

plot_radar <- function(df) {
  metrics <- c("hbond_count", "contact_count", "binding_site_size", "ligand_atoms")
  available <- metrics[metrics %in% names(df)]
  if (length(available) < 3) return(NULL)
  radar_df <- df[, c("name", available)]
  radar_df[is.na(radar_df)] <- 0

  # Normalise to 0-1
  for (col in available) {
    rng <- range(radar_df[[col]], na.rm = TRUE)
    if (diff(rng) > 0)
      radar_df[[col]] <- (radar_df[[col]] - rng[1]) / diff(rng)
  }

  mat <- as.data.frame(radar_df[, available])
  rownames(mat) <- radar_df$name

  # fmsb expects max/min rows first
  plot_data <- rbind(rep(1, ncol(mat)), rep(0, ncol(mat)), mat)

  colours <- hue_pal()(nrow(mat))
  op <- par(mar = c(1, 1, 2, 1))
  on.exit(par(op))
  fmsb::radarchart(
    plot_data,
    pcol  = colours,
    pfcol = adjustcolor(colours, 0.25),
    plwd  = 2,
    cglcol = "grey70", axiscol = "grey50",
    vlcex = 0.9,
    title = "Interaction Radar (normalised)"
  )
  legend("topright", legend = rownames(mat), col = colours,
         lty = 1, lwd = 2, cex = 0.8, bty = "n")
}

plot_heatmap <- function(df) {
  rows <- lapply(seq_len(nrow(df)), function(i) {
    r <- trimws(strsplit(df$binding_site_residues[i], ";")[[1]])
    r[nzchar(r)]
  })
  all_res <- unique(unlist(rows))
  if (length(all_res) == 0) return(NULL)

  mat_df <- expand.grid(complex = df$name, residue = all_res,
                        stringsAsFactors = FALSE)
  mat_df$present <- mapply(function(cx, res) {
    idx <- which(df$name == cx)
    res %in% trimws(strsplit(df$binding_site_residues[idx], ";")[[1]])
  }, mat_df$complex, mat_df$residue)

  ggplot(mat_df, aes(x = residue, y = complex, fill = present)) +
    geom_tile(colour = "white", linewidth = 0.5) +
    scale_fill_manual(values = c("TRUE" = "#1A5276", "FALSE" = "#D6EAF8"),
                      labels = c("TRUE" = "Present", "FALSE" = "Absent"),
                      name = "") +
    labs(title = "Binding-Site Residue Fingerprint Heatmap",
         x = "Residue", y = "Complex") +
    theme_mol() +
    theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 8))
}

plot_property_space <- function(df) {
  ggplot(df, aes(x = ligand_atoms, y = hbond_count,
                 colour = contact_count, label = name)) +
    geom_point(size = 4, alpha = 0.85) +
    geom_text_repel(size = 3, show.legend = FALSE) +
    scale_colour_gradient(low = "#AED6F1", high = "#922B21",
                          name = "Contact\nCount") +
    labs(title = "Ligand Property Space",
         x = "Ligand Heavy Atoms", y = "H-Bond Count") +
    theme_mol()
}

plot_bfactor <- function(df) {
  d <- df[!is.na(df$bfactor_mean), ]
  if (nrow(d) == 0) return(NULL)
  ggplot(d, aes(x = reorder(name, -bfactor_mean), y = bfactor_mean,
                ymin = bfactor_mean - bfactor_sd,
                ymax = bfactor_mean + bfactor_sd)) +
    geom_col(fill = "#1F618D", width = 0.65, alpha = 0.85) +
    geom_errorbar(width = 0.25, colour = "#922B21", linewidth = 0.8) +
    labs(title = "Mean B-Factor ± SD per Complex",
         x = "Complex", y = "B-Factor (Å²)") +
    theme_mol()
}

# ── UI ────────────────────────────────────────────────────────────────────────

ui <- dashboardPage(
  skin = "blue",
  dashboardHeader(title = "Molecular Viewer & Analytics"),

  dashboardSidebar(
    sidebarMenu(
      menuItem("📁 Load Structures",  tabName = "load",     icon = icon("folder-open")),
      menuItem("📊 Analysis & Plots", tabName = "analysis", icon = icon("chart-bar")),
      menuItem("🔬 Structure Detail", tabName = "detail",   icon = icon("microscope")),
      menuItem("⚙ Settings",          tabName = "settings", icon = icon("cog"))
    )
  ),

  dashboardBody(
    tags$head(tags$style(HTML("
      .content-wrapper { background: #f4f6f9; }
      .box { border-radius: 6px; }
      pre { background: #1e1e2e; color: #cdd6f4; padding: 10px;
            border-radius: 4px; font-size: 12px; }
    "))),

    tabItems(

      # ── Tab 1: Load ──────────────────────────────────────────────────────────
      tabItem(tabName = "load",
        fluidRow(
          box(title = "Upload PDB Files", width = 6, status = "primary",
            fileInput("pdb_files", "Select PDB file(s)",
                      multiple = TRUE, accept = ".pdb",
                      buttonLabel = "Browse…", placeholder = "No files selected"),
            textInput("labels", "Complex labels (comma-separated, optional)",
                      placeholder = "Complex A, Complex B, …"),
            numericInput("contact_cutoff", "Contact cutoff (Å)", 4.0, 1, 10, 0.5),
            numericInput("hbond_cutoff",   "H-bond cutoff (Å)", 3.5, 1,  6, 0.5),
            numericInput("site_radius",    "Binding-site radius (Å)", 5.0, 1, 15, 0.5),
            actionBttn("analyse_btn", "🔬 Analyse Structures",
                       style = "gradient", color = "primary", block = TRUE)
          ),
          box(title = "Loaded Complexes", width = 6, status = "info",
            DTOutput("summary_table")
          )
        ),
        fluidRow(
          box(title = "Export", width = 12, status = "success",
            downloadButton("download_csv", "📥 Download CSV"),
            tags$span("  "),
            downloadButton("download_rds", "💾 Save Analysis (.rds)")
          )
        )
      ),

      # ── Tab 2: Analysis ──────────────────────────────────────────────────────
      tabItem(tabName = "analysis",
        fluidRow(
          box(title = "Plot Settings", width = 3, status = "primary",
            selectInput("plot_type", "Plot template",
              choices = list(
                "H-Bond Bar Chart"               = "hbond_bar",
                "H-Bonds vs. Contacts Scatter"   = "scatter",
                "Binding-Site Residue Frequency" = "bsite_freq",
                "Interaction Radar Chart"        = "radar",
                "Fingerprint Heatmap"            = "heatmap",
                "Ligand Property Space"          = "property",
                "B-Factor Comparison"            = "bfactor"
              )),
            selectInput("plot_format", "Export format",
                        choices = c("PNG" = "png", "PDF" = "pdf", "SVG" = "svg")),
            numericInput("plot_dpi", "DPI (PNG only)", 300, 72, 1200, 1),
            numericInput("plot_w",   "Width (cm)",  18, 5, 60, 1),
            numericInput("plot_h",   "Height (cm)", 14, 5, 60, 1),
            actionBttn("render_btn", "▶ Render Plot",
                       style = "gradient", color = "success", block = TRUE),
            tags$br(),
            downloadButton("download_plot", "💾 Save Plot")
          ),
          box(title = "Plot Output", width = 9, status = "info",
            plotOutput("main_plot", height = "520px")
          )
        ),
        fluidRow(
          box(title = "Interactive Plot (Plotly)", width = 12, status = "warning",
            plotlyOutput("interactive_plot", height = "400px")
          )
        )
      ),

      # ── Tab 3: Structure Detail ───────────────────────────────────────────────
      tabItem(tabName = "detail",
        fluidRow(
          box(title = "Select Complex", width = 3, status = "primary",
            uiOutput("complex_selector"),
            tags$hr(),
            verbatimTextOutput("structure_summary")
          ),
          box(title = "Binding-Site Residues", width = 9, status = "info",
            DTOutput("bsite_table")
          )
        ),
        fluidRow(
          box(title = "Full Atom Data (first 500 atoms)", width = 12, status = "success",
            DTOutput("atom_table")
          )
        )
      ),

      # ── Tab 4: Settings ───────────────────────────────────────────────────────
      tabItem(tabName = "settings",
        fluidRow(
          box(title = "Application Settings", width = 6, status = "primary",
            h4("Plot Aesthetics"),
            colourpicker::colourInput("bg_colour", "Plot background", "#FFFFFF"),
            selectInput("base_font_size", "Base font size",
                        choices = c("10", "12", "13", "14", "16"), selected = "13"),
            tags$hr(),
            h4("Analysis Defaults"),
            numericInput("def_contact", "Default contact cutoff (Å)", 4.0, 1, 10, 0.5),
            numericInput("def_hbond",   "Default H-bond cutoff (Å)",  3.5, 1,  6, 0.5),
            numericInput("def_radius",  "Default binding-site radius (Å)", 5.0, 1, 15, 0.5),
            tags$hr(),
            h4("Package Information"),
            verbatimTextOutput("pkg_info")
          ),
          box(title = "About", width = 6, status = "info",
            tags$h4("Molecular Viewer & Analytics"),
            tags$p("A pure-R / Shiny dashboard for protein–ligand docking analysis."),
            tags$ul(
              tags$li(tags$b("bio3d"), " — structure I/O and distance calculations"),
              tags$li(tags$b("ggplot2"), " — publication-quality static plots"),
              tags$li(tags$b("plotly"), " — interactive visualisation"),
              tags$li(tags$b("fmsb"), " — radar / spider charts"),
              tags$li(tags$b("shinydashboard"), " — dashboard layout")
            ),
            tags$hr(),
            tags$p("Supported file format: ", tags$code(".pdb")),
            tags$p(tags$b("Analyses performed:")),
            tags$ul(
              tags$li("H-bond count (N/O atom proximity)"),
              tags$li("Close-contact count"),
              tags$li("Binding-site residue identification"),
              tags$li("Ligand heavy-atom count"),
              tags$li("B-factor statistics")
            )
          )
        )
      )
    )
  )
)

# ── Server ────────────────────────────────────────────────────────────────────

server <- function(input, output, session) {

  # -- Reactive: analysis results --------------------------------------------
  results <- reactiveVal(NULL)

  observeEvent(input$analyse_btn, {
    req(input$pdb_files)
    files  <- input$pdb_files
    labels <- trimws(strsplit(input$labels, ",")[[1]])

    withProgress(message = "Analysing structures…", value = 0, {
      res_list <- lapply(seq_len(nrow(files)), function(i) {
        incProgress(1 / nrow(files),
                    detail = paste("Processing", files$name[i]))
        lbl <- if (i <= length(labels) && nzchar(labels[i])) labels[i] else NULL
        analyse_complex(
          path           = files$datapath[i],
          label          = lbl,
          contact_cutoff = input$contact_cutoff,
          hbond_cutoff   = input$hbond_cutoff,
          site_radius    = input$site_radius
        )
      })
    })

    res_list <- Filter(Negate(is.null), res_list)
    if (length(res_list) == 0) {
      showNotification("No valid PDB files found.", type = "error")
      return()
    }

    df <- bind_rows(lapply(res_list, as.data.frame, stringsAsFactors = FALSE))
    results(df)
    showNotification(paste(nrow(df), "complex(es) analysed."), type = "message")
  })

  # -- Summary table ----------------------------------------------------------
  output$summary_table <- renderDT({
    req(results())
    cols <- c("name","filename","hbond_count","contact_count",
              "binding_site_size","ligand_atoms","bfactor_mean","n_residues","n_chains")
    show <- intersect(cols, names(results()))
    datatable(results()[, show], rownames = FALSE,
              options = list(pageLength = 10, scrollX = TRUE),
              caption = "Summary statistics per complex")
  })

  # -- Reactive plot ----------------------------------------------------------
  current_plot <- reactiveVal(NULL)

  make_plot <- reactive({
    req(results(), input$plot_type)
    df <- results()
    pt <- input$plot_type
    if (pt == "hbond_bar")  return(plot_hbond_bar(df))
    if (pt == "scatter")    return(plot_scatter_hb_contacts(df))
    if (pt == "bsite_freq") return(plot_bsite_freq(df))
    if (pt == "property")   return(plot_property_space(df))
    if (pt == "bfactor")    return(plot_bfactor(df))
    if (pt == "heatmap")    return(plot_heatmap(df))
    NULL  # radar handled separately (base graphics)
  })

  observeEvent(input$render_btn, {
    req(results())
    if (input$plot_type == "radar") {
      current_plot(NULL)   # signal base-graphic path
    } else {
      p <- make_plot()
      current_plot(p)
    }
  })

  output$main_plot <- renderPlot({
    req(input$render_btn > 0, results())
    if (input$plot_type == "radar") {
      plot_radar(results())
    } else {
      req(current_plot())
      print(current_plot())
    }
  }, bg = "white")

  # -- Interactive (plotly) ---------------------------------------------------
  output$interactive_plot <- renderPlotly({
    req(results())
    df <- results()
    if (input$plot_type %in% c("hbond_bar","scatter","property","bfactor","heatmap")) {
      p <- make_plot()
      if (!is.null(p)) return(ggplotly(p, tooltip = "all"))
    }
    # default: H-bond bar
    ggplotly(plot_hbond_bar(df), tooltip = "all")
  })

  # -- Downloads --------------------------------------------------------------
  output$download_csv <- downloadHandler(
    filename = function() paste0("molecular_analysis_", Sys.Date(), ".csv"),
    content  = function(file) {
      req(results())
      write.csv(results(), file, row.names = FALSE)
    }
  )

  output$download_rds <- downloadHandler(
    filename = function() paste0("molecular_analysis_", Sys.Date(), ".rds"),
    content  = function(file) {
      req(results())
      saveRDS(results(), file)
    }
  )

  output$download_plot <- downloadHandler(
    filename = function() paste0("plot_", Sys.Date(), ".", input$plot_format),
    content  = function(file) {
      req(results())
      fmt  <- input$plot_format
      w_in <- input$plot_w / 2.54
      h_in <- input$plot_h / 2.54

      if (input$plot_type == "radar") {
        grDevices::png(file, width = w_in, height = h_in, units = "in",
                       res = input$plot_dpi)
        plot_radar(results())
        dev.off()
        return()
      }

      p <- make_plot()
      req(p)
      if (fmt == "png") {
        ggplot2::ggsave(file, plot = p, width = w_in, height = h_in,
                        units = "in", dpi = input$plot_dpi)
      } else {
        ggplot2::ggsave(file, plot = p, width = w_in, height = h_in,
                        units = "in", device = fmt)
      }
    }
  )

  # -- Structure Detail -------------------------------------------------------
  output$complex_selector <- renderUI({
    req(results())
    selectInput("selected_complex", "Complex", choices = results()$name)
  })

  output$structure_summary <- renderText({
    req(results(), input$selected_complex)
    row <- results()[results()$name == input$selected_complex, ]
    paste0(
      "Complex : ", row$name, "\n",
      "File    : ", row$filename, "\n",
      "Residues: ", row$n_residues, "\n",
      "Chains  : ", row$n_chains, "\n",
      "Lig atoms: ", row$ligand_atoms, "\n",
      "H-bonds : ", row$hbond_count, "\n",
      "Contacts: ", row$contact_count, "\n",
      "B-factor: ", row$bfactor_mean, " ± ", row$bfactor_sd, " Å²"
    )
  })

  output$bsite_table <- renderDT({
    req(results(), input$selected_complex)
    row  <- results()[results()$name == input$selected_complex, ]
    res  <- trimws(strsplit(row$binding_site_residues, ";")[[1]])
    res  <- res[nzchar(res)]
    data.frame(Residue = res)
  })

  output$atom_table <- renderDT({
    req(input$selected_complex, results())
    row   <- results()[results()$name == input$selected_complex, ]
    fpath <- row$datapath[1]
    pdb   <- read_structure(fpath)
    req(!is.null(pdb))
    head(pdb$atom, 500)
  }, options = list(scrollX = TRUE, pageLength = 15))

  # -- Settings ---------------------------------------------------------------
  output$pkg_info <- renderText({
    pkgs <- c("shiny","bio3d","ggplot2","plotly","dplyr","fmsb")
    paste(sapply(pkgs, function(p) {
      v <- tryCatch(as.character(packageVersion(p)), error = function(e) "not installed")
      paste0(p, " ", v)
    }), collapse = "\n")
  })
}

# ── Launch ────────────────────────────────────────────────────────────────────
shinyApp(ui = ui, server = server)
