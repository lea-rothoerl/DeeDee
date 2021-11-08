#' DeeDee App
#'
#' @description `deedee_app` opens the DeeDee Shiny web application, combining
#' the functionalities of all other DeeDee functions with a user-friendly
#' graphical user interface.
#'
#' @param deedee_obj A list of DeeDee tables to be analyzed.
#'
#' @return A shiny app
#' @export
#'
#' @examples
#'
#' data(DE_results_IFNg_naive, package = "DeeDee")
#' IFNg_naive <- deedee_prepare(IFNg_naive, "DESeq2")
#'
#' data(DE_results_IFNg_both, package = "DeeDee")
#' IFNg_both <- deedee_prepare(IFNg_both, "DESeq2")
#'
#' data(DE_results_Salm_naive, package = "DeeDee")
#' Salm_naive <- deedee_prepare(Salm_naive, "DESeq2")
#'
#' data(DE_results_Salm_both, package = "DeeDee")
#' Salm_both <- deedee_prepare(Salm_both, "DESeq2")
#'
#' dd_list <- list(
#'   IFNg_naive = IFNg_naive, IFNg_both = IFNg_both,
#'   Salm_naive = Salm_naive, Salm_both = Salm_both
#' )
#'
#' if (interactive()) {
#'   deedee_app(dd_list)
#' }
deedee_app <- function(deedee_obj = NULL) {

  # ----------------------------------------------------------------------------
  # --------------------------------- U I --------------------------------------
  # ----------------------------------------------------------------------------

  deedee_ui <- shiny::navbarPage("DeeDee",
    theme = shinythemes::shinytheme("flatly"),


    # ----------------------------- data input ---------------------------------
    shiny::tabPanel(
      "Input",
      shiny::fluidRow(
        shiny::column(
          8,
          shiny::fileInput("inp", "Upload your DEA results or DeeDee objects",
            multiple = TRUE,
            accept = c(".rds", ".txt", ".xlsx")
          ),
          shiny::tableOutput("inp_infobox")
        ),
        shiny::column(
          4,
          shiny::selectInput("organism", "Organism",
            choices = list(
              "Human" = "org.Hs.eg.db",
              "Mouse" = "org.Mm.eg.db",
              "Fly" = "org.Dm.eg.db",
              "Rat" = "org.Rn.eg.db"
            )
          ),
          shiny::uiOutput("key_inp"),
          shiny::uiOutput("datasets"),
          shiny::conditionalPanel(
            "output.inp_infobox",
            shiny::downloadButton(
              "inp_download",
              "Download DeeDee object (.RDS)"
            )
          )
        )
      ),
      shinyBS::bsCollapse(
        shinyBS::bsCollapsePanel("INFO",
          shiny::includeMarkdown(system.file("extdata",
            "input.md",
            package = "DeeDee"
          )),
          style = "primary"
        )
      ),

      # shiny::downloadButton("vignette",
      #                "Download DeeDee Package vignette (.html)")
    ),


    # ------------------------------- scatter ----------------------------------
    shiny::tabPanel(
      "Scatterplot",
      shiny::fluidRow(
        shiny::column(
          4,
          shiny::uiOutput("scatter_choices1"),
          shiny::uiOutput("scatter_choices2"),
          shiny::selectInput("scatter_color_by", "Color by",
            choices = list(
              "1st p-value" = "pval1",
              "2nd p-value" = "pval2"
            ),
            selected = "pval1"
          ),
          shiny::numericInput("scatter_pthresh", "P-value threshold",
            value = 0.05, min = 0.01, max = 1, step = 0.01
          ),
          shiny::actionButton("ora_button", "Over-representation analysis")
        ),
        shiny::column(
          8,
          shinycssloaders::withSpinner(
            shiny::plotOutput("scatter",
              # dblclick = "scatter_dblclick",
              brush = shiny::brushOpts(
                id = "scatter_brush",
                resetOnNew = FALSE
              )
            )
          )
        )
      ),
      shinyBS::bsCollapse(
        shinyBS::bsCollapsePanel("INFO",
          shiny::includeMarkdown(system.file("extdata",
            "scatter.md",
            package = "DeeDee"
          )),
          style = "primary"
        )
      ),
      shinyBS::bsModal(
        id = "modalExample",
        title = "Gene Ontology over-representation analysis",
        trigger = "ora_button",
        size = "large",
        shinycssloaders::withSpinner(shiny::plotOutput("scatter_ora")),
        shiny::downloadButton(
          "ora_download",
          "Download enrichment result object (.RDS)"
        )
      ),
      shiny::downloadButton(
        "scatter_brush_download",
        "Download brushed genes (.xlsx)"
      ),
      shiny::tableOutput("scatter_brush_info")
    ),


    # ------------------------------- heatmap ----------------------------------
    shiny::tabPanel(
      "Heatmap",
      shiny::fluidRow(
        shiny::column(
          4,
          shiny::numericInput("heatmap_show_first",
            "Show first",
            value = 25,
            min = 1
          ),
          shiny::checkboxInput("heatmap_show_gene_names",
            "Show gene names",
            value = FALSE
          ),
          shiny::checkboxInput("heatmap_showNA",
            "Show NA",
            value = FALSE
          ),
          shiny::selectInput("heatmap_dist", "Distance measure",
            choices = list(
              "Euclidean" = "euclidean",
              "Manhattan" = "manhattan",
              "Pearson" = "pearson",
              "Spearman" = "spearman"
            ),
            selected = "euclidean"
          ),
          shiny::selectInput("heatmap_clust", "Clustering method",
            choices = list(
              "Single" = "single",
              "Complete" = "complete",
              "Average" = "average",
              "Centroid" = "centroid"
            ),
            selected = "average"
          ),
          shiny::numericInput("heatmap_pthresh", "P-value threshold",
            value = 0.05, min = 0.01, max = 1, step = 0.01
          )
        ),
        shiny::column(
          8,
          shiny::textOutput("heatmap_errors"),
          shiny::conditionalPanel(
            "output.heatmap_errors == ''",
            shinycssloaders::withSpinner(
              InteractiveComplexHeatmap::
              InteractiveComplexHeatmapOutput()
            )
          )
        )
      ),
      shinyBS::bsCollapse(
        shinyBS::bsCollapsePanel("INFO",
          shiny::includeMarkdown(system.file("extdata",
            "heatmap.md",
            package = "DeeDee"
          )),
          style = "primary"
        )
      )
    ),


    # -------------------------------- venn ------------------------------------
    shiny::tabPanel(
      "Venn Diagram",
      shiny::fluidRow(
        shiny::column(
          4,
          shiny::selectInput("venn_mode", "Mode",
            choices = list(
              "Up" = "up",
              "Down" = "down",
              "Both" = "both"
            ),
            selected = "both"
          ),
          shiny::numericInput("venn_pthresh", "P-value threshold",
            value = 0.05, min = 0.01, max = 1, step = 0.01
          )
        ),
        shiny::column(
          8,
          shinycssloaders::withSpinner(
            shiny::plotOutput("venn")
          )
        )
      ),
      shinyBS::bsCollapse(
        shinyBS::bsCollapsePanel("INFO",
          shiny::includeMarkdown(system.file("extdata",
            "venn.md",
            package = "DeeDee"
          )),
          style = "primary"
        )
      )
    ),


    # -------------------------------- upSet -----------------------------------
    shiny::tabPanel(
      "UpSet Plot",
      shiny::fluidRow(
        shiny::column(
          4,
          shiny::selectInput("upset_mode", "Mode",
            choices = list(
              "Up" = "up",
              "Down" = "down",
              "Both" = "both"
            ),
            selected = "both"
          ),
          shiny::conditionalPanel(
            condition = "input.upset_mode == 'both'",
            shiny::checkboxInput(
              "upset_colored",
              "Coloring",
              TRUE
            )
          ),
          shiny::numericInput("upset_minset", "Minimum set size",
            value = 10, min = 0, step = 1
          ),
          shiny::numericInput("upset_pthresh", "P-value threshold",
            value = 0.05, min = 0.01, max = 1, step = 0.01
          )
        ),
        shiny::column(
          8,
          shinycssloaders::withSpinner(
            shiny::plotOutput("upset")
          )
        )
      ),
      shinyBS::bsCollapse(
        shinyBS::bsCollapsePanel("INFO",
          shiny::includeMarkdown(system.file("extdata",
            "upset.md",
            package = "DeeDee"
          )),
          style = "primary"
        )
      )
    ),


    # --------------------------------- qq -------------------------------------
    shiny::tabPanel(
      "Quantile-Quantile Plot",
      shiny::fluidRow(
        shiny::column(
          4,
          shiny::checkboxInput("qq_multiple", "Multiple", value = FALSE),
          shiny::conditionalPanel(
            condition = "!input.qq_multiple",
            shiny::uiOutput("qq_choices1"),
            shiny::uiOutput("qq_choices2"),
            shiny::selectInput("qq_color_by", "Color by",
              choices = list(
                "1st p-value" = "pval1",
                "2nd p-value" = "pval2"
              ),
              selected = "pval1"
            ),
            shiny::checkboxInput("qq_line", "As line", value = FALSE)
          ),
          shiny::conditionalPanel(
            condition = "input.qq_multiple",
            shiny::uiOutput("qq_ref"),
          ),
          shiny::numericInput("qq_pthresh", "P-value threshold",
            value = 0.05, min = 0.01, max = 1, step = 0.01
          )
        ),
        shiny::column(
          8,
          shiny::conditionalPanel(
            condition = "!input.qq_multiple",
            shinycssloaders::withSpinner(
              shiny::plotOutput("qq",
                brush = "qq_brush"
              )
            )
          ),
          shiny::conditionalPanel(
            condition = "input.qq_multiple",
            shinycssloaders::withSpinner(
              shiny::plotOutput("qq_mult")
            )
          )
        )
      ),
      shinyBS::bsCollapse(
        shinyBS::bsCollapsePanel("INFO",
          shiny::includeMarkdown(system.file("extdata",
            "qq.md",
            package = "DeeDee"
          )),
          style = "primary"
        )
      ),
      shiny::conditionalPanel(
        condition = "!input.qq_multiple",
        shiny::downloadButton(
          "qq_brush_download",
          "Download brushed genes (.xlsx)"
        ),
        shiny::tableOutput("qq_brush_info")
      )
    ),


    # --------------------------------- cat ------------------------------------
    shiny::tabPanel(
      "Concordance At the Top Plot",
      shiny::fluidRow(
        shiny::column(
          4,
          shiny::selectInput("cat_mode", "Mode",
            choices = list(
              "Up" = "up",
              "Down" = "down",
              "Both" = "both"
            ),
            selected = "up"
          ),
          shiny::numericInput("cat_maxrank",
            "Max rank",
            value = 1000,
            min = 1
          ),
          shiny::uiOutput("cat_choice"),
          shiny::numericInput("cat_pthresh", "P-value threshold",
            value = 0.05, min = 0.01, max = 1, step = 0.01
          )
        ),
        shiny::column(
          8,
          shinycssloaders::withSpinner(
            shiny::plotOutput("cat")
          )
        )
      ),
      shinyBS::bsCollapse(
        shinyBS::bsCollapsePanel("INFO",
          shiny::includeMarkdown(system.file("extdata",
            "cat.md",
            package = "DeeDee"
          )),
          style = "primary"
        )
      )
    ),

    # ------------------------------- summary ----------------------------------
    shiny::tabPanel(
      "Summary",
      shiny::numericInput("sum_pthresh", "P-value threshold",
        value = 0.05, min = 0.01, max = 1, step = 0.01
      ),
      shiny::actionButton("sum_button", "Create Summary"),
      shinyBS::bsModal(
        id = "sum_modal",
        title = "DeeDee Summary",
        trigger = "sum_button",
        size = "large",
        shiny::downloadButton(
          "sum_download",
          "Download your DeeDee Summary (.html)"
        ),
        shinycssloaders::withSpinner(shiny::uiOutput("show_html_summary"))
      ),

      # shinyBS::bsCollapse(
      #   shinyBS::bsCollapsePanel("INFO",
      #                            shiny::includeMarkdown(system.file("extdata",
      #                                                               "cat.md",
      #                                                               package = "DeeDee"
      #                            )),
      #                            style = "primary"
      #   )
      # ),

      shinydashboard::tabBox(
        title = "",
        width = 12,
        id = "sum_params", height = "250px",
        shiny::tabPanel(
          "Scatterplot",
          shiny::uiOutput("sum_scatter_choices1"),
          shiny::uiOutput("sum_scatter_choices2"),
          shiny::selectInput("sum_scatter_color_by", "Color by",
            choices = list(
              "1st p-value" = "pval1",
              "2nd p-value" = "pval2"
            ),
            selected = "pval1"
          )
        ),
        shiny::tabPanel(
          "Heatmap",
          shiny::numericInput("sum_heatmap_show_first",
            "Show first",
            value = 25,
            min = 1
          ),
          shiny::checkboxInput("sum_heatmap_show_gene_names",
            "Show gene names",
            value = FALSE
          ),
          shiny::checkboxInput("sum_heatmap_show_na",
            "Show NA",
            value = FALSE
          ),
          shiny::selectInput("sum_heatmap_dist", "Distance measure",
            choices = list(
              "Euclidean" = "euclidean",
              "Manhattan" = "manhattan",
              "Pearson" = "pearson",
              "Spearman" = "spearman"
            ),
            selected = "euclidean"
          ),
          shiny::selectInput("sum_heatmap_clust", "Clustering method",
            choices = list(
              "Single" = "single",
              "Complete" = "complete",
              "Average" = "average",
              "Centroid" = "centroid"
            ),
            selected = "average"
          )
        ),
        shiny::tabPanel(
          "Venn Diagram",
          shiny::selectInput("sum_venn_mode", "Mode",
            choices = list(
              "Up" = "up",
              "Down" = "down",
              "Both" = "both"
            ),
            selected = "both"
          )
        ),
        shiny::tabPanel(
          "UpSet Plot",
          shiny::selectInput("sum_upset_mode", "Mode",
            choices = list(
              "Up" = "up",
              "Down" = "down",
              "Both" = "both"
            ),
            selected = "both"
          ),
          shiny::conditionalPanel(
            condition = "input.sum_upset_mode == 'both'",
            shiny::checkboxInput(
              "sum_upset_colored",
              "Coloring",
              TRUE
            )
          ),
          shiny::numericInput("sum_upset_min_setsize", "Minimum set size",
            value = 10, min = 0, step = 1
          )
        ),
        shiny::tabPanel(
          "Quantile-Quantile Plot",
          shiny::uiOutput("sum_qq_ref"),
        ),
        shiny::tabPanel(
          "Concordance At the Top Plot",
          shiny::selectInput("sum_cat_mode", "Mode",
            choices = list(
              "Up" = "up",
              "Down" = "down",
              "Both" = "both"
            ),
            selected = "up"
          ),
          shiny::numericInput("sum_cat_maxrank",
            "Max rank",
            value = 1000,
            min = 1
          ),
          shiny::uiOutput("sum_cat_choice")
        )
      )
    )
  )



  # ----------------------------------------------------------------------------
  # ----------------------------- S E R V E R ----------------------------------
  # ----------------------------------------------------------------------------

  deedee_server <- function(input, output, session) {

    # ----------------------------- data input ---------------------------------
    output$key_inp <- shiny::renderUI({
      shiny::req(input$organism)
      anno <- input$organism
      require(anno, character.only = TRUE)
      shiny::selectInput("key_type",
        "Key type of gene IDs",
        choices = keytypes(get(anno))
      )
    })

    mydata <- shiny::reactive({
      shiny::req(shiny::isTruthy(input$inp) || shiny::isTruthy(deedee_obj))

      ext <- c()
      res <- list()

      if (!is.null(deedee_obj)) {
        ext[1] <- "arg"
        res[[1]] <- deedee_obj
      }

      k <- 0

      # reading out input files
      if (length(input$inp[, 1] > 0)) {
        for (i in (length(ext) + 1):(length(input$inp[, 1]) + length(ext))) {
          k <- k + 1
          ext[i] <- tools::file_ext(input$inp[k, "datapath"])
          shiny::validate(shiny::need(
            ext[[i]] == "rds" ||
              ext[[i]] == "RDS" ||
              ext[[i]] == "xlsx" ||
              ext[[i]] == "txt",
            "Please upload only .RDS, .xlsx or .txt files"
          ))

          # .RDS input
          if (ext[[i]] == "rds" || ext[[i]] == "RDS") {
            res[[i]] <- readRDS(input$inp[[k, "datapath"]])
            if (class(res[[i]]) == "DESeqResults") {
              res[[i]] <- deedee_prepare(res[[i]], "DESeq2")
              res[[i]] <- list(res[[i]])
              names(res[[i]]) <- unlist(strsplit(input$inp[k, "name"],
                split = ".",
                fixed = TRUE
              ))[1]
            } else if (class(res[[i]]) == "DGEExact") {
              res[[i]] <- deedee_prepare(res[[i]], "edgeR")
              res[[i]] <- list(res[[i]])
              names(res[[i]]) <- unlist(strsplit(input$inp[k, "name"],
                split = ".",
                fixed = TRUE
              ))[1]
            } else if (class(res[[i]]) == "list") {
              for (j in length(res[[i]])) {
                if (checkmate::test_subset(
                  names(res[[i]][[j]]),
                  c("logFC", "pval")
                ) == FALSE) {
                  return(NULL)
                }
              }
            } else if (class(res[[i]]) == "data.frame") {
              if (length(res[[i]]) == 2) {
                if (checkmate::test_subset(
                  names(res[[i]]),
                  c("logFC", "pval")
                ) == FALSE) {
                  return(NULL)
                }
                res[[i]] <- list(res[[i]])
                names(res[[i]]) <- unlist(strsplit(input$inp[k, "name"],
                  split = ".",
                  fixed = TRUE
                ))[1]
              } else if (length(res[[i]]) == 6) {
                if (checkmate::test_subset(names(res[[i]]), c(
                  "logFC",
                  "AveExpr",
                  "t",
                  "P.Value",
                  "adj.P.Val",
                  "B"
                )) == FALSE) {
                  return(NULL)
                }
                res[[i]] <- deedee_prepare(res[[i]], "limma")
                res[[i]] <- list(res[[i]])
                names(res[[i]]) <- unlist(strsplit(input$inp[k, "name"],
                  split = ".",
                  fixed = TRUE
                ))[1]
              } else {
                return(NULL)
              }
            }
            # .xlsx input
          } else if (ext[[i]] == "xlsx") {
            sheets <- readxl::excel_sheets(input$inp[[k, "datapath"]])
            if (length(sheets) > 1) {
              res[[i]] <- lapply(sheets,
                readxl::read_excel,
                path = input$inp[[k, "datapath"]]
              )
              names(res[[i]]) <- sheets
              for (j in 1:length(sheets)) {
                res[[i]][[sheets[j]]] <- as.data.frame(res[[i]][[sheets[j]]])
                res[[i]][[sheets[j]]] <- tibble::column_to_rownames(
                  res[[i]][[sheets[j]]], "rowname"
                )
                if (checkmate::test_subset(
                  names(res[[i]][[j]]) == FALSE,
                  c("logFC", "pval")
                )) {
                  return(NULL)
                }
              }
            } else {
              res[[i]] <- readxl::read_excel(path = input$inp[[k, "datapath"]])
              res[[i]] <- tibble::column_to_rownames(res[[i]], "rowname")
              res[[i]] <- list(res[[i]])
              names(res[[i]]) <- unlist(strsplit(input$inp[k, "name"],
                split = ".",
                fixed = TRUE
              ))[1]
            }
            # .txt input
          } else if (ext[[i]] == "txt") {
            temp <- utils::read.table(input$inp[[k, "datapath"]])
            res[[i]] <- list(temp)
            names(res[[i]]) <- unlist(strsplit(input$inp[k, "name"],
              split = ".",
              fixed = TRUE
            ))[1]
          } else {
            return(NULL)
          }
        }
      }

      # merging input data structures
      dat <- list()
      for (i in 1:length(res)) {
        for (j in 1:length(res[[i]])) {
          checkmate::assert_subset(names(res[[i]][[j]]), c("logFC", "pval"))
          dat[[names(res[[i]])[[j]]]] <- res[[i]][[j]]
        }
      }

      return(dat)
    })

    output$datasets <- shiny::renderUI({
      shiny::req(shiny::isTruthy(input$inp) || shiny::isTruthy(deedee_obj))
      shiny::checkboxGroupInput("select_datasets",
        "Select datasets to be used",
        choices = names(mydata()),
        selected = names(mydata())
      )
    })

    mydata_use <- shiny::reactive({
      shiny::req(input$inp)
      use <- input$select_datasets
      dat2 <- list()
      for (i in use) {
        dat2[i] <- mydata()[i]
      }
      return(dat2)
    })

    output$inp_download <- shiny::downloadHandler(
      filename = "DeeDee_object.RDS",
      content = function(file) {
        shiny::req(input$inp)
        saveRDS(mydata_use(), file)
      }
    )

    output$inp_infobox <- shiny::renderTable({
      shiny::validate(shiny::need(
        !is.null(mydata()),
        "Faulty input data provided."
      ))

      shiny::req(shiny::isTruthy(input$inp) || shiny::isTruthy(deedee_obj))

      ext <- c()
      filename <- c()
      res <- list()
      type <- c()
      contrast <- c()
      genes <- c()
      count <- 0

      if (!is.null(deedee_obj)) {
        for (j in 1:length(deedee_obj)) {
          count <- count + 1
          type[count] <- "DeeDee object"
          filename[count] <- "input as argument"
          contrast[count] <- names(deedee_obj)[j]
          genes[count] <- length(deedee_obj[j]
          [[contrast[count]]][["logFC"]])
        }
      }

      # reading out input files
      if (length(input$inp[, 1] > 0)) {
        for (i in 1:(length(input$inp[, 1]))) {
          ext[i] <- tools::file_ext(input$inp[i, "datapath"])

          if (ext[[i]] == "rds" || ext[[i]] == "RDS") {
            res[[i]] <- readRDS(input$inp[[i, "datapath"]])
          } else if (ext[[i]] == "xlsx") {
            sheets <- readxl::excel_sheets(input$inp[[i, "datapath"]])
            res[[i]] <- lapply(sheets,
              readxl::read_excel,
              path = input$inp[[i, "datapath"]]
            )
            names(res[[i]]) <- sheets
            for (j in 1:length(sheets)) {
              res[[i]][[sheets[j]]] <- as.data.frame(res[[i]][[sheets[j]]])
              res[[i]][[sheets[j]]] <- tibble::column_to_rownames(
                res[[i]][[sheets[j]]], "rowname"
              )
            }
          } else if (ext[[i]] == "txt") {
            temp <- utils::read.table(input$inp[[i, "datapath"]])
            res[[i]] <- list(temp)
            names(res[[i]]) <- unlist(strsplit(input$inp[i, "name"],
              split = ".",
              fixed = TRUE
            ))[1]
          }

          if (class(res[[i]]) == "DESeqResults" ||
            class(res[[i]]) == "DGEExact" ||
            length(names(res[[i]])) == 6) {
            count <- count + 1
            type[count] <- class(res[[i]])
            filename[count] <- input$inp[i, "name"]
            contrast[count] <- unlist(strsplit(filename[count],
              split = ".",
              fixed = TRUE
            ))[1]
            genes[count] <- length(mydata()
            [[contrast[count]]][["logFC"]])
          } else {
            if (class(res[[i]]) == "data.frame") {
              res[[i]] <- list(res[[i]])
              names(res[[i]]) <- unlist(strsplit(input$inp[i, "name"],
                split = ".",
                fixed = TRUE
              ))[1]
            }
            for (j in 1:length(res[[i]])) {
              count <- count + 1
              type[count] <- "DeeDee object"
              filename[count] <- input$inp[i, "name"]
              contrast[count] <- names(res[[i]])[j]
              genes[count] <- length(res[[i]][j]
              [[contrast[count]]][["logFC"]])
            }
          }
        }
      }

      df <- data.frame(filename, type, contrast, genes)
      return(df)
    })


    # ------------------------------- scatter ----------------------------------
    # --- selectors ---
    output$scatter_choices1 <- shiny::renderUI({
      shiny::req(input$inp)
      shiny::selectInput("scatter_select1",
        "1st data set",
        choices = names(mydata_use())
      )
    })

    output$scatter_choices2 <- shiny::renderUI({
      shiny::req(input$inp)
      shiny::selectInput("scatter_select2",
        "2nd data set",
        selected = names(mydata_use())[2],
        choices = names(mydata_use())
      )
    })

    # --- plot output ---
    ranges <- shiny::reactiveValues(x = NULL, y = NULL)

    output$scatter <- shiny::renderPlot({
      shiny::validate(
        shiny::need(
          input$inp,
          "Please input at least two contrasts."
        )
      )
      shiny::validate(
        shiny::need(
          length(mydata()) >= 2,
          "Please upload at least two contrasts."
        )
      )
      shiny::validate(
        shiny::need(
          length(mydata_use()) >= 2,
          "Please select at least two contrasts."
        )
      )

      sel1 <- match(input$scatter_select1, names(mydata_use()))
      sel2 <- match(input$scatter_select2, names(mydata_use()))
      shiny::req(sel1)
      shiny::req(sel2)
      res <- deedee_scatter(mydata_use(),
        select1 = sel1,
        select2 = sel2,
        color_by = input$scatter_color_by,
        pthresh = input$scatter_pthresh
      )
      shiny::validate(
        shiny::need(!is.null(res), "No common genes in input datasets.")
      )
      res +
        ggplot2::coord_cartesian(
          xlim = ranges$x,
          ylim = ranges$y,
          expand = FALSE
        )
    })

    # --- brushing ---
    # scatter_download_button <- shiny::renderUI(
    #     conditionalPanel("output.scatter_brushed()",
    #     downloadButton("scatter_brush_download",
    #                    "Download brushed genes (.txt)")))

    scatter_brushed <- shiny::reactive({
      shiny::req(input$scatter_brush)
      df <- data.frame(
        x = mydata_use()[[input$scatter_select1]],
        y = mydata_use()[[input$scatter_select2]]
      )

      df <- subset(df, x.pval < 0.05 & y.pval < 0.05)

      names(df) <- c(
        paste(input$scatter_select1, ".logFC", sep = ""),
        paste(input$scatter_select1, ".pval", sep = ""),
        paste(input$scatter_select2, ".logFC", sep = ""),
        paste(input$scatter_select2, ".pval", sep = "")
      )
      shiny::brushedPoints(df,
        input$scatter_brush,
        xvar = paste(input$scatter_select1, ".logFC", sep = ""),
        yvar = paste(input$scatter_select2, ".logFC", sep = "")
      )
    })

    output$scatter_brush_info <- shiny::renderTable(
      {
        shiny::req(scatter_brushed())
        scatter_brushed()
      },
      rownames = TRUE
    )

    output$scatter_brush_download <- shiny::downloadHandler(
      filename = "scatter_brushed_genes.xlsx",
      content = function(file) {
        shiny::req(input$scatter_brush)
        bru <- tibble::rownames_to_column(scatter_brushed())
        first <- data.frame(bru[1], bru[2], bru[3])
        nm1 <- unlist(strsplit(names(first)[2],
          split = ".",
          fixed = TRUE
        ))[1]
        names(first) <- c("rowname", "logFC", "pval")
        second <- data.frame(bru[1], bru[4], bru[5])
        nm2 <- unlist(strsplit(names(second)[2],
          split = ".",
          fixed = TRUE
        ))[1]
        names(second) <- c("rowname", "logFC", "pval")
        l <- list(first, second)
        names(l) <- c(nm1, nm2)
        writexl::write_xlsx(l, file, col_names = TRUE)
      }
    )

    # shiny::observeEvent(input$scatter_dblclick, {
    #     brush <- input$scatter_brush
    #     if (!is.null(brush)) {
    #         ranges$x <- c(brush$xmin, brush$xmax)
    #         ranges$y <- c(brush$ymin, brush$ymax)
    #
    #     } else {
    #         ranges$x <- NULL
    #         ranges$y <- NULL
    #     }
    # })

    # --- enrich ---
    enrich <- shiny::reactive({
      shiny::req(input$inp)
      shiny::validate(shiny::need(
        scatter_brushed(),
        "No brushed genes."
      ))
      sel1 <- match(input$scatter_select1, names(mydata_use()))
      sel2 <- match(input$scatter_select2, names(mydata_use()))
      shiny::req(sel1)
      shiny::req(sel2)
      data <- list(mydata_use()[[sel1]], mydata_use()[[sel2]])
      res <- ora(
        geneList = scatter_brushed(),
        universe = data,
        orgDB = input$organism,
        key_type = input$key_type
      )
      shiny::validate(shiny::need(
        class(res) == "enrichResult",
        "Not working."
      ))
      # temp <- as.data.frame(res)
      # if (dim(temp)[1] == 0) {
      #     return(NULL)
      # }
      return(res)
    })

    output$scatter_ora <- shiny::renderPlot({
      shiny::validate(shiny::need(
        !is.null(enrich()),
        "Something went wrong..."
      ))
      en <- enrich()
      en_df <- as.data.frame(en)
      shiny::validate(shiny::need(
        nrow(en_df) > 0,
        "No enriched terms found."
      ))
      options(ggrepel.max.overlaps = Inf)
      plt <- enrichplot::emapplot(enrichplot::pairwise_termsim(en))
      shiny::validate(shiny::need(
        !is.null(plt),
        "No enriched terms found."
      ))
      print(plt)
    })

    output$ora_download <- shiny::downloadHandler(
      filename = "enrichment_results.RDS",
      content = function(file) {
        shiny::req(!is.null(enrich()))
        saveRDS(enrich(), file)
      }
    )

    # ------------------------------- heatmap ----------------------------------
    output$heatmap_errors <- shiny::renderText({
      shiny::validate(
        shiny::need(
          input$inp,
          "Please input at least two contrasts."
        )
      )
      shiny::validate(
        shiny::need(
          length(mydata()) >= 2,
          "Please upload at least two contrasts."
        )
      )
      shiny::validate(
        shiny::need(
          length(mydata_use()) >= 2,
          "Please select at least two contrasts."
        )
      )
      shiny::validate(
        shiny::need(
          !is.null(heatmap_output()),
          "No common genes in input datasets."
        )
      )
      ""
    })

    heatmap_output <- shiny::reactive({
      shiny::req(input$inp)
      shiny::req(input$heatmap_show_first)
      shiny::req(mydata_use())
      res <- deedee_heatmap(mydata_use(),
        show_first = input$heatmap_show_first,
        show_gene_names = input$heatmap_show_gene_names,
        dist = input$heatmap_dist,
        clust = input$heatmap_clust,
        pthresh = input$heatmap_pthresh,
        show_na = input$heatmap_showNA
      )
      shiny::validate(
        shiny::need(!is.null(res), "No common genes in input datasets.")
      )

      res <- ComplexHeatmap::draw(res)

      return(res)
    })

    listen <- shiny::reactive({
      list(
        input$heatmap_show_first,
        input$heatmap_show_gene_names,
        input$heatmap_dist,
        input$heatmap_clust,
        input$heatmap_pthresh,
        input$heatmap_showNA,
        mydata_use()
      )
    })

    shiny::observeEvent(listen(), {
      shiny::req(length(mydata_use()) >= 2)
      InteractiveComplexHeatmap::makeInteractiveComplexHeatmap(
        input,
        output,
        session,
        heatmap_output()
      )
    })


    # -------------------------------- venn ------------------------------------
    output$venn <- shiny::renderPlot({
      shiny::validate(
        shiny::need(
          input$inp,
          "Please input at least two contrasts."
        )
      )
      shiny::validate(
        shiny::need(
          length(mydata()) >= 2,
          "Please upload at least two contrasts."
        )
      )
      shiny::validate(
        shiny::need(
          length(mydata_use()) >= 2,
          "Please select at least two contrasts."
        )
      )

      shiny::req(input$inp)
      res <- deedee_venn(mydata_use(),
        mode = input$venn_mode,
        pthresh = input$venn_pthresh
      )

      shiny::validate(
        shiny::need(
          !is.null(res),
          "No genes in your datasets. Maybe your specified p-value threshold is too low?"
        )
      )

      res
    })


    # -------------------------------- upset -----------------------------------
    output$upset <- shiny::renderPlot({
      shiny::validate(
        shiny::need(
          input$inp,
          "Please input at least two contrasts."
        )
      )
      shiny::validate(
        shiny::need(
          length(mydata()) >= 2,
          "Please upload at least two contrasts."
        )
      )
      shiny::validate(
        shiny::need(
          length(mydata_use()) >= 2,
          "Please select at least two contrasts."
        )
      )

      shiny::req(input$inp)
      if (input$upset_mode == "both" && input$upset_colored) {
        mode <- "both_colored"
      } else {
        mode <- input$upset_mode
      }
      res <- deedee_upset(mydata_use(),
        mode = mode,
        pthresh = input$upset_pthresh,
        min_setsize = input$upset_minset
      )

      shiny::validate(
        shiny::need(
          !is.null(res),
          "No genes in your datasets. Maybe your specified p-value threshold is too low?"
        )
      )
      res
    })


    # --------------------------------- qq -------------------------------------
    output$qq_choices1 <- shiny::renderUI({
      shiny::req(input$inp)
      shiny::selectInput("qq_select1",
        "1st data set",
        choices = names(mydata_use())
      )
    })

    output$qq_choices2 <- shiny::renderUI({
      shiny::req(input$inp)
      shiny::selectInput("qq_select2",
        "2nd data set",
        selected = names(mydata_use())[2],
        choices = names(mydata_use())
      )
    })

    output$qq_ref <- shiny::renderUI({
      shiny::req(input$inp)
      shiny::selectInput("qq_reference",
        "Reference",
        choices = names(mydata_use())
      )
    })

    output$qq <- shiny::renderPlot({
      shiny::validate(
        shiny::need(
          input$inp,
          "Please input at least two contrasts."
        )
      )
      shiny::validate(
        shiny::need(
          length(mydata()) >= 2,
          "Please upload at least two contrasts."
        )
      )
      shiny::validate(
        shiny::need(
          length(mydata_use()) >= 2,
          "Please select at least two contrasts."
        )
      )

      shiny::req(input$inp)
      sel1 <- match(input$qq_select1, names(mydata_use()))
      sel2 <- match(input$qq_select2, names(mydata_use()))
      shiny::req(sel1)
      shiny::req(sel2)
      res <- deedee_qq(mydata_use(),
        select1 = sel1,
        select2 = sel2,
        color_by = input$qq_color_by,
        pthresh = input$qq_pthresh,
        as_line = input$qq_line
      )

      shiny::validate(
        shiny::need(
          !is.null(res),
          "No genes in your datasets. Maybe your specified p-value threshold is too low?"
        )
      )

      res
    })

    output$qq_mult <- shiny::renderPlot({
      shiny::validate(
        shiny::need(
          input$inp,
          "Please input at least two contrasts."
        )
      )
      shiny::validate(
        shiny::need(
          length(mydata()) >= 2,
          "Please upload at least two contrasts."
        )
      )
      shiny::validate(
        shiny::need(
          length(mydata_use()) >= 2,
          "Please select at least two contrasts."
        )
      )

      shiny::req(input$inp)
      ref <- match(input$qq_reference, names(mydata_use()))
      shiny::req(ref)
      res <- deedee_qqmult(mydata_use(),
        ref = ref,
        pthresh = input$qq_pthresh
      )

      shiny::validate(
        shiny::need(
          !is.null(res),
          "No genes in your datasets. Maybe your specified p-value threshold is too low?"
        )
      )

      res
    })

    qq_brushed <- shiny::reactive({
      shiny::req(input$qq_brush)
      x <- mydata()[[input$qq_select1]]$logFC
      y <- mydata()[[input$qq_select2]]$logFC
      pval1 <- mydata()[[input$qq_select2]]$pval
      pval2 <- mydata()[[input$qq_select2]]$pval
      names(x) <- row.names(mydata()[[input$qq_select2]])
      names(y) <- row.names(mydata()[[input$qq_select2]])

      sx_idx <- order(x)
      sy_idx <- order(y)

      sx <- x[sx_idx]
      sy <- y[sy_idx]
      pval1 <- pval1[sx_idx]
      pval2 <- pval2[sy_idx]

      lenx <- length(sx)
      leny <- length(sy)

      if (leny < lenx) {
        sx <- stats::approx(1L:lenx, sx, n = leny)$y
        pval1 <- stats::approx(1L:lenx, pval1, n = leny)$y
      }
      if (leny > lenx) {
        sy <- stats::approx(1L:leny, sy, n = lenx)$y
        pval2 <- stats::approx(1L:leny, pval2, n = lenx)$y
      }

      sx <- tibble::rownames_to_column(as.data.frame(sx))
      sx[3] <- pval1
      sy <- tibble::rownames_to_column(as.data.frame(sy))
      sy[3] <- pval2

      qq <- data.frame(
        x = sx,
        y = sy
      )

      names(qq) <- c(
        paste(input$qq_select1, "gene", sep = "."),
        paste(input$qq_select1, "logFC", sep = "."),
        paste(input$qq_select1, "pval", sep = "."),
        paste(input$qq_select2, "gene", sep = "."),
        paste(input$qq_select2, "logFC", sep = "."),
        paste(input$qq_select2, "pval", sep = ".")
      )

      temp1 <- paste(input$qq_select1, "logFC", sep = ".")
      temp2 <- paste(input$qq_select2, "logFC", sep = ".")

      shiny::brushedPoints(qq,
        input$qq_brush,
        xvar = temp1,
        yvar = temp2
      )
    })

    output$qq_brush_info <- shiny::renderTable(
      {
        shiny::req(qq_brushed())
        qq_brushed()
      },
      rownames = FALSE
    )

    output$qq_brush_download <- shiny::downloadHandler(
      filename = "qq_brushed_genes.xlsx",
      content = function(file) {
        shiny::req(input$qq_brush)
        bru <- qq_brushed()
        first <- data.frame(bru[1], bru[2], bru[3])
        nm1 <- unlist(strsplit(names(first)[2],
          split = ".",
          fixed = TRUE
        ))[1]
        names(first) <- c("rowname", "logFC", "pval")
        second <- data.frame(bru[4], bru[5], bru[6])
        nm2 <- unlist(strsplit(names(second)[2],
          split = ".",
          fixed = TRUE
        ))[1]
        names(second) <- c("rowname", "logFC", "pval")
        l <- list(first, second)
        names(l) <- c(nm1, nm2)
        writexl::write_xlsx(l, file, col_names = TRUE)
      }
    )


    # --------------------------------- cat ------------------------------------
    output$cat_choice <- shiny::renderUI({
      shiny::req(input$inp)
      shiny::selectInput("cat_ref",
        "Reference contrast",
        selected = names(mydata_use())[1],
        choices = names(mydata_use())
      )
    })

    output$cat <- shiny::renderPlot({
      shiny::validate(
        shiny::need(
          input$inp,
          "Please input at least two contrasts."
        )
      )
      shiny::validate(
        shiny::need(
          length(mydata()) >= 2,
          "Please upload at least two contrasts."
        )
      )
      shiny::validate(
        shiny::need(
          length(mydata_use()) >= 2,
          "Please select at least two contrasts."
        )
      )

      shiny::req(input$inp)
      shiny::req(input$cat_maxrank)
      shiny::req(input$cat_ref)

      ref <- match(input$cat_ref, names(mydata_use()))
      res <- deedee_cat(mydata_use(),
        ref = ref,
        maxrank = input$cat_maxrank,
        mode = input$cat_mode,
        pthresh = input$cat_pthresh
      )

      shiny::validate(
        shiny::need(
          !is.null(res),
          "No genes in your datasets. Maybe your specified p-value threshold is too low?"
        )
      )

      res
    })

    # ------------------------------- summary ----------------------------------
    output$sum_scatter_choices1 <- shiny::renderUI({
      shiny::req(input$inp)
      shiny::selectInput("sum_scatter_select1",
        "1st data set",
        choices = names(mydata_use())
      )
    })

    output$sum_scatter_choices2 <- shiny::renderUI({
      shiny::req(input$inp)
      shiny::selectInput("sum_scatter_select2",
        "2nd data set",
        selected = names(mydata_use())[2],
        choices = names(mydata_use())
      )
    })

    output$sum_qq_ref <- shiny::renderUI({
      shiny::req(input$inp)
      shiny::selectInput("sum_qqmult_ref",
        "Reference",
        selected = names(mydata_use())[1],
        choices = names(mydata_use())
      )
    })

    output$sum_cat_choice <- shiny::renderUI({
      shiny::req(input$inp)
      shiny::selectInput("sum_cat_ref",
        "Reference contrast",
        selected = names(mydata_use())[1],
        choices = names(mydata_use())
      )
    })

    summary <- reactive({
      shiny::req(mydata_use())

      outfile <- tempfile(fileext = ".html")

      sc_sel1 <- match(input$sum_scatter_select1, names(mydata_use()))
      if (is.null(sc_sel1)) {
        sc_sel1 <- 1
      }
      sc_sel2 <- match(input$sum_scatter_select2, names(mydata_use()))
      if (is.null(sc_sel2)) {
        sc_sel2 <- 2
      }

      qq_ref <- match(input$sum_qqmult_ref, names(mydata_use()))
      if (length(qq_ref) == 0) {
        qq_ref <- 1
      }

      cat_ref <- match(input$sum_cat_ref, names(mydata_use()))
      if (length(cat_ref) == 0) {
        cat_ref <- 1
      }

      if (input$sum_upset_mode == "both" && input$sum_upset_colored) {
        ups_mode <- "both_colored"
      } else {
        ups_mode <- input$sum_upset_mode
      }

      shiny::req(sc_sel1)
      shiny::req(sc_sel2)
      shiny::req(qq_ref)
      shiny::req(cat_ref)

      deedee_summary(mydata_use(),
        output_path = outfile,
        overwrite = TRUE,
        pthresh = input$sum_pthresh,
        scatter_select1 = sc_sel1,
        scatter_select2 = sc_sel2,
        scatter_color_by = input$sum_scatter_color_by,
        heatmap_show_first = input$sum_heatmap_show_first,
        heatmap_show_gene_names = input$sum_heatmap_show_gene_names,
        heatmap_dist = input$sum_heatmap_dist,
        heatmap_clust = input$sum_heatmap_clust,
        heatmap_show_na = input$sum_heatmap_show_na,
        venn_mode = input$sum_venn_mode,
        upset_mode = ups_mode,
        upset_min_setsize = input$sum_upset_min_setsize,
        qqmult_ref = qq_ref,
        cat_ref = cat_ref,
        cat_maxrank = input$sum_cat_maxrank,
        cat_mode = input$sum_cat_mode,
        silent = TRUE,
        open_file = FALSE
      )

      return(outfile)
    })

    output$sum_download <- shiny::downloadHandler(
      filename = "DeeDee_Summary.html",
      content = function(file) {
        shiny::req(summary())

        file.copy(summary(), file)
      }
    )

    output$show_html_summary <- shiny::renderUI({
      shiny::req(summary())

      shiny::includeHTML(summary())
    })
  }


  # ------------------------------ run application -------------------------------
  shiny::shinyOptions(deedee_obj = deedee_obj)
  shiny::shinyApp(ui = deedee_ui, server = deedee_server)
}
