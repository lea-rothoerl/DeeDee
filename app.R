# ------------------------------- source packages ------------------------------
library(shiny)
source("~/Development/DeeDee/R/deedee_venn.R")
source("~/Development/DeeDee/R/deedee_upSet.R")
source("~/Development/DeeDee/R/deedee_scatter.R")
source("~/Development/DeeDee/R/deedee_heatmap.R")
source("~/Development/DeeDee/R/deedee_qq.R")
source("~/Development/DeeDee/R/deedee_cat.R")
source("~/Development/DeeDee/R/deedee_prepare.R")
source("~/Development/DeeDee/R/gsea.R")
library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(org.Dm.eg.db)
library(org.Rn.eg.db)
library(shinythemes)
library(shinyBS)

# ------------------------------------------------------------------------------
# --------------------------------- U I ----------------------------------------
# ------------------------------------------------------------------------------

ui <- navbarPage("DeeDee", theme = shinytheme("flatly"),

    # ----------------------------- data input ---------------------------------
    tabPanel("Input",
             fluidRow(
                 column(8,
                    fileInput("inp", "Upload your DEA results or DeeDee objects",
                           multiple = TRUE,
                           accept = c(".rds", ".txt", ".xlsx")),
                        tableOutput("inp_infobox")),

                 column(4,
                    selectInput("organism", "Organism",
                                choices = list("Human" = "org.Hs.eg.db",
                                        "Mouse" = "org.Mm.eg.db",
                                        "Fly" = "org.Dm.eg.db",
                                        "Rat" = "org.Rn.eg.db")),
                        uiOutput("key_inp"),
                        uiOutput("datasets"),
                        conditionalPanel("output.inp_infobox",
                             downloadButton("inp_download",
                                    "Download DeeDee object (.RDS)")))),
             # downloadButton("vignette",
             #                "Download DeeDee Package vignette (.html)")
             ),



    # ------------------------------- scatter ----------------------------------
    tabPanel("Scatterplot",
             fluidRow(
                 column(4,
                     uiOutput("scatter_choices1"),

                     uiOutput("scatter_choices2"),

                     selectInput("scatter_color_by", "Color by",
                                 choices = list("1st p-value" = "pval1",
                                                "2nd p-value" = "pval2"),
                                                selected = "pval1"),

                     numericInput("scatter_pthresh" , "P-value threshold",
                                  value = 0.05, min = 0.01, max = 1, step = 0.01),

                     actionButton("ora_button", "Over-representation analysis")),

             column(8,
                shinycssloaders::withSpinner(
                    plotOutput("scatter",
                       dblclick = "scatter_dblclick",
                       brush = brushOpts(id = "scatter_brush",
                                         resetOnNew = FALSE))))),

             shinyBS::bsCollapse(
                 shinyBS::bsCollapsePanel("INFO",
                                          includeMarkdown("scatter.md"),
                                          style = "primary")),


                shinyBS::bsModal(id = "modalExample",
                                 title = "Gene Ontology over-enrichment analysis",
                                 trigger = "ora_button",
                                 size = "large",
                        shinycssloaders::withSpinner(plotOutput("scatter_ora")),
                        downloadButton('ora_download',
                                       'Download enrichment result object (.RDS)')),

             downloadButton("scatter_brush_download",
                            "Download brushed genes (.txt)"),
             tableOutput("scatter_brush_info")),


    # ------------------------------- heatmap ----------------------------------
    tabPanel("Heatmap",
             fluidRow(
                 column(4,
                     numericInput("heatmap_show_first",
                                  "Show first",
                                  value = 25,
                                  min = 1),

                     checkboxInput("heatmap_show_gene_names",
                                   "Show gene names",
                                   value = FALSE),

                     checkboxInput("heatmap_showNA",
                                   "Show NA",
                                   value = FALSE),

                     selectInput("heatmap_dist", "Distance measure",
                                 choices = list("Euclidean" = "euclidean",
                                                "Manhattan" = "manhattan",
                                                "Pearson" = "pearson",
                                                "Spearman" = "spearman"),
                                 selected = "euclidean"),

                     selectInput("heatmap_clust", "Clustering method",
                                 choices = list("Single" = "single",
                                                "Complete" = "complete",
                                                "Average" = "average",
                                                "Centroid" = "centroid"),
                                 selected = "average"),

                     numericInput("heatmap_pthresh" , "P-value threshold",
                                  value = 0.05, min = 0.01, max = 1, step = 0.01)),

                 column(8,
                        textOutput("heatmap_errors"),
                           conditionalPanel("output.heatmap_errors == ''",
                                            shinycssloaders::withSpinner(
                                                InteractiveComplexHeatmap::
                                                    InteractiveComplexHeatmapOutput())))),

                 shinyBS::bsCollapse(
                     shinyBS::bsCollapsePanel("INFO",
                            includeMarkdown("heatmap.md"),
                            style = "primary"))),


    # -------------------------------- venn ------------------------------------
    tabPanel("Venn Diagram",
             fluidRow(
                 column(4,
                     selectInput("venn_mode", "Mode",
                                 choices = list("Up" = "up",
                                                "Down" = "down",
                                                "Both" = "both"),
                                 selected = "both"),

                     numericInput("venn_pthresh" , "P-value threshold",
                                  value = 0.05, min = 0.01, max = 1, step = 0.01)),

                 column(8,
                     shinycssloaders::withSpinner(
                        plotOutput("venn")))),

                 shinyBS::bsCollapse(
                     shinyBS::bsCollapsePanel("INFO",
                        includeMarkdown("venn.md"),
                        style = "primary"))),


    # -------------------------------- upSet -----------------------------------
    tabPanel("UpSet Plot",
             fluidRow(
                 column(4,
                     selectInput("upSet_mode", "Mode",
                                 choices = list("Up" = "up",
                                                "Down" = "down",
                                                "Both" = "both"),
                                 selected = "both"),

                     conditionalPanel(
                         condition = "input.upSet_mode == 'both'",
                         checkboxInput("upSet_colored",
                                       "Coloring",
                                       TRUE)),

                     numericInput("upSet_pthresh" , "P-value threshold",
                                  value = 0.05, min = 0.01, max = 1, step = 0.01),

                     numericInput("upSet_minset", "Minimum set size",
                                  value = 10, min = 0, step = 1)),

                    column(8,
                        shinycssloaders::withSpinner(
                            plotOutput("upSet")))),

                 shinyBS::bsCollapse(
                     shinyBS::bsCollapsePanel("INFO",
                        includeMarkdown("upSet.md"),
                        style = "primary"))),


    # --------------------------------- qq -------------------------------------
    tabPanel("Quantile-Quantile Plot",
             fluidRow(
                 column(4,
                     uiOutput("qq_choices1"),

                     uiOutput("qq_choices2"),

                     selectInput("qq_color_by", "Color by",
                                 choices = list("1st p-value" = "pval1",
                                                "2nd p-value" = "pval2"),
                                 selected = "pval1"),

                     numericInput("qq_pthresh" , "P-value threshold",
                                  value = 0.05, min = 0.01, max = 1, step = 0.01)),

                 column(8,
                     shinycssloaders::withSpinner(
                         plotOutput("qq",
                                      brush = "qq_brush")))),

             shinyBS::bsCollapse(
                 shinyBS::bsCollapsePanel("INFO",
                                          includeMarkdown("qq.md"),
                                          style = "primary")),

             downloadButton("qq_brush_download",
                            "Download brushed genes (.txt)"),

             tableOutput("qq_brush_info")),


    # --------------------------------- cat ------------------------------------
    tabPanel("Concordance At the Top Plot",
             fluidRow(
                 column(4,
                    numericInput("cat_maxrank",
                                "Max rank",
                                value = 1000,
                                min = 1),

                    numericInput("cat_pthresh" , "P-value threshold",
                                 value = 0.05, min = 0.01, max = 1, step = 0.01),

                    uiOutput("cat_choice")),

                 column(8,
                     shinycssloaders::withSpinner(
                         plotOutput("cat")))),

                shinyBS::bsCollapse(
                    shinyBS::bsCollapsePanel("INFO",
                        includeMarkdown("cat.md"),
                        style = "primary"))))

# ------------------------------------------------------------------------------
# ----------------------------- S E R V E R ------------------------------------
# ------------------------------------------------------------------------------

server <- function(input, output, session) {

    # ----------------------------- data input ---------------------------------
    output$key_inp <- renderUI({
        req(input$organism)
        selectInput("key_type",
                    "Key type of gene IDs",
                    choices = keytypes(get(input$organism)))
    })

    mydata <- reactive({
        req(input$inp)

        ext <- vector(mode = "numeric", length = length(input$inp[,1]))
        res <- list()

        # reading out input files
        for(i in 1:length(input$inp[,1])) {
            ext[i] <- tools::file_ext(input$inp[i, "datapath"])
            validate(need(ext[[i]] == "rds"||
                              ext[[i]] == "RDS" ||
                              ext[[i]] == "xlsx" ||
                              ext[[i]] == "txt",
                          "Please upload only .RDS, .xlsx or .txt files"))

            # .RDS input
            if (ext[[i]] == "rds"|| ext[[i]] == "RDS") {
                res[[i]] <- readRDS(input$inp[[i, "datapath"]])
                if (class(res[[i]]) == "DESeqResults") {
                    res[[i]] <- deedee_prepare(res[[i]], "DESeq2")
                    res[[i]] <- list(res[[i]])
                    names(res[[i]]) <- unlist(strsplit(input$inp[i, "name"],
                                                       split=".",
                                                       fixed=TRUE))[1]
                }
                else if (class(res[[i]]) == "DGEExact") {
                    res[[i]] <- deedee_prepare(res[[i]], "edgeR")
                    res[[i]] <- list(res[[i]])
                    names(res[[i]]) <- unlist(strsplit(input$inp[i, "name"],
                                                       split=".",
                                                       fixed=TRUE))[1]
                }
                else if (class(res[[i]]) == "list") {
                    for (j in length(res[[i]])) {
                        if (checkmate::test_subset(names(res[[i]][[j]]),
                                                c("logFC", "pval")) == FALSE) {
                            return(NULL)
                        }
                    }
                }
                else if (class(res[[i]]) == "data.frame") {
                    if (length(res[[i]]) == 2) {
                        if (checkmate::test_subset(names(res[[i]]),
                                            c("logFC", "pval")) == FALSE) {
                            return(NULL)
                        }
                        res[[i]] <- list(res[[i]])
                        names(res[[i]]) <- unlist(strsplit(input$inp[i, "name"],
                                                           split=".",
                                                           fixed=TRUE))[1]
                    }
                    else if (length(res[[i]]) == 6) {
                        if (checkmate::test_subset(names(res[[i]]), c("logFC",
                                                           "AveExpr",
                                                           "t",
                                                           "P.Value",
                                                           "adj.P.Val",
                                                           "B")) == FALSE) {
                            return(NULL)
                        }
                        res[[i]] <- deedee_prepare(res[[i]], "limma")
                        res[[i]] <- list(res[[i]])
                        names(res[[i]]) <- unlist(strsplit(input$inp[i, "name"],
                                                       split=".",
                                                       fixed=TRUE))[1]
                    }
                    else {
                        return(NULL)
                    }
                }
            # .xlsx input
            } else if (ext[[i]] == "xlsx") {
                sheets <- readxl::excel_sheets(input$inp[[i, "datapath"]])
                if (length(sheets) > 1) {
                    res[[i]] <- lapply(sheets,
                                       readxl::read_excel,
                                       path=input$inp[[i, "datapath"]])
                    names(res[[i]]) <- sheets
                    for (j in 1:length(sheets)) {
                        res[[i]][[sheets[j]]] <- as.data.frame(res[[i]][[sheets[j]]])
                        res[[i]][[sheets[j]]] <- tibble::column_to_rownames(
                            res[[i]][[sheets[j]]], "rowname")
                        if (checkmate::test_subset(names(res[[i]][[j]]),
                                                   c("logFC", "pval"))) {
                            return(NULL)
                        }
                    }
                }
                else {
                    res[[i]] <- readxl::read_excel(path=input$inp[[i, "datapath"]])
                    res[[i]] <- tibble::column_to_rownames(res[[i]], "rowname")
                    res[[i]] <- list(res[[i]])
                    names(res[[i]]) <- unlist(strsplit(input$inp[i, "name"],
                                                       split=".",
                                                       fixed=TRUE))[1]
                }
            # .txt input
            } else if (ext[[i]] == "txt") {
                temp <- read.table(input$inp[[i, "datapath"]])
                res[[i]] <- list(temp)
                names(res[[i]]) <- unlist(strsplit(input$inp[i, "name"],
                                                   split=".",
                                                   fixed=TRUE))[1]
            }
            else {
                return(NULL)
            }
        }

        # merging input data structures
        dat <- list()
        for(i in 1:length(res)) {
            for(j in 1:length(res[[i]])) {
                checkmate::assert_subset(names(res[[i]][[j]]), c("logFC", "pval"))
                dat[[names(res[[i]])[[j]]]] <- res[[i]][[j]]
            }
        }

        return(dat)
    })

    output$datasets <- renderUI({
        req(input$inp)
        checkboxGroupInput("select_datasets",
                           "Select datasets to be used",
                           choices = names(mydata()),
                           selected = names(mydata()))
    })

    mydata_use <- reactive({
        req(input$inp)
        use <- input$select_datasets
        dat2 <- list()
        for (i in use) {
            dat2[i] <- mydata()[i]
        }
        return(dat2)
    })

    output$inp_download <- downloadHandler(
        filename = "DeeDee_object.RDS",
        content = function(file) {
            req(input$inp)
            saveRDS(mydata_use(), file)
        })

    output$inp_infobox <- renderTable({
        validate(need(!is.null(mydata()),
                      'Faulty input data provided.'))

        req(input$inp)

        ext <- c()
        filename <- c()
        res <- list()
        type <- c()
        contrast <- c()
        genes <- c()
        count <- 0

        for(i in 1:length(input$inp[,1])) {
            ext[i] <- tools::file_ext(input$inp[i, "datapath"])

            if (ext[[i]] == "rds"|| ext[[i]] == "RDS") {
                res[[i]] <- readRDS(input$inp[[i, "datapath"]])

            } else if (ext[[i]] == "xlsx") {
                sheets <- readxl::excel_sheets(input$inp[[i, "datapath"]])
                res[[i]] <- lapply(sheets,
                                   readxl::read_excel,
                                   path=input$inp[[i, "datapath"]])
                names(res[[i]]) <- sheets
                for (j in 1:length(sheets)) {
                    res[[i]][[sheets[j]]] <- as.data.frame(res[[i]][[sheets[j]]])
                    res[[i]][[sheets[j]]] <- tibble::column_to_rownames(
                        res[[i]][[sheets[j]]], "rowname")
                }

            } else if (ext[[i]] == "txt") {
                temp <- read.table(input$inp[[i, "datapath"]])
                temp2 <- list()
                nm <- c()
                for (j in 0:((length(temp)/2)-1)) {
                    a <- 2*j+1
                    b <- 2*j+2
                    temp2[[j+1]] <- c(temp[a], temp[b])
                    temp2[[j+1]] <- as.data.frame(temp2[[j+1]])
                    nm[j+1] <- unlist(strsplit(names(temp)[2*j+1],
                                               split=".",
                                               fixed=TRUE))[1]
                }
                res[[i]] <- c(temp2[1:length(temp2)])
                names(res[[i]]) <- nm
                for (j in 1:((length(temp)/2))) {
                }
                for (j in 1:length(res[[i]])) {
                    names(res[[i]][[j]]) <- c("logFC", "pval")
                    row.names(res[[i]][[j]]) <- row.names(temp)
                }
            }
            if (class(res[[i]]) == "DESeqResults" ||
                class(res[[i]]) == "DGEExact" ||
                length(names(res[[i]])) == 6) {

                count <- count + 1
                type[count] <- class(res[[i]])
                filename[count] <- input$inp[i,"name"]
                contrast[count] <- unlist(strsplit(filename[count],
                                                   split=".",
                                                   fixed=TRUE))[1]
                genes[count] <- length(mydata()
                                       [[contrast[count]]][["logFC"]])

            } else {
                if (class(res[[i]]) == "data.frame") {
                    res[[i]] <- list(res[[i]])
                    names(res[[i]]) <- unlist(strsplit(input$inp[i, "name"],
                                                       split=".",
                                                       fixed=TRUE))[1]
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

        df <- data.frame(filename, type, contrast, genes)
        return(df)
    })


    # ------------------------------- scatter ----------------------------------
    # --- selectors ---
    output$scatter_choices1 <- renderUI({
        req(input$inp)
        selectInput("scatter_select1",
                    "1st dataset",
                    choices = names(mydata_use()))
    })

    output$scatter_choices2 <- renderUI({
        req(input$inp)
        selectInput("scatter_select2",
                    "2nd dataset",
                    selected = names(mydata_use())[2],
                    choices = names(mydata_use()))
    })

    # --- plot output ---
    ranges <- reactiveValues(x = NULL, y = NULL)

    output$scatter <- renderPlot({
        validate(
            need(input$inp,
                 'Please input at least two contrasts.')
        )
        validate(
            need(length(mydata()) >= 2,
                 'Please upload at least two contrasts.')
        )
        validate(
            need(length(mydata_use()) >= 2,
                 'Please select at least two contrasts.')
        )

        sel1 <- match(input$scatter_select1, names(mydata_use()))
        sel2 <- match(input$scatter_select2, names(mydata_use()))
        req(sel1)
        req(sel2)
        res <- deedee_scatter(mydata_use(),
                       select1 = sel1,
                       select2 = sel2,
                       color_by = input$scatter_color_by,
                       pthresh = input$scatter_pthresh)
        validate(
            need(!is.null(res), "No common genes in input datasets.")
        )
         res +
             ggplot2::coord_cartesian(xlim = ranges$x,
                                      ylim = ranges$y,
                                      expand = FALSE)
    })

    # --- brushing ---
    scatter_brushed <- reactive ({
        req(input$scatter_brush)
        df <- data.frame(x = mydata_use()[[input$scatter_select1]],
                         y = mydata_use()[[input$scatter_select2]])

        df <- subset(df, x.pval < 0.05 & y.pval < 0.05)

        names(df) <- c(paste(input$scatter_select1, ".logFC", sep = ""),
                       paste(input$scatter_select1, ".pval", sep = ""),
                       paste(input$scatter_select2, ".logFC", sep = ""),
                       paste(input$scatter_select2, ".pval", sep = ""))
        brushedPoints(df,
                      input$scatter_brush,
                      xvar = paste(input$scatter_select1, ".logFC", sep = ""),
                      yvar = paste(input$scatter_select2, ".logFC", sep = ""))
    })

    output$scatter_brush_info <- renderTable({
        req(scatter_brushed())
        scatter_brushed()}, rownames = TRUE)

    output$scatter_brush_download <- downloadHandler(
        filename = "scatter_brushed_genes.txt",
        content = function(file) {
            req(input$scatter_brush)
            write.table(scatter_brushed(), file)
    })

    observeEvent(input$scatter_dblclick, {
        brush <- input$scatter_brush
        if (!is.null(brush)) {
            ranges$x <- c(brush$xmin, brush$xmax)
            ranges$y <- c(brush$ymin, brush$ymax)

        } else {
            ranges$x <- NULL
            ranges$y <- NULL
        }
    })

    # --- enrich ---
    enrich <- reactive({
        req(input$inp)
        validate(need(scatter_brushed(),
                      "No brushed genes."))
        sel1 <- match(input$scatter_select1, names(mydata_use()))
        sel2 <- match(input$scatter_select2, names(mydata_use()))
        req(sel1)
        req(sel2)
        data <- list(mydata_use()[[sel1]], mydata_use()[[sel2]])
        res <- gsea(geneList = scatter_brushed(),
                    universe = data,
                    orgDB = input$organism,
                    select = 1)
        validate(need(class(res) == "enrichResult",
                      "Not working."))

        temp <- as.data.frame(res)
        if (dim(temp)[1] == 0) {
            return(NULL)
        }
        return(res)
    })

    output$scatter_ora <- renderPlot({
        req(!is.null(enrich()))
        en <- enrich()
        validate(need(class(en) == "enrichResult",
                      "No enriched terms found."))
        enrichplot::emapplot(enrichplot::pairwise_termsim(en))
    })

    output$ora_download <- downloadHandler(
        filename = "enrichment_results.RDS",
        content = function(file) {
            req(!is.null(enrich()))
            readr::writeRDS(enrich(), file)
    })

    # ------------------------------- heatmap ----------------------------------
    output$heatmap_errors <- renderText({
        validate(
            need(input$inp,
                 'Please input at least two contrasts.')
        )
        validate(
            need(length(mydata()) >= 2,
                 'Please upload at least two contrasts.')
        )
        validate(
            need(length(mydata_use()) >= 2,
                 'Please select at least two contrasts.')
        )
        validate(
            need(!is.null(heatmap_output()),
                 "No common genes in input datasets.")
        )
        ''
    })

    heatmap_output <- reactive({
        req(input$inp)
        req(input$heatmap_show_first)
        req(mydata_use())
        res <- deedee_heatmap(mydata_use(),
                              show_first = input$heatmap_show_first,
                              show_gene_names = input$heatmap_show_gene_names,
                              dist = input$heatmap_dist,
                              clust = input$heatmap_clust,
                              pthresh = input$heatmap_pthresh,
                              show_na = input$heatmap_showNA)
        validate(
            need(!is.null(res), "No common genes in input datasets.")
        )

        res <- ComplexHeatmap::draw(res)

        return(res)
    })

    listen <- reactive({
        list(input$heatmap_show_first,
          input$heatmap_show_gene_names,
          input$heatmap_dist,
          input$heatmap_clust,
          input$heatmap_pthresh,
          input$heatmap_showNA,
          mydata_use())
    })

    observeEvent(listen(), {
        req(length(mydata_use()) >= 2)
        InteractiveComplexHeatmap::makeInteractiveComplexHeatmap(
             input,
             output,
             session,
             heatmap_output())})


    # -------------------------------- venn ------------------------------------
    output$venn <- renderPlot({
        validate(
            need(input$inp,
                 'Please input at least two contrasts.')
        )
        validate(
            need(length(mydata()) >= 2,
                 'Please upload at least two contrasts.')
        )
        validate(
            need(length(mydata_use()) >= 2,
                 'Please select at least two contrasts.')
        )

        req(input$inp)
        res <- deedee_venn(mydata_use(),
                    mode = input$venn_mode,
                    pthresh = input$venn_pthresh)

        validate(
            need(!is.null(res),
                 "No genes in your datasets. Maybe your specified p-value threshold is too low?")
        )

        res
    })


    # -------------------------------- upSet -----------------------------------
    output$upSet <- renderPlot({
        validate(
            need(input$inp,
                 'Please input at least two contrasts.')
        )
        validate(
            need(length(mydata()) >= 2,
                 'Please upload at least two contrasts.')
        )
        validate(
            need(length(mydata_use()) >= 2,
                 'Please select at least two contrasts.')
        )

        req(input$inp)
        if (input$upSet_mode == "both" && input$upSet_colored) {
             mode = "both_colored"
        } else {
            mode = input$upSet_mode
        }
        res <- deedee_upSet(mydata_use(),
                     mode = mode,
                     pthresh = input$upSet_pthresh,
                     min_setsize = input$upSet_minset)

        validate(
            need(!is.null(res),
                 "No genes in your datasets. Maybe your specified p-value threshold is too low?")
        )

        res
    })


    # --------------------------------- qq -------------------------------------
    output$qq_choices1 <- renderUI ({
        req(input$inp)
        selectInput("qq_select1",
                    "1st dataset",
                    choices = names(mydata_use()))
    })

    output$qq_choices2 <- renderUI ({
        req(input$inp)
        selectInput("qq_select2",
                    "2nd dataset",
                    selected = names(mydata_use())[2],
                    choices = names(mydata_use()))
    })

    output$qq <- renderPlot({
        validate(
            need(input$inp,
                 'Please input at least two contrasts.')
        )
        validate(
            need(length(mydata()) >= 2,
                 'Please upload at least two contrasts.')
        )
        validate(
            need(length(mydata_use()) >= 2,
                 'Please select at least two contrasts.')
        )

        req(input$inp)
        sel1 <- match(input$qq_select1, names(mydata_use()))
        sel2 <- match(input$qq_select2, names(mydata_use()))
        req(sel1)
        req(sel2)
        res <- deedee_qq(mydata_use(),
                  select1 = sel1,
                  select2 = sel2,
                  color_by = input$qq_color_by,
                  pthresh = input$qq_pthresh)

        validate(
            need(!is.null(res),
                 "No genes in your datasets. Maybe your specified p-value threshold is too low?")
        )

        res
    })

    qq_brushed <- reactive ({
        req(input$qq_brush)
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
            sx <- approx(1L:lenx, sx, n = leny)$y
            pval1 <- approx(1L:lenx, pval1, n = leny)$y
        }
        if (leny > lenx) {
            sy <- approx(1L:leny, sy, n = lenx)$y
            pval2 <- approx(1L:leny, pval2, n = lenx)$y
        }

        sx <- tibble::rownames_to_column(as.data.frame(sx))
        sx[3] <- pval1
        sy <- tibble::rownames_to_column(as.data.frame(sy))
        sy[3] <- pval2

        qq <- data.frame(
            x = sx,
            y = sy)

        names(qq) <- c(
            paste(input$qq_select1, "gene", sep = "."),
            paste(input$qq_select1, "logFC", sep = "."),
            paste(input$qq_select1, "pval", sep = "."),
            paste(input$qq_select2, "gene", sep = "."),
            paste(input$qq_select2, "logFC", sep = "."),
            paste(input$qq_select2, "pval", sep = "."))

        temp1 <- paste(input$qq_select1, "logFC", sep = ".")
        temp2 <- paste(input$qq_select2, "logFC", sep = ".")

        brushedPoints(qq,
                      input$qq_brush,
                      xvar = temp1,
                      yvar = temp2)
    })

    output$qq_brush_info <- renderTable({
        req(qq_brushed())
        qq_brushed()}, rownames = FALSE)

    output$qq_brush_download <- downloadHandler(
        filename = "scatter_brushed_genes.txt",
        content = function(file) {
            req(input$qq_brush)
            write.table(qq_brushed(), file)
        })


    # --------------------------------- cat ------------------------------------
    output$cat_choice <- renderUI({
        req(input$inp)
        selectInput("cat_ref",
                    "Reference contrast",
                    selected = names(mydata_use())[1],
                    choices = names(mydata_use()))
    })

    output$cat <- renderPlot({
        validate(
            need(input$inp,
                 'Please input at least two contrasts.')
        )
        validate(
            need(length(mydata()) >= 2,
                 'Please upload at least two contrasts.')
        )
        validate(
            need(length(mydata_use()) >= 2,
                 'Please select at least two contrasts.')
        )

        req(input$inp)
        req(input$cat_maxrank)
        req(input$cat_ref)

        ref <- match(input$cat_ref, names(mydata_use()))
        res <- deedee_cat(mydata_use(),
                   ref = ref,
                   maxrank = input$cat_maxrank,
                   pthresh = input$cat_pthresh)

        validate(
            need(!is.null(res),
                 "No genes in your datasets. Maybe your specified p-value threshold is too low?")
        )

        res
    })
}


# ------------------------------ run application -------------------------------
shinyApp(ui = ui, server = server)
