# ---------------------------- --- source packages ------------------------------
library(shiny)
source("~/Development/DeeDee_wip/deedee_venn.R")
source("~/Development/DeeDee_wip/deedee_upSet.R")
source("~/Development/DeeDee_wip/deedee_scatter.R")
source("~/Development/DeeDee_wip/deedee_heatmap.R")
source("~/Development/DeeDee_wip/deedee_qq.R")
source("~/Development/DeeDee_wip/deedee_cat2.R")
source("~/Development/DeeDee_wip/deedee_prepare.R")

# ------------------------------------------------------------------------------
# --------------------------------- U I ----------------------------------------
# ------------------------------------------------------------------------------

ui <- navbarPage("DeeDee",

    # ----------------------------- data input ---------------------------------
    tabPanel("Input",
             fileInput("inp", "Upload",
                       multiple = TRUE,
                       accept = c(".rds", ".txt", ".xlsx")),
             uiOutput("datasets"),
             downloadButton("inp_download", "Download DeeDee object (.RDS)")),


    # ------------------------------- scatter ----------------------------------
    tabPanel("Scatterplot",
             sidebarPanel(
                 uiOutput("scatter_choices1"),

                 uiOutput("scatter_choices2"),

                 selectInput("scatter_color_by", h3("Color by"),
                             choices = list("1st p-value" = "pval1",
                                            "2nd p-value" = "pval2",
                                            "mean of p-values" = "pval_mean"),
                                            selected = "pval1")),

             mainPanel(plotOutput("scatter",
                       click = "scatter_click")),
                 verbatimTextOutput("scatter_click_info")),


    # ------------------------------- heatmap ----------------------------------
    tabPanel("Heatmap",
             sidebarPanel(
                 numericInput("heatmap_show_first",
                              h3("Show first"),
                              value = 25,
                              min = 1),

                 checkboxInput("heatmap_show_gene_names",
                               "Show gene names",
                               value = FALSE)),

             mainPanel(plotOutput("heatmap"))),


    # -------------------------------- venn ------------------------------------
    tabPanel("Venn Diagram",
             sidebarPanel(
                 selectInput("venn_mode", h3("Mode"),
                             choices = list("Up" = "up",
                                            "Down" = "down",
                                            "Both" = "both"),
                             selected = "both")),

             mainPanel(plotOutput("venn"))),


    # -------------------------------- upSet -----------------------------------
    tabPanel("UpSet Plot",
             sidebarPanel(
                 selectInput("upSet_mode", h3("Mode"),
                             choices = list("Up" = "up",
                                            "Down" = "down",
                                            "Both" = "both"),
                             selected = "both"),
                 conditionalPanel(
                     condition = "input$upSet_mode == both",
                     checkboxInput("upSet_colored",
                                   "Coloring",
                                   TRUE))),

             mainPanel(plotOutput("upSet"))),


    # --------------------------------- qq -------------------------------------
    tabPanel("Quantile-Quantile Plot",
             sidebarPanel(
                 uiOutput("qq_choices1"),

                 uiOutput("qq_choices2"),

                 selectInput("qq_color_by", h3("Color by"),
                             choices = list("1st p-value" = "pval1",
                                            "2nd p-value" = "pval2"),
                             selected = "pval1")),

             mainPanel(plotOutput("qq"))),


    # --------------------------------- cat ------------------------------------
    tabPanel("Concordance At the Top Plot",
             sidebarPanel(
                numericInput("cat_maxrank",
                            h3("Max rank"),
                            value = 1000,
                            min = 1)),
             mainPanel(plotOutput("cat")))
)

# ------------------------------------------------------------------------------
# ----------------------------- S E R V E R ------------------------------------
# ------------------------------------------------------------------------------

server <- function(input, output) {

    # ----------------------------- data input ---------------------------------
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
                    res[[i]] <- list(as.data.frame(c(res[[i]])))
                    names(res[[i]]) <- unlist(strsplit(input$inp[i, "name"],
                                                       split=".",
                                                       fixed=TRUE))[1]
                } else if (class(res[[i]]) == "edgeR") {
                    res[[i]] <- deedee_prepare(res[[i]], "edgeR")
                    res[[i]] <- list(as.data.frame(c(res[[i]])))
                    names(res[[i]]) <- unlist(strsplit(input$inp[i, "name"],
                                                       split=".",
                                                       fixed=TRUE))[1]
                } else if (class(res[[i]]) == "limma") {
                    res[[i]] <- deedee_prepare(res[[i]], "limma")
                    res[[i]] <- list(as.data.frame(c(res[[i]])))
                    names(res[[i]]) <- unlist(strsplit(input$inp[i, "name"],
                                                       split=".",
                                                       fixed=TRUE))[1]
                }

            # .xlsx input
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

            # .txt input
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
        }

        # merging input data structures
        dat <- list()
        for(i in 1:length(res)) {
            for(j in 1:length(res[[i]])) {
                assert_subset(names(res[[i]][[j]]), c("logFC", "pval"))
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


    # ------------------------------- scatter ----------------------------------
    output$scatter_choices1 <- renderUI({
        req(input$inp)
        selectInput("scatter_select1",
                    h3("1st dataset"),
                    choices = names(mydata_use()))
    })

    output$scatter_choices2 <- renderUI({
        req(input$inp)
        selectInput("scatter_select2",
                    h3("2nd dataset"),
                    selected = names(mydata_use())[2],
                    choices = names(mydata_use()))
    })

    output$scatter <- renderPlot({
        req(input$inp)
        sel1 <- match(input$scatter_select1, names(mydata_use()))
        sel2 <- match(input$scatter_select2, names(mydata_use()))
        req(sel1)
        req(sel2)
        deedee_scatter(mydata_use(),
                       select1 = sel1,
                       select2 = sel2,
                       color_by = input$scatter_color_by)
    })

    output$scatter_click_info <- renderPrint({
        req(input$scatter_click)
        df <- data.frame(x = mydata_use()[[input$scatter_select1]],
                         y = mydata_use()[[input$scatter_select2]])
        names(df) <- c(paste(input$scatter_select1, ".logFC", sep = ""),
                       paste(input$scatter_select1, ".pval", sep = ""),
                       paste(input$scatter_select2, ".logFC", sep = ""),
                       paste(input$scatter_select2, ".pval", sep = ""))
        nearPoints(df,
                   input$scatter_click,
                   xvar = paste(input$scatter_select1, ".logFC", sep = ""),
                   yvar = paste(input$scatter_select2, ".logFC", sep = ""),
                   addDist = FALSE)
    })


    # ------------------------------- heatmap ----------------------------------
    output$heatmap <- renderPlot({
        req(input$inp)
        deedee_heatmap(mydata_use(),
                       show_first = input$heatmap_show_first,
                       show_gene_names = input$heatmap_show_gene_names)
    })


    # -------------------------------- venn ------------------------------------
    output$venn <- renderPlot({
        req(input$inp)
        deedee_venn(mydata_use(),
                    mode = input$venn_mode)
    })


    # -------------------------------- upSet -----------------------------------
    output$upSet <- renderPlot({
        req(input$inp)
        if (input$upSet_mode == "both" && input$upSet_colored) {
             mode = "both_colored"
        } else {
            mode = input$upSet_mode
        }
        deedee_upSet(mydata_use(),
                     mode = mode)
    })


    # --------------------------------- qq -------------------------------------
    output$qq_choices1 <- renderUI ({
        req(input$inp)
        selectInput("qq_select1",
                    h3("1st dataset"),
                    choices = names(mydata_use()))
    })

    output$qq_choices2 <- renderUI ({
        req(input$inp)
        selectInput("qq_select2",
                    h3("2nd dataset"),
                    selected = names(mydata_use())[2],
                    choices = names(mydata_use()))
    })

    output$qq <- renderPlot({
        req(input$inp)
        sel1 <- match(input$qq_select1, names(mydata_use()))
        sel2 <- match(input$qq_select2, names(mydata_use()))
        req(sel1)
        req(sel2)
        deedee_qq(mydata_use(),
                  select1 = sel1,
                  select2 = sel2,
                  color_by = input$qq_color_by)
    })


    # --------------------------------- cat ------------------------------------
    output$cat <- renderPlot({
        req(input$inp)
        deedee_cat2(mydata_use(),
                    maxrank = input$cat_maxrank)
    })
}


# ------------------------------ run application -------------------------------
shinyApp(ui = ui, server = server)
