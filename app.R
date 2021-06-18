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

    tabPanel("Input",
             fileInput("inp", "Upload DeeDee Input object (.RDS)",
                       multiple = TRUE,
                       accept = c(".rds", ".txt", ".xlsx")),
             tableOutput("files")),

    tabPanel("Scatterplot",
             sidebarPanel(
                 uiOutput("scatter_slider1"),

                 uiOutput("scatter_slider2"),

                 selectInput("scatter_color_by", h3("Color by"),
                             choices = list("1st p-value" = "pval1",
                                            "2nd p-value" = "pval2",
                                            "mean of p-values" = "pval_mean"),
                                            selected = "pval1")),

             mainPanel(plotOutput("scatter"))),

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

    tabPanel("Venn Diagram",
             sidebarPanel(
                 selectInput("venn_mode", h3("Mode"),
                             choices = list("Up" = "up",
                                            "Down" = "down",
                                            "Both" = "both"),
                             selected = "both")),

             mainPanel(plotOutput("venn"))),

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

    tabPanel("Quantile-Quantile Plot",
             sidebarPanel(
                 uiOutput("qq_slider1"),

                 uiOutput("qq_slider2"),

                 selectInput("qq_color_by", h3("Color by"),
                             choices = list("1st p-value" = "pval1",
                                            "2nd p-value" = "pval2"),
                             selected = "pval1")),

             mainPanel(plotOutput("qq"))),

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

    output$files <- renderTable(input$inp)

    mydata <- reactive({
        req(input$inp)

        ext <- vector(mode = "numeric", length = length(input$inp[,1]))
        res <- list()

        for(i in 1:length(input$inp[,1])) {
            ext[i] <- tools::file_ext(input$inp[i, "datapath"])
            validate(need(ext[[i]] == "rds"||
                              ext[[i]] == "RDS" ||
                              ext[[i]] == "xlsx" ||
                              ext[[i]] == "txt",
                          "Please upload only .RDS, .xlsx or .txt files"))
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
        }

        dat <- list()
        for(i in 1:length(res)) {
            for(j in 1:length(res[[i]])) {
                assert_subset(names(res[[i]][[j]]), c("logFC", "pval"))
                dat[[names(res[[i]])[[j]]]] <- res[[i]][[j]]
            }
        }

        return(dat)})

    output$scatter_slider1 <- renderUI({
        numericInput("scatter_select1",
                     h3("1st dataset"),
                     value = 1,
                     min = 1,
                     max = length(mydata()))
    })

    output$scatter_slider2 <- renderUI({
        numericInput("scatter_select2",
                     h3("2nd dataset"),
                     value = 2,
                     min = 1,
                     max = length(mydata()))
    })


    output$scatter <- renderPlot({
                deedee_scatter(mydata(),
                               select1 = input$scatter_select1,
                               select2 = input$scatter_select2,
                               color_by = input$scatter_color_by)
    })

    output$heatmap <- renderPlot({
        deedee_heatmap(mydata(),
                       show_first = input$heatmap_show_first,
                       show_gene_names = input$heatmap_show_gene_names)
    })

    output$venn <- renderPlot({
        deedee_venn(mydata(),
                    mode = input$venn_mode)
    })

    output$upSet <- renderPlot({
        if (input$upSet_mode == "both" && input$upSet_colored) {
             mode = "both_colored"
        } else {
            mode = input$upSet_mode
        }
        deedee_upSet(mydata(),
                     mode = mode)
    })

    output$qq_slider1 <- renderUI ({
        numericInput("qq_select1",
                     h3("1st dataset"),
                     value = 1,
                     min = 1,
                     max = length(mydata()))
    })

    output$qq_slider2 <- renderUI ({
        numericInput("qq_select2",
                     h3("2nd dataset"),
                     value = 2,
                     min = 1,
                     max = length(mydata()))
    })

    output$qq <- renderPlot({
        deedee_qq(mydata(),
                  select1 = input$qq_select1,
                  select2 = input$qq_select2,
                  color_by = input$qq_color_by)
    })

    output$cat <- renderPlot({
        deedee_cat2(mydata(),
                    maxrank = input$cat_maxrank)
    })
}


# Run the application
shinyApp(ui = ui, server = server)
