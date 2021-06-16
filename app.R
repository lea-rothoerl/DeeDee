library(shiny)
source("~/Development/DeeDee_wip/deedee_venn.R")
source("~/Development/DeeDee_wip/deedee_upSet.R")
source("~/Development/DeeDee_wip/deedee_scatter.R")
source("~/Development/DeeDee_wip/deedee_heatmap.R")
source("~/Development/DeeDee_wip/deedee_qq.R")
source("~/Development/DeeDee_wip/deedee_cat2.R")

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
                 numericInput("scatter_select1",
                              h3("1st dataset"),
                              value = 1,
                              min = 1),

                 numericInput("scatter_select2",
                              h3("2nd dataset"),
                              value = 2,
                              min = 1),

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
                                            "Both" = "both",
                                            "Both colored" = "both_colored"),
                             selected = "both")),

             mainPanel(plotOutput("upSet"))),

    tabPanel("Quantile-Quantile Plot",
             sidebarPanel(
                 numericInput("qq_select1",
                              h3("1st dataset"),
                              value = 1,
                              min = 1),

                 numericInput("qq_select2",
                              h3("2nd dataset"),
                              value = 2,
                              min = 1),

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

            } else if (ext[[i]] == "xlsx") {
                sheets <- readxl::excel_sheets(input$inp[[i, "datapath"]])
                res[[i]] <- lapply(sheets,
                                    readxl::read_excel,
                                    path=input$inp[[i, "datapath"]])
                names(res[[i]]) <- sheets
                for (j in 1:length(sheets)) {
                    res[[i]][[sheets[j]]] <- as.data.frame(res[[i]][[sheets[j]]])
                }

            } else if (ext[[i]] == "txt") {
                res[[i]] <- read.table(input$inp[[i, "datapath"]])
                View(res[[i]])
            }
        }

        dat <- list()
        for(i in 1:length(res)) {
            for(j in 1:length(res[[i]])) {
                dat[[names(res[[i]])[[j]]]] <- res[[i]][[j]]
            }
        }

        return(dat)})

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
        deedee_upSet(mydata(),
                     mode = input$upSet_mode)
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
