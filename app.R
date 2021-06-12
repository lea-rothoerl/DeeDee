library(shiny)
source("~/Development/DeeDee_wip/deedee_venn.R")
source("~/Development/DeeDee_wip/deedee_upSet.R")
source("~/Development/DeeDee_wip/deedee_scatter.R")
source("~/Development/DeeDee_wip/deedee_heatmap.R")
source("~/Development/DeeDee_wip/deedee_qq.R")
source("~/Development/DeeDee_wip/deedee_cat2.R")

# Define UI for application that draws a histogram
ui <- navbarPage("DeeDee",

    tabPanel("Input",
             fileInput("inp", "Upload DeeDee Input object (.RDS)",
                       multiple = TRUE,
                       accept = c(".rds"))),

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
                              value = 500,
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
                                            "2nd p-value" = "pval2",
                                            "1st logFC" = "logFC1",
                                            "2nd logFC" = "logFC2"),
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

# Define server logic required to draw a histogram
server <- function(input, output) {

    mydata <- reactive({
        file <- input$inp
        ext <- tools::file_ext(file$datapath)

        req(file)
        validate(need(ext == "rds"|| ext == "RDS", "Please upload a .RDS file"))
        mydata <- readRDS(file$datapath)

        return(mydata)})

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
