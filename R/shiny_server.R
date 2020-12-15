library(shiny)

# Define UI for SangeR app
ui <- fluidPage(

  # App title
  titlePanel("SangeR"),

  # Sidebar layout with input and output definitions
  sidebarLayout(



    # Sidebar panel for inputs
    sidebarPanel(

      # Input: Select a file
      fileInput("file1", "Choose ab1 File",
                multiple = FALSE),

      # Input: Select a POI
      fileInput("file2", "Choose POI File",
                multiple = FALSE),

      # Horizontal line
      tags$hr(),
      tags$label("Read in parameter"),

      # Input for offset
      sliderInput(inputId = "offset",
                  label = "offset:",
                  min = 1,
                  max = 50,
                  value = 33),

      # Input for minimal sequence length
      sliderInput(inputId = "min_seq_len",
                  label = "minimum sequence length:",
                  min = 1,
                  max = 50,
                  value = 20),

      # Input for cutoff
      sliderInput(inputId = "cutoff",
                  label = "cutoff:",
                  min = 0.001,
                  max = 0.3,
                  value = 0.05),

      # Input for upstram region
      sliderInput(inputId = "upstream",
                  label = "upstream:",
                  min = 0,
                  max = 2000,
                  value = 500),

      # Horizontal line ----
      tags$hr(),
      tags$label("Host parameter"),

      #Input for host
      selectInput(
        inputId = "host",
        label = "host",
        choices = c("grch37.ensembl.org","www.ensembl.org"),
      ),

      #donwload button
      downloadButton('downloadPlot', 'Download Plot')

    ),

    # Main panel for displaying outputs
    mainPanel(

      # Output: Histogram ----
      plotOutput(outputId = "distPlot")

    )
  )
)

# Define server logic required to draw a histogram
server <- function(input, output) {

  # Histogram

  plotInput <- reactive({

    if(!is.null(input$file1)){
      histograms <- sangeR::plot_hist(sangeR::allign(sangeR::get_ref(sangeR::read.ab1(filename = input$file1$datapath, cutoff = input$cutoff, min_seq_len = input$min_seq_len, offset = input$offset),upstream = input$upstream, host = input$host)),POI = input$file2$datapath)

      do.call("grid.arrange", c(histograms, ncol=1))
    }

  })

  output$distPlot <- renderPlot({

    if(!is.null(input$file1)){
      print(plotInput())
    }
  })

  #download of png
  output$downloadPlot <- downloadHandler(
    filename = "Shinyplot.png",
    content = function(file) {
      png(file)
      print(plotInput())
      dev.off()
    })

}

# Create Shiny app ----
shinyApp(ui = ui, server = server)
