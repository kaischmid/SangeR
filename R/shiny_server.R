library(gridExtra)

#' @name ui
#' @title User Interface
#' @export
#'



# Define UI for SangeR app
ui <- shiny::fluidPage(

  # App title
  shiny::titlePanel("SangeR"),

  # Sidebar layout with input and output definitions
  shiny::sidebarLayout(



    # Sidebar panel for inputs
    shiny::sidebarPanel(

      # Input: Select a file
      shiny::fileInput("file1", "Choose ab1 File",
                multiple = FALSE),

      #Input: delimiter
      shiny::textInput(inputId = "delimiter",
                       label = "delimiter in ab1 file", value = "_", placeholder = "_"),

      # Input position of genename
      shiny::sliderInput(inputId = "genename_pos",
                  label = "genename Position:",
                  min = 1,
                  max = 5,
                  value = 1),

      # Input: Select a POI
      shiny::fileInput("file2", "Choose POI File",
                multiple = FALSE),

      # Horizontal line
      shiny::tags$hr(),
      shiny::tags$label("Read in parameter"),

      # Input for offset
      shiny::sliderInput(inputId = "offset",
                  label = "offset:",
                  min = 1,
                  max = 50,
                  value = 33),

      # Input for minimal sequence length
      shiny::sliderInput(inputId = "min_seq_len",
                  label = "minimum sequence length:",
                  min = 1,
                  max = 50,
                  value = 20),

      # Input for cutoff
      shiny::sliderInput(inputId = "cutoff",
                  label = "cutoff:",
                  min = 0.001,
                  max = 0.3,
                  value = 0.05),

      # Input for upstram region
      shiny::sliderInput(inputId = "upstream",
                  label = "upstream:",
                  min = 0,
                  max = 2000,
                  value = 500),

      # Horizontal line ----
      shiny::tags$hr(),
      shiny::tags$label("Host parameter"),

      #Input for host
      shiny::selectInput(
        inputId = "host",
        label = "host",
        choices = c("grch37.ensembl.org","www.ensembl.org"),
      ),

      #donwload button
      shiny::downloadButton('downloadPlot', 'Download Plot')

    ),

    # Main panel for displaying outputs
    shiny::mainPanel(

      # Output: Histogram ----
      shiny::plotOutput(outputId = "distPlot")

    )
  )
)

#' @name server
#' @title Server
#' @param input input for shiny ui
#' @param output output from shiny ui
#' @export
#'

# Define server logic required to draw a histogram
server <- function(input, output) {

  # Histogram

  plotInput <- shiny::reactive({

    if(!is.null(input$file1)){
      histograms <- sangeR::plot_hist(sangeR::allign(sangeR::get_ref(sangeR::read.ab1(filename = input$file1$datapath, delimiter = input$delimiter, genename_pos = input$genename_pos, cutoff = input$cutoff, min_seq_len = input$min_seq_len, offset = input$offset),upstream = input$upstream, host = input$host)))

      if(!is.null(histograms$PNG_list)){
        do.call("grid.arrange", c(histograms$PNG_list, ncol=1))
      }
    }
    shiny::validate(
      shiny::need(!is.null(histograms$PNG_list), "Could not find any mutations. Check in for a POI or less strict parameters.")
    )

  })



  output$distPlot <- shiny::renderPlot({

    if(!is.null(input$file1)){
      print(plotInput())
    }
  })

  #download of png
  output$downloadPlot <- shiny::downloadHandler(
    filename = "Shinyplot.png",
    content = function(file) {
      grDevices::png(file)
      print(plotInput())
      grDevices::dev.off()
    })

}

# Create Shiny app ----
shiny::shinyApp(ui = ui, server = server)
