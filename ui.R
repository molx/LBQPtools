library(shiny)

subtitlesStyle <- "font-weight: bold; color: #000000;"

# Define UI for application that draws a histogram
shinyUI(fluidPage(
  tags$head(tags$style("#preview{font-family: monospace;}"
                         )
  ),
  titlePanel("Fasta Tools"),
  
  
  
  # Sidebar with a slider input for number of bins 
  sidebarLayout(
    sidebarPanel(
      fileInput('file1', 'Choose FASTA File #1',
                accept=c('FASTA', 
                         '.fasta')),
      selectInput("seqtype", "Type of sequence", 
                   choices = c(Aminoacids = "AA", Nucleotides = "DNA"), selected = "AA"),
      helpText("Number of Sequences", style = subtitlesStyle),
      fluidRow(
        column(6,
               uiOutput("nBefore")
        ),
        column(6,
               uiOutput("nAfter")
        )
      ),
      tags$hr(),
      helpText("Filters", style = subtitlesStyle),
      helpText("Sequence Length (set 0 to disable)"),
      fluidRow(
        column(6,
               numericInput("minSeqFilter", "Min", value = 0,
                            min = 0, max = 1e5, step = 1)
               ),
        column(6,
               numericInput("maxSeqFilter", "Max", value = 0,
                            min = 0, max = 1e5, step = 1)
        )
      ),
      helpText("Merge", style = subtitlesStyle),
      fileInput('file2', 'Choose FASTA File #2',
                accept=c('FASTA', 
                         '.fasta')),
      checkboxInput("cb_RemoveRedundancies", label = "Remove Redundancies",
                    value = TRUE),
      checkboxInput("cb_RemoveReverses", label = "Remove Reverses",
                    value = TRUE),
      checkboxInput("cb_RemoveContaminants", label = "Remove Contaminants",
                    value = TRUE),
      fluidRow(
        column(6,
               checkboxInput("cb_RemoveTag", label = "Remove with tag:",
                             value = TRUE)
               ),
        column(6,
               textInput("tx_tagRemove", label = "Tag:")
               )
      ),
      downloadButton("bt_doMerge", label = "Merge Fasta files")
    ),
    # Show a plot of the generated distribution
    mainPanel(
      tags$h3("Output Preview"),
      fluidRow(
        column(6,
               sliderInput(inputId =  "nPreview", label = "Number of sequences in preview", 
                           min = 1, max = 20, value = 10)
               ),
        column(6,
               sliderInput(inputId =  "nbcharPreview", label = "Characters per line", 
                           min = 20, max = 100, value = 60, step = 5)
        )
      ),
      hr(),
      htmlOutput("preview")
    )
  )
))
