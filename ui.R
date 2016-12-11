library(shiny)
library(shinyjs)

subtitlesStyle <- "font-weight: bold; color: #000000;"

numericInput2 <- function(inputId, label, value, min = NA, max = NA, step = NA, 
                          width = NULL) {
  value <- restoreInput(id = inputId, default = value)
  inputTag <- tags$input(id = inputId, type = "number", class = "form-control", 
                         value = formatNoSci(value))
  if (!is.na(min)) 
    inputTag$attribs$min = min
  if (!is.na(max)) 
    inputTag$attribs$max = max
  if (!is.na(step)) 
    inputTag$attribs$step = step
  div(class = "form-group shiny-input-container", style = if (!is.null(width)) 
    paste0("width: ", validateCssUnit(width), ";"), label %AND% 
      tags$label(label, `for` = inputId), inputTag)
}

# Define UI for application that draws a histogram
shinyUI(
  tagList(useShinyjs(),
          navbarPage("LBQP Tools",
                     tabPanel("Fasta Tools",
                              tags$head(tags$style("#preview{font-family: monospace;}"
                              )
                              ),
                              titlePanel("Fasta Tools"),
                              sidebarLayout(
                                sidebarPanel(
                                  fileInput('file1', 'Choose FASTA file',
                                            accept=c('FASTA', 
                                                     '.fasta')),
                                  selectInput("seqtype", "Type of sequence", 
                                              choices = c(Aminoacids = "AA", Nucleotides = "DNA"), selected = "AA"),
                                  uiOutput("hp_NumOfSeq")
                                ),
                                mainPanel(
                                  tabsetPanel(
                                    tabPanel("Filters",
                                             tags$br(),
                                             helpText("Header List"),
                                             fileInput('fileHeaderFilter', 'Choose Header list file',
                                                       accept=c('.txt')),
                                             actionButton("bt_resetHeaderList", label = "Reset"),
                                             helpText("Sequence Length (set 0 to disable)"),
                                             div(style="display:inline-block;",
                                                 numericInput("minSeqFilter", "Min", value = 0,
                                                              min = 0, max = 1e5, step = 1,
                                                              width = 100)),
                                             div(style="display:inline-block;",
                                                 numericInput("maxSeqFilter", "Max", value = 0,
                                                              min = 0, max = 1e5, step = 1,
                                                              width = 100)),
                                             tags$br(),tags$br(),
                                             downloadButton("bt_doFilter", label = "Filter and save")
                                    ),
                                    tabPanel("Merge",
                                             tags$br(),
                                             fileInput('file2', 'Choose FASTA file #2',
                                                       accept=c('FASTA', 
                                                                '.fasta')),
                                             checkboxInput("cb_RemoveRedundancies", label = "Remove Redundancies",
                                                           value = TRUE),
                                             checkboxInput("cb_RemoveReverses", label = "Remove Reverses",
                                                           value = TRUE),
                                             checkboxInput("cb_RemoveContaminants", label = "Remove Contaminants",
                                                           value = TRUE),
                                             div(style="display:inline-block;width:150px;",
                                                 checkboxInput("cb_RemoveTag", label = "Remove with tag:",
                                                               value = FALSE)),
                                             div(style="display:inline-block;",
                                                 tags$input(id = "tx_tagRemove", class="input-small")),
                                             tags$br(),tags$br(),
                                             downloadButton("bt_doMerge", label = "Merge and save")
                                    )
                                  )
                                  # tags$h3("Output Preview"),
                                  # actionButton("bt_genPreview", label = "Generate Preview"),
                                  # fluidRow(
                                  #   column(6,
                                  #          sliderInput(inputId =  "nPreview", label = "Number of sequences in preview", 
                                  #                      min = 1, max = 20, value = 10)
                                  #   ),
                                  #   column(6,
                                  #          sliderInput(inputId =  "nbcharPreview", label = "Characters per line", 
                                  #                      min = 20, max = 100, value = 60, step = 5)
                                  #   )
                                  # ),
                                  # hr(),
                                  # htmlOutput("preview")
                                )
                              )
                     ),
                     tabPanel("GO Tools"),
                     tabPanel(HTML("iTraq Ratio Check</a></li><li><a href='http://lbqp.unb.br/NetWheels/' target = '_blank'>NetWheels")),
                     tabPanel("About")
          )
  )
)
