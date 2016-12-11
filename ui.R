library(shiny)
library(shinyjs)

subtitlesStyle <- "font-weight: bold; color: #000000;"

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
                                  fileInput('file_Fasta1', 'Choose FASTA file',
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
                                             fileInput('file_Fasta2', 'Choose FASTA file #2',
                                                       accept=c('FASTA', 
                                                                '.fasta')),
                                             checkboxInput("cb_rmFastaDups", label = "Remove Redundancies",
                                                           value = TRUE),
                                             checkboxInput("cb_rmFastaRevs", label = "Remove Reverses",
                                                           value = TRUE),
                                             checkboxInput("cb_rmFastaConts", label = "Remove Contaminants",
                                                           value = TRUE),
                                             div(style="display:inline-block;width:150px;",
                                                 checkboxInput("cb_rmFastaTag", label = "Remove with tag:",
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
                     tabPanel("GO Tools",
                              titlePanel("GO Tools"),
                              sidebarLayout(
                                sidebarPanel(
                                  fileInput('file_b2gTable', 'Choose Blast2GO table file',
                                            accept=c('Text file', 
                                                     '.txt')),
                                  uiOutput("hp_NumGOLines"),
                                  helpText("File should be a .txt tab-delimited table, with GO IDs on a column named 'GO IDs list'")
                                ),
                                mainPanel(
                                  tabsetPanel(
                                    tabPanel("Extract IDs from table",
                                             checkboxInput("cb_rmGOcat", "Remove Category from GO IDs ('C:', 'F:' or 'P:')", 
                                                           value = TRUE, width = "500"),
                                             checkboxInput("cb_rmGOdups", "Remove Redundancies", 
                                                           value = FALSE, width = "500"),
                                             downloadButton("bt_extractGO_all", "Extract All GO IDs"), tags$br(), tags$br(),
                                             downloadButton("bt_extractGO_C", "Extract Cellular Component IDs"), tags$br(), tags$br(),
                                             downloadButton("bt_extractGO_F", "Extract Molecular Function IDs"), tags$br(), tags$br(),
                                             downloadButton("bt_extractGO_P", "Extract Biological Process IDs")
                                             ),
                                    tabPanel("GO Statistics",
                                             plotOutput("plot_GOpie"))
                                  )
                                )
                              )),
                     tabPanel(HTML("iTraq Ratio Check</a></li><li><a href='http://lbqp.unb.br/NetWheels/' target = '_blank'>NetWheels")),
                     tabPanel("About")
          )
  )
)
