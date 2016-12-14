library(shiny)
library(shinyjs)

subtitlesStyle <- "font-weight: bold; color: #000000;"

# Define UI for application that draws a histogram
shinyUI(
  tagList(useShinyjs(),
          tags$head(
            tags$style((HTML("
                              #shiny-notification-panel {
                                  position: fixed;
                                  bottom: 50%;
                                  right: 50%;
                              }"
            )))
          ),
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
                                  fileInput("file_b2gTable", "Choose GO IDs file",
                                            accept=c(".csv", ".tsv", ".txt")),
                                  selectInput("cb_GOfileType", "File type", 
                                              choices = c(`Comma separated values (US .csv)` = "csv1",
                                                          `Semicolon separated values (BR .csv)` = "csv2",
                                                          `Tab separated values` = "tsv"),
                                              selected = "tsv"),
                                  fluidRow(
                                    column(6,
                                           textInput("tx_GOcol", label = "GO IDs column name", value = "GO IDs list")
                                    ),
                                    column(6,
                                           textInput("tx_GOsep", label = "GO IDs separator", value = "; ")
                                    )
                                  ),
                                  uiOutput("hp_NumGOLines"),
                                  helpText(""),
                                  helpText("For example, one cell of the column could be 'P:GO:0051297; C:GO:0031252; C:GO:0036064; P:GO:0035020'")
                                ),
                                mainPanel(
                                  tabsetPanel(
                                    tabPanel("Descriptive Stastitics",
                                             tags$br(),
                                             actionButton("bt_GOdesc", "Load file"),
                                             fluidRow(
                                               column(4,
                                                      plotOutput("plot_GOpie")
                                               ),
                                               column(3,
                                                      tableOutput("df_GOcounts")
                                               ),
                                               column(5,
                                                      DT::dataTableOutput("dt_GOall")
                                               )
                                             )
                                    ),
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
                                    tabPanel("Enrichment test",
                                             fileInput("file_b2gTableRef", "Choose reference GO IDs file",
                                                       accept=c(".csv", ".tsv", ".txt")),
                                             fluidRow(
                                               column(4,
                                                      selectInput("cb_GOfileTypeRef", "File type", 
                                                                  choices = c(`Comma separated values (US .csv)` = "csv1",
                                                                              `Semicolon separated values (BR .csv)` = "csv2",
                                                                              `Tab separated values` = "tsv"),
                                                                  selected = "tsv")),
                                               column(4,
                                                      textInput("tx_GOcolRef", label = "GO IDs column name", value = "GO IDs list")
                                               ),
                                               column(4,
                                                      textInput("tx_GOsepRef", label = "GO IDs separator", value = "; ")
                                               )
                                             ),
                                             fluidRow(
                                               column(2,
                                                      actionButton("bt_GOFisher" ,"Run Fisher's Exact Test")
                                               ),
                                               column(10,
                                                      helpText("(Wait until the table appears below to continue)"))
                                             ),
                                             tags$br(),
                                             fluidRow(
                                               column(2,
                                                      downloadButton("bt_writeGOstat", "Download Results"),
                                                      tags$br(),
                                                      selectInput("cb_GOstatFormat", "Output format",
                                                                  choices = c(`csv (US)`= "csv1",
                                                                              `csv (BR)`= "csv2"))),
                                               column(5,
                                                      DT::dataTableOutput("tb_GOstat"))
                                             ),
                                             column(5)
                                    )
                                  )
                                )
                              )),
                     tabPanel(HTML("iTraq Ratio Check</a></li><li><a href='http://lbqp.unb.br/NetWheels/' target = '_blank'>NetWheels")),
                     tabPanel("About")
          )
  )
)
