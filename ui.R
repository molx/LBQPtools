library(shiny)
#library(shinyjs)

subtitlesStyle <- "font-weight: bold; color: #000000;"

# Define UI for application that draws a histogram
shinyUI(
  tagList(shinyjs::useShinyjs(),
          tags$head(
            tags$style((HTML("
                              #shiny-notification-panel {
                                  position: fixed;
                                  bottom: 50%;
                                  right: 50%;
                              }
                              #progress-text {
                                  position: fixed;
                                  top: 50%;
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
                                            accept = c('FASTA', 
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
                                                       accept = ('.txt')),
                                             actionButton("bt_resetHeaderList", label = "Reset"),
                                             helpText("Sequence Length (set 0 to disable)"),
                                             div(style = "display:inline-block;",
                                                 numericInput("minSeqFilter", "Min", value = 0,
                                                              min = 0, max = 1e5, step = 1,
                                                              width = 100)),
                                             div(style = "display:inline-block;",
                                                 numericInput("maxSeqFilter", "Max", value = 0,
                                                              min = 0, max = 1e5, step = 1,
                                                              width = 100)),
                                             tags$br(),tags$br(),
                                             downloadButton("bt_doFilter", label = "Filter and save")
                                    ),
                                    tabPanel("Merge",
                                             tags$br(),
                                             fileInput('file_Fasta2', 'Choose FASTA file #2',
                                                       accept = c('FASTA', 
                                                                '.fasta')),
                                             checkboxInput("cb_rmFastaDups", label = "Remove Redundancies",
                                                           value = TRUE),
                                             checkboxInput("cb_rmFastaRevs", label = "Remove Reverses",
                                                           value = TRUE),
                                             checkboxInput("cb_rmFastaConts", label = "Remove Contaminants",
                                                           value = TRUE),
                                             div(style = "display:inline-block;width:150px;",
                                                 checkboxInput("cb_rmFastaTag", label = "Remove with tag:",
                                                               value = FALSE)),
                                             div(style = "display:inline-block;",
                                                 tags$input(id = "tx_tagRemove", class = "input-small")),
                                             tags$br(), tags$br(),
                                             downloadButton("bt_doMerge", label = "Merge and save")
                                    )
                                  )
                                )
                              )
                     ),
                     tabPanel("GO Tools",
                              titlePanel("GO Tools"),
                              sidebarLayout(
                                sidebarPanel(
                                  selectInput("cb_GOfileType", "File type", 
                                              choices = c(`Select one` = "wait",
                                                          `Comma separated values (US .csv)` = "csv1",
                                                          `Semicolon separated values (EU/BR .csv)` = "csv2",
                                                          `Tab separated values` = "tsv",
                                                          `Text file` = "txt"),
                                              selected = "wait"),
                                  fileInput("file_b2gTable", "Choose GO IDs file",
                                            accept = c(".csv", ".tsv", ".txt")),
                                  fluidRow(
                                    column(6,
                                           selectInput("cb_GOcol", label = "GO IDs column name", 
                                                       choices = c(`Load file first` = "wait"), 
                                                       selected = "wait")
                                    ),
                                    column(6#,
                                           #textInput("tx_GOsep", label = "GO IDs separator", value = "")
                                    )
                                  ),
                                  uiOutput("hp_NumGOLines"),
                                  helpText("Instructions:", style = subtitlesStyle),
                                  htmlOutput("vb_GOhelp")
                                ),
                                mainPanel(
                                  tabsetPanel(
                                    tabPanel("Descriptive Stastitics",
                                             tags$br(),
                                             actionButton("bt_GOdesc", "Analyze file"),
                                             #helpText("This function is currently only useful for GO IDs which have the category indicator at the beginning."),
                                             fluidRow(
                                               column(7,
                                                       tableOutput("df_GOcounts")
                                               ),
                                               column(5,
                                                     plotOutput("plot_GObar")
                                               )
                                             ),
                                             DT::dataTableOutput("dt_GOall")
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
                                             fluidRow(
                                               column(4,
                                                      selectInput("cb_GOfileTypeRef", "File type", 
                                                                  choices = c(`Select one` = "wait",
                                                                              `Comma separated values (US .csv)` = "csv1",
                                                                              `Semicolon separated values (EU/BR .csv)` = "csv2",
                                                                              `Tab separated values` = "tsv"),
                                                                  selected = "wait")
                                               ),
                                               column(4,
                                                      fileInput("file_b2gTableRef", "Choose reference GO IDs file",
                                                                accept = c(".csv", ".tsv", ".txt"))
                                               ), 
                                               column(4)),
                                             fluidRow(
                                               column(4,
                                                      selectInput("cb_GOcolRef", label = "GO IDs column name", 
                                                                  choices = c(`Load file first` = "wait"), 
                                                                  selected = "wait")
                                               ),
                                               column(4,
                                                      textInput("tx_GOsepRef", label = "GO IDs separator", value = "; ")
                                               ),
                                               column(4,
                                                      uiOutput("hp_NumGOLinesRef"))
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
                                                                  choices = c(`csv (US)` = "csv1",
                                                                              `csv (EU/BR)` = "csv2"))),
                                               column(5,
                                                      DT::dataTableOutput("tb_GOstat"))
                                             ),
                                             column(5)
                                    )
                                  )
                                )
                              )),
                     tabPanel("Enrichment Tester",
                              titlePanel("Enrichment Tester"),
                              sidebarLayout(
                                sidebarPanel(
                                  selectInput("se_ETfileType", "File type", 
                                              choices = c(`Select one` = "wait",
                                                          `Comma separated values (US .csv)` = "csv1",
                                                          `Semicolon separated values (EU/BR .csv)` = "csv2",
                                                          `Tab separated values` = "tsv"),
                                              selected = "wait"),
                                  fileInput("file_ET", "Choose file",
                                            accept = c(".csv", ".tsv", ".txt")),
                                  uiOutput("hp_NumETLines"),
                                  fluidRow(
                                    column(6,
                                           selectInput("se_ETcol", label = "IDs column name", 
                                                       choices = c(`Load file first` = "wait"), 
                                                       selected = "wait")
                                    ),
                                    column(6
                                    )
                                  ),
                                  fluidRow(
                                    column(4,
                                           radioButtons("rb_ETsep", label = NULL,
                                                        choices = c(`Multiple IDs per line` = "mult",
                                                                    `One ID per line` = "single"))
                                           ),
                                    column(6,
                                           textInput("tx_ETsep", label = "IDs separator", value = "")
                                           )
                                  ),
                                  fluidRow(
                                    column(4,
                                           radioButtons("rb_ETcountedData", label = NULL,
                                                        choices = c(`Count IDs` = "count",
                                                                    `IDs already counted` = "nocount"))
                                           ),
                                    column(6,
                                           selectInput("se_ETcountCol", label = "Select the count column", 
                                                       choices = c(`Load file first` = "wait"), 
                                                       selected = "wait")
                                           )
                                  ),
                                  radioButtons("rb_ETrefStyle", label = NULL,
                                               choices = c(`Reference data on the other file` = "other",
                                                           `Reference data on the same file` = "same"),
                                               inline = TRUE),
                                  conditionalPanel("input.rb_ETrefStyle == 'other'",
                                                   fluidRow(
                                                     column(6,
                                                            selectInput("se_ETfileTypeRef", "Reference File type", 
                                                               choices = c(`Select one` = "wait",
                                                                           `Comma separated values (US .csv)` = "csv1",
                                                                           `Semicolon separated values (EU/BR .csv)` = "csv2",
                                                                           `Tab separated values` = "tsv"),
                                                               selected = "wait")
                                                            ),
                                                     column(6,
                                                            fileInput("file_ETRef", "Choose reference file",
                                                                      accept = c(".csv", ".tsv", ".txt"))
                                                            )
                                                   ),
                                                   uiOutput("hp_NumETLinesRef")
                                  ),
                                  fluidRow(
                                    column(6,
                                           selectInput("se_ETcolRef", label = "Reference IDs column name", 
                                                       choices = c(`Load file first` = "wait"), 
                                                       selected = "wait")
                                    ),
                                    column(6
                                    )
                                  ),
                                  fluidRow(
                                    column(4,
                                           radioButtons("rb_ETsepRef", label = NULL,
                                                        choices = c(`Multiple IDs per line` = "mult",
                                                                    `One ID per line` = "single"))
                                    ),
                                    column(6,
                                           textInput("tx_ETsepRef", label = "Reference IDs separator", value = "")
                                    )
                                  ),
                                  fluidRow(
                                    column(4,
                                           radioButtons("rb_ETcountedDataRef", label = NULL,
                                                        choices = c(`Count IDs` = "count",
                                                                    `IDs already counted` = "nocount"))
                                    ),
                                    column(6,
                                           selectInput("se_ETcountColRef", label = "Select the reference count column", 
                                                       choices = c(`Load file first` = "wait"), 
                                                       selected = "wait")
                                    )
                                  )
                                ),
                                mainPanel(
                                  tabsetPanel(
                                    tabPanel("Results",
                                             tags$br(),
                                             fluidRow(
                                               column(2,
                                                      actionButton("bt_ETFisher", "Run tests!")
                                                      ),
                                               column(3,
                                                      downloadButton("bt_writeETstat", "Download results")
                                                      ),
                                               column(3,
                                                      selectInput("se_ETGOstatFormat", NA,
                                                                  choices = c(`csv (US)` = "csv1",
                                                                              `csv (EU/BR)` = "csv2",
                                                                              `Tab separated` = "tsv")))
                                             ),
                                             uiOutput("ET_colsCB", inline = TRUE),
                                             tags$br(),
                                             DT::dataTableOutput("tb_ETstat")),
                                    tabPanel("Instructions",
                                             htmlOutput("vb_EThelp"))
                                  )
                                )
                              )),
                     tabPanel(HTML("iTraq Ratio Check</a></li><li><a href='http://lbqp.unb.br/NetWheels/' target = '_blank'>NetWheels")),
                     tabPanel("About")
          )
  )
)
