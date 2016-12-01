library(shiny)
#library(seqinr)

# Define server logic required to draw a histogram
shinyServer(function(input, output) {
   
  inFile1 <- reactive(input$file1)
  
  inFile2 <- reactive(input$file2)
  
# output$preview <- renderUI({
#     
#     if (is.null(inFile1()))
#       return(NULL)
#     
#     #fastaFile <- read.fasta("../proteins grandiflora.fasta")
#     fastaFile <- inFile1()
#     
#     output$nBefore <- renderUI(helpText(paste("Original:", length(fastaFile))))
#     
#     sequences <- seqinr::read.fasta(fastaFile$datapath, seqtype = input$seqtype,
#                             as.string = FALSE)
#     
#     sequences <- filter.fasta(sequences, 
#                               minRes = input$minSeqFilter, maxRes = input$maxSeqFilter)
#     
#     output$nAfter <- renderUI(helpText(paste("Filtered:", length(sequences))))
#     
#     if (length(sequences) != 0) {
#       sequences <- sequences[1:(min(length(sequences), input$nPreview))]
#       
#       header <- unlist(getAnnot(sequences))
#       
#       if (!is.list(sequences)) {
#         HTML(write.oneseq(pep = sequences, name = header, nbchar = input$nbcharPreview))
#       } else {
#         n.seq <- length(sequences)
#         HTML(paste0(sapply(seq_len(n.seq), function(x) {
#           write.oneseq(pep = as.character(sequences[[x]]), 
#                        name = header[x], nbchar = input$nbcharPreview)
#         }), collapse = "<br>"))
#       }
#     } else {
#       return("") 
#     }
#   })
  
  output$bt_doMerge <- downloadHandler(
    filename = function() {
      "MergedFasta.fasta"
    },
    content = function(file) {
      f1 <- seqinr::read.fasta(inFile1()$datapath, seqtype = input$seqtype,
                               as.string = TRUE)
      #f1 <- seqinr::read.fasta("WF1.fasta", seqtype = "AA",
                               #as.string = TRUE)
      
      f2 <- seqinr::read.fasta(inFile2()$datapath, seqtype = input$seqtype,
                               as.string = TRUE)
      #f2 <- seqinr::read.fasta("WF2.fasta", seqtype = "AA",
                               #as.string = TRUE)
      
      headers <- clearHeader(seqinr::getAnnot(c(f1, f2)))
      
      dups <- if (input$cb_RemoveRedundancies) {
        duplicated(c(as.character(f1), as.character(f2)))
      } else {
        FALSE
      }
      
      revs <- if (input$cb_RemoveReverses) {
        grepl("REVERSED", headers)
      } else {
        FALSE
      }
      
      conts <- if (input$cb_RemoveContaminants) {
        grepl("CONTAMINANT", headers)
      } else {
        FALSE
      }
      
      tagsr <- if (input$cb_RemoveTag) {
        grepl(input$tx_tagRemove, headers, fixed = TRUE)
      } else {
        FALSE
      }
      
      unqs <- c(f1, f2)[!(dups | revs | conts | tagsr)]
      
      seqinr::write.fasta(sequences = unqs, 
                          names = clearHeader(seqinr::getAnnot(unqs)),
                          file.out = file,
                          as.string = TRUE)
    }
  )
  
})

clearHeader <- function(x) {
  gsub("^>", "", unlist(x))
}

#Adaptado de write.fasta do pacote seqinr
write.oneseq <- function(pep, name, nbchar = 60) {
  start <- paste(name, "<br>", sep = "")
  l <- length(pep)
  q <- floor(l/nbchar)
  r <- l - nbchar * q
  aas <- ""
  if (q > 0) {
    for (x in seq_len(q)) {
      aas <- paste0(aas, seqinr::c2s(pep[(nbchar * (x - 1) + 1):(nbchar * x)]), "<br>", sep = "")
    }
  }
  if (r > 0) {
    aas <- paste0(aas, seqinr::c2s(pep[(nbchar * q + 1):l]), "<br>", sep = "")
  }
  paste0(start, aas)
}

filter.fasta <- function(sequences, minRes, maxRes) {
  nRes <- sapply(fastaFile, length)
  filtered <- sequences
  minResTF <- if(minRes > 0) {
    nRes >= minRes
  } else TRUE
  maxResTF <- if(maxRes > 0) {
    nRes <= maxRes
  } else TRUE
  filtered <- sequences[minResTF & maxResTF]
  return(filtered)
}





