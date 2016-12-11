library(shiny)
#library(seqinr)

# Define server logic required to draw a histogram
shinyServer(function(input, output) {
  
  inFile1 <- reactive(input$file1)
  
  inFile2 <- reactive(input$file2)
  
  Fasta1 <- reactive(seqinr::read.fasta(input$file1$datapath, seqtype = input$seqtype,
                                        as.string = TRUE))
  
  output$hp_NumOfSeq <- renderUI({
    
    if (is.null(inFile1())) {
      helpText("Number of sequences:", style = "font-weight: bold; color: #000000;")
    } else {
      sequences <- Fasta1()
      helpText(paste("Number of sequences:", length(sequences)), style = "font-weight: bold; color: #000000;")
    }
  })
  
  output$bt_doMerge <- downloadHandler(
    filename = function() {
      "MergedFasta.fasta"
    },
    content = function(file) {
      f1 <- Fasta1() #seqinr::read.fasta(inFile1()$datapath, seqtype = input$seqtype,
                     #             as.string = TRUE)
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
      
      tagsr <- if (input$cb_RemoveTag && (input$tx_tagRemove != "")) {
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
  
  output$bt_doFilter <- downloadHandler(
    filename = function () {
      "FilteredFasta.fasta"
    },
    content = function(file) {
      fastaIn <- seqinr::read.fasta(inFile1()$datapath, seqtype = input$seqtype,
                               as.string = TRUE)
      
      # fastaIn <- seqinr::read.fasta("2seqFasta.fasta", seqtype = "AA",
      #                               as.string = TRUE)
      
      out <- filter.fasta(sequences = fastaIn,
                          headers = input$fileHeaderFilter,
                          minRes = input$minSeqFilter, maxRes = input$maxSeqFilter)
      
      seqinr::write.fasta(sequences = out, 
                          names = clearHeader(seqinr::getAnnot(out)),
                          file.out = file,
                          as.string = TRUE)
    }
  )
  
  observeEvent(input$bt_resetHeaderList, {
    reset("fileHeaderFilter")
  })
  
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

filter.fasta <- function(sequences, headers, minRes, maxRes) {
  nRes <- sapply(sequences, nchar)
  minResTF <- if(!is.na(minRes) && minRes > 0) {
    nRes >= minRes
  } else TRUE
  maxResTF <- if(!is.na(maxRes) && maxRes > 0) {
    nRes <= maxRes
  } else TRUE
  headersTF <- if(!is.null(headers)) {
    headersText <- readLines(headers$datapath)
    names(sequences) %in% headersText
  } else TRUE
  filtered <- sequences[minResTF & maxResTF & headersTF]
  return(filtered)
}





