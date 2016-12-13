library(shiny)
#library(seqinr)

# Define server logic required to draw a histogram
shinyServer(function(input, output) {
  
  ##### FASTA TOOLS #####
  
  #inFile1 <- reactive(input$file_Fasta1)
  #inFile2 <- reactive(input$file_Fasta2)
  
  Fasta1 <- reactive(seqinr::read.fasta(input$file_Fasta1$datapath, seqtype = input$seqtype,
                                        as.string = TRUE))
  
  output$hp_NumOfSeq <- renderUI({
    
    if (is.null(input$file_Fasta1)) {
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
      
      f2 <- seqinr::read.fasta(input$file_Fasta2$datapath, seqtype = input$seqtype,
                               as.string = TRUE)
      #f2 <- seqinr::read.fasta("WF2.fasta", seqtype = "AA",
      #as.string = TRUE)
      
      headers <- clearHeader(seqinr::getAnnot(c(f1, f2)))
      
      dups <- if (input$cb_rmFastaDups) {
        duplicated(c(as.character(f1), as.character(f2)))
      } else {
        FALSE
      }
      
      revs <- if (input$cb_rmFastaRevs) {
        grepl("REVERSED", headers)
      } else {
        FALSE
      }
      
      conts <- if (input$cb_rmFastaConts) {
        grepl("CONTAMINANT", headers)
      } else {
        FALSE
      }
      
      tagsr <- if (input$cb_rmFastaTag && (input$tx_tagRemove != "")) {
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
      fastaIn <- Fasta1()#seqinr::read.fasta(inFile1()$datapath, seqtype = input$seqtype,
                         #      as.string = TRUE)
      
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
  
  
  
  ##### GO TOOLS #####
  
  b2gTable <- reactive(switch(input$cb_GOfileType,
                              tsv = readr::read_tsv(input$file_b2gTable$datapath),
                              csv1 = readr::read_csv(input$file_b2gTable$datapath),
                              csv2 = readr::read_csv2(input$file_b2gTable$datapath)))
  
  b2gTableRef <- reactive(switch(input$cb_GOfileTypeRef,
                              tsv = readr::read_tsv(input$file_b2gTableRef$datapath),
                              csv1 = readr::read_csv(input$file_b2gTableRef$datapath),
                              csv2 = readr::read_csv2(input$file_b2gTableRef$datapath)))
  
  output$hp_NumGOLines <- renderUI({
    if (is.null(input$file_b2gTable)) {
      helpText("Number of lines:", style = "font-weight: bold; color: #000000;")
    } else {
      tab <- b2gTable()
      helpText(paste("Number of lines:", nrow(tab)), style = "font-weight: bold; color: #000000;")
    }
  })
  
  output$bt_extractGO_all <-  downloadHandler(
    filename = function() {
      "GO_all_list.txt"
    },
    content = function(file) {
      out <- extract.go(GOvec = b2gTable()[[input$tx_GOcol]],
                        type = "all", removeCat = input$cb_rmGOcat,
                        removeDups = input$cb_rmGOdups,
                        GOsep = input$tx_GOsep)$vec
      
      writeLines(out, con = file)
    })
  
  output$bt_extractGO_C <-  downloadHandler(
    filename = function() {
      "GO_C_list.txt"
    },
    content = function(file) {
      out <- extract.go(GOvec = b2gTable()[[input$tx_GOcol]],
                        type = "C", removeCat = input$cb_rmGOcat,
                        removeDups = input$cb_rmGOdups,
                        GOsep = input$tx_GOsep)$vec
      
      writeLines(out, con = file)
    })
  
  output$bt_extractGO_F <-  downloadHandler(
    filename = function() {
      "GO_F_list.txt"
    },
    content = function(file) {
      out <- extract.go(GOvec = b2gTable()[[input$tx_GOcol]],
                        type = "F", removeCat = input$cb_rmGOcat,
                        removeDups = input$cb_rmGOdups,
                        GOsep = input$tx_GOsep)$vec
      
      writeLines(out, con = file)
    })
  
  output$bt_extractGO_P <-  downloadHandler(
    filename = function() {
      "GO_P_list.txt"
    },
    content = function(file) {
      out <- extract.go(GOvec = b2gTable()[[input$tx_GOcol]],
                        type = "P", removeCat = input$cb_rmGOcat,
                        removeDups = input$cb_rmGOdups,
                        GOsep = input$tx_GOsep)$vec
      
      writeLines(out, con = file)
    })
  
  
  FisherOut <- eventReactive(input$bt_GOFischer, {
    GOvec <- extract.go(GOvec = b2gTable()[[input$tx_GOcol]],
                                type = "all", removeCat = FALSE,
                                removeDups = FALSE,
                                GOsep = input$tx_GOsep)$vec
      
    GOvecRef <- extract.go(GOvec = b2gTable()[[input$tx_GOcolRef]],
                           type = "all", removeCat = FALSE,
                           removeDups = FALSE,
                           GOsep = input$tx_GOsep)$vec
    
    GOvec <- GOvec[GOvec %in% unique(GOvecRef)]
    
    Fisher <- lapply(seq_along(unique(GOvec)), function(GO) {
      target <- GOvec[GO][1]
      
      target_exp <- sum(GOvec == target)
      target_ref <- sum(GOvecRef == target)
      
      rest_exp <- length(GOvec) - target_exp
      rest_ref <- length(GOvecRef) - target_ref
      
      mat <- matrix(c(target_exp, target_ref, rest_exp, rest_ref), ncol = 2)
      
      f <- fisher.test(mat, alternative = "greater")
      
      data.frame(ID = target, p.value = f$p.value)
    })
    
    do.call(rbind, Fisher)
  })
  
  
  output$tb_GOstat <- DT::renderDataTable(FisherOut(),
                                          rownames= FALSE,
                                          options = list(
                                            autoWidth = TRUE,
                                            columnDefs = list(list(width = '50%', targets = "_all"))))
  
  output$bt_writeGOstat <-  downloadHandler(
    filename = function() {
      "GO_Enrichment_Test.csv"
    },
    content = function(file) {
      if (input$cb_GOfileType == "csv1") {
        write.csv(x = FisherOut(), file = file, row.names = FALSE)
      } else {
        write.csv2(x = FisherOut(), file = file, row.names = FALSE)
      }
    })
  
  output$plot_GOpie <- renderPlot({
    print(input$file_b2gTable)
    if (!is.null(input$file_b2gTable)) {
      cats <- extract.go(GOvec = b2gTable()[[input$tx_GOcol]],
                               type = "all", removeCat = input$cb_rmGOcat,
                               removeDups = input$cb_rmGOdups,
                       GOsep = input$tx_GOsep)$cats
    
    pie(table(cats))
    } else {
      NULL
    }
  })
  
  
  
  
  
  
  
  
  ##### iTraq Ratio Check #####
  
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

extract.go <- function(GOvec, type = c("all", "C", "F", "P"),
                       removeCat = TRUE, removeDups = FALSE,
                       GOsep) {
  vec <- unlist(strsplit(na.omit(GOvec), split = GOsep))
  
  if (removeDups) {
    vec <- unique(vec)
  }
  
  cats <- gsub(":GO:\\d{7}", "", vec)
  
  if (removeCat) {
    vec <- gsub("[CFP]:", "", vec)
  }
  
  if (type[1] != "all") {
    vec <- vec[cats == type[1]]
  }
  
  list(vec = vec, cats = cats)
}




