library(shiny)
#library(seqinr)

# Define server logic required to draw a histogram
shinyServer(function(input, output, session) {
  
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
  
  output$vb_GOhelp <- renderUI(HTML(paste("1. Select the file type of the GO IDs table to be analyzed. 
The file must be a table with column names, even if it has only one column.",
                                       "2. Select and upload the file.",
                                       "3. From the dropdown list, select the column with the GO IDs.",
                                       "4. In case there are multiple IDs per line (say multiple IDs related to one protein), 
the app will try to identify the string between IDs. If the guess isn't correct, 
type in the correct separator.",
                                       sep = "<BR><BR>")))
  
  observe(shinyjs::toggleState("file_b2gTable", condition = input$cb_GOfileType != "wait"))
  observe(shinyjs::toggleState("cb_GOcol", condition = !is.null(input$file_b2gTable$datapath)))
  
  b2gTable <- eventReactive(input$file_b2gTable, {
    out <- switch(input$cb_GOfileType,
                  tsv = readr::read_tsv(input$file_b2gTable$datapath),
                  csv1 = readr::read_csv(input$file_b2gTable$datapath),
                  csv2 = readr::read_csv2(input$file_b2gTable$datapath))
    
    updateSelectInput(session, "cb_GOcol", choices = c("Select one", colnames(out)))
    return(out)
  })
  
  # These output and output options make sure file upload is tracked and the call to the eventReactive above is executed when necessary
  output$file_GOUploaded <- reactive({
    return(!is.null(b2gTable()))
  })
  outputOptions(output, "file_GOUploaded", suspendWhenHidden = FALSE)
  
  observeEvent(input$cb_GOcol, {
    GOs <- as.character(b2gTable()[[input$cb_GOcol]])
    cand <- GOs[nchar(GOs, keepNA = FALSE) > 12][1]
    sep <- regmatches(cand, regexpr("(?<=\\d{7})([^CFPGO]*)", cand, perl = TRUE))
    updateTextInput(session, "tx_GOsep", value = sep)
  })
  
  
  
  ### Descriptive 
  
  observeEvent(input$bt_GOdesc, {
    
    if (!is.null(input$file_b2gTable) && input$cb_GOcol != "wait") {
      GOvec <- b2gTable()[[input$cb_GOcol]]
      
      extracted <- extract.go(GOvec = GOvec,
                              type = "all", removeCat = FALSE,
                              removeDups = FALSE,
                              GOsep = input$tx_GOsep)
      
      cats <- extracted$cats 
      
      if (length(unique(cats)) != 3) {
        hide("plot_GOpie")
        hide("df_GOcounts")
        hide("dt_GOall")
        return (NULL)
      }
      vec <- extracted$vec
      
      
      cats_names <- c(F = "Molecular\nFunction", P = "Biological\nProcess", C = "Cellular\nComponent")
      tab_cats <- table(cats)
      tab_ids <- table(vec)
      
      vec_C <- vec[cats == "C"]
      vec_F <- vec[cats == "F"]
      vec_P <- vec[cats == "P"]
      
      output$plot_GOpie <- renderPlot({
        par(mar = c(4,3,2,1))
        bp <- barplot(tab_cats, names.arg = cats_names[names(tab_cats)],
                      col = "lightblue", xlab = "Category", ylab = "Counts")
        text(x = bp, y = tab_cats, adj = c(0.5, 1), labels = tab_cats)
      })
      
      counts_table <- data.frame(Type = c("Cellular component", "Molecular function", 
                                          "Biological process", "All"),
                                 Total = c(tab_cats["C"], tab_cats["F"], 
                                           tab_cats["P"], length(vec)),
                                 Unique = c(length(unique(vec_C)), length(unique(vec_F)),
                                            length(unique(vec_P)), length(unique(vec))))
      
      output$df_GOcounts <- renderTable(counts_table)
      
      GO_counted <- as.data.frame(tab_ids)
      colnames(GO_counted) <- c("ID", "Counts")
      GO_counted$Category <- cats_names[gsub(":GO:\\d{7}", "", GO_counted$ID)]
      GO_counted <- GO_counted[order(GO_counted$Counts, GO_counted$ID, decreasing = TRUE),]
      
      output$dt_GOall <- DT::renderDataTable(GO_counted, 
                                             rownames = FALSE,
                                             class = "compact",
                                             options = list(
                                               autoWidth = TRUE,
                                               columnDefs = list(list(width = "80%", targets = 2),
                                                                 list(className = 'dt-center', targets = "_all"))))
      
    } else {
      NULL
    }
    
  })
  
  
  #### Extract IDs
  
  observe(shinyjs::toggleState("file_b2gTableRef", condition = input$cb_GOfileTypeRef != "wait"))
  observe(shinyjs::toggleState("cb_GOcolRef", condition = !is.null(input$file_b2gTableRef$datapath)))
  
  b2gTableRef <- eventReactive(input$file_b2gTableRef, {
    out <- switch(input$cb_GOfileTypeRef,
                  tsv = readr::read_tsv(input$file_b2gTableRef$datapath),
                  csv1 = readr::read_csv(input$file_b2gTableRef$datapath),
                  csv2 = readr::read_csv2(input$file_b2gTableRef$datapath))
    
    updateSelectInput(session, "cb_GOcolRef", choices = c("Select one", colnames(out)))
    return(out)
  })
  
  # These output and output options make sure file upload is tracked and the call to the eventReactive above is executed when necessary
  output$file_GOUploadedRef <- reactive({
    return(!is.null(b2gTableRef()))
  })
  outputOptions(output, "file_GOUploadedRef", suspendWhenHidden = FALSE)
  
  observeEvent(input$cb_GOcolRef, {
    GOs <- as.character(b2gTableRef()[[input$cb_GOcolRef]])
    cand <- GOs[nchar(GOs, keepNA = FALSE) > 12][1]
    sep <- regmatches(cand, regexpr("(?<=\\d{7})([^CFPGO]*)", cand, perl = TRUE))
    updateTextInput(session, "tx_GOsepRef", value = sep)
  })
  
  output$hp_NumGOLines <- renderUI({
    if (is.null(input$file_b2gTable)) {
      helpText("Number of lines:", style = "font-weight: bold; color: #000000; visibility: hidden;")
    } else {
      tab <- b2gTable()
      helpText(paste("Number of lines:", nrow(tab)), style = "font-weight: bold; color: #000000;")
    }
  })
  
  output$hp_NumGOLinesRef <- renderUI({
    if (is.null(input$file_b2gTableRef)) {
      helpText(" ", style = "font-weight: bold; color: #000000; visibility: hidden;")
    } else {
      tab <- b2gTableRef()
      helpText(paste("Number of lines:", nrow(tab)), style = "font-weight: bold; color: #000000;")
    }
  })
  
  output$bt_extractGO_all <-  downloadHandler(
    filename = function() {
      "GO_all_list.txt"
    },
    content = function(file) {
      out <- extract.go(GOvec = b2gTable()[[input$cb_GOcol]],
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
      out <- extract.go(GOvec = b2gTable()[[input$cb_GOcol]],
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
      out <- extract.go(GOvec = b2gTable()[[input$cb_GOcol]],
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
      out <- extract.go(GOvec = b2gTable()[[input$cb_GOcol]],
                        type = "P", removeCat = input$cb_rmGOcat,
                        removeDups = input$cb_rmGOdups,
                        GOsep = input$tx_GOsep)$vec
      
      writeLines(out, con = file)
    })
  
  
  #### Enrichment
  
  FisherOut <- eventReactive(input$bt_GOFisher, {
    GOvec <- extract.go(GOvec = b2gTable()[[input$cb_GOcol]],
                        type = "all", removeCat = FALSE,
                        removeDups = FALSE,
                        GOsep = input$tx_GOsep)$vec
    
    GOvecRef <- extract.go(GOvec = b2gTableRef()[[input$cb_GOcolRef]],
                           type = "all", removeCat = FALSE,
                           removeDups = FALSE,
                           GOsep = input$tx_GOsepRef)$vec
    
    GOvec <- GOvec[GOvec %in% unique(GOvecRef)]
    
    if (is.null(GOvec)) return(NULL)
    
    # progress <- shiny::Progress$new()
    
    # on.exit(progress$close())
    # https://shiny.rstudio.com/articles/progress.html
    # progress$set(message = "Performing tests", value = 0)
    
    nGOs <- length(unique(GOvec))
    
    prog.interval <- seq(0, nGOs, length.out = 100)
    
    withProgress(message = "Performing tests", value = 0,
                 {
                   Fisher <- lapply(seq_along(unique(GOvec)), function(GO) {
                     target <- GOvec[GO][1]
                     
                     target_exp <- sum(GOvec == target)
                     target_ref <- sum(GOvecRef == target)
                     
                     rest_exp <- length(GOvec) - target_exp
                     rest_ref <- length(GOvecRef) - target_ref
                     
                     mat <- matrix(c(target_exp, target_ref, rest_exp, rest_ref), ncol = 2)
                     
                     f <- fisher.test(mat, alternative = "greater")
                     
                     setProgress(findInterval(GO, prog.interval)/100,
                                 detail = paste(GO, "out of", nGOs))
                     #incProgress(1/nGOs, detail = paste(GO, "out of", nGOs))
                     
                     data.frame(ID = target, p.value = f$p.value)
                   })
                 })
    #do.call(rbind, Fisher)
    out <- do.call(rbind, Fisher)
    out$q.value <- p.adjust(out$p.value, method = "fdr")
    out
    
  })
  
  output$tb_GOstat <- DT::renderDataTable({
    fisher <- FisherOut()
    if (!is.null(fisher)) {
      dt <- DT::datatable(fisher, rownames = FALSE, 
                          options = list(autoWidth = TRUE,
                                         columnDefs = list(list(width = '30%', targets = 0))))
      DT::formatSignif(dt, c("p.value", "q.value"), 5)
    } else {
      NULL
    }
  })
  
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
  
  ##### Enrichment Tester #####
  
  # Main file UI behaviour 
  
  observe(shinyjs::toggleState("file_ET", condition = input$se_ETfileType != "wait"))
  observe(shinyjs::toggleState("se_ETcol", condition = !is.null(input$file_ET$datapath)))
  
  ETmainfile <- eventReactive(input$file_ET, {
    out <- switch(input$se_ETfileType,
                  tsv = readr::read_tsv(input$file_ET$datapath),
                  csv1 = readr::read_csv(input$file_ET$datapath),
                  csv2 = readr::read_csv2(input$file_ET$datapath))
    
    updateSelectInput(session, "se_ETcol", choices = c("Select one", colnames(out)))
    return(out)
  })
  
  # These output and output options make sure file upload is tracked and the call to the eventReactive above is executed when necessary
  output$file_ETuploaded <- reactive({
    return(!is.null(ETmainfile()))
  })
  outputOptions(output, "file_ETuploaded", suspendWhenHidden = FALSE)
  
  observe(shinyjs::toggleState("tx_ETsep", condition = input$rb_ETsep == "mult" && !is.null(input$file_ET$datapath)))
  
  observe({
    #shinyjs::toggleState("se_ETcountCol", condition = input$rb_ETcountedData == "nocount")
    if (input$rb_ETcountedData == "nocount" && !is.null(input$file_ET$datapath)) {
      shinyjs::enable("se_ETcountCol")
      numCols <- colnames(ETmainfile())[sapply(ETmainfile(), class) == "integer"]
      updateSelectInput(session, "se_ETcountCol",
                        choices = numCols)
    } else {
      shinyjs::disable("se_ETcountCol")
      updateSelectInput(session, "se_ETcountCol",
                        choices = "")
    }
  })

  output$hp_NumETLines <- renderUI({
    if (is.null(input$file_ET)) {
      helpText("Number of lines:", style = "font-weight: bold; color: #000000; visibility: hidden;")
    } else {
      tab <- ETmainfile()
      helpText(paste("Number of lines:", nrow(tab)), style = "font-weight: bold; color: #000000;")
    }
  })
  
  # Ref file UI behaviour 
  
  observe({
    if (input$rb_ETrefStyle == "same") {
      updateSelectInput(session, "se_ETcolRef", choices = c("Select one", colnames(ETmainfile())))
    } 
  })
  
  observe(shinyjs::toggleState("file_ETRef", condition = input$se_ETfileTypeRef != "wait"))
  observe(shinyjs::toggleState("se_ETcolRef", 
                                condition = !((input$rb_ETrefStyle == "other" && is.null(input$file_ETRef$datapath)) ||
                                 (input$rb_ETrefStyle == "same" && is.null(input$file_ET$datapath)))))
  
  
  ETReffile <- reactive({
    if (input$rb_ETrefStyle == "other" && !is.null(input$file_ETRef$datapath)) {
      out <- switch(input$se_ETfileTypeRef,
                  tsv = readr::read_tsv(input$file_ETRef$datapath),
                  csv1 = readr::read_csv(input$file_ETRef$datapath),
                  csv2 = readr::read_csv2(input$file_ETRef$datapath))
    
    updateSelectInput(session, "se_ETcolRef", choices = c("Select one", colnames(out)))
    return(out)
    } else if (input$rb_ETrefStyle == "same") {
      ETmainfile()
    } else {
      NULL
    }
  })
  
  # These output and output options make sure file upload is tracked and the call to the eventReactive above is executed when necessary
  output$file_ETuploadedRef <- reactive({
    return(!is.null(ETReffile()))
  })
  outputOptions(output, "file_ETuploadedRef", suspendWhenHidden = FALSE)
  
  observe(shinyjs::toggleState("tx_ETsepRef", condition = input$rb_ETsepRef == "mult"))
  observe({
    #shinyjs::toggleState("se_ETcountColRef", condition = input$rb_ETcountedDataRef == "nocount")
    dataFile <- if (input$rb_ETrefStyle == "same") ETmainfile() else ETReffile()
    if (input$rb_ETcountedDataRef == "nocount") {
      shinyjs::enable("se_ETcountColRef")
      numCols <- colnames(dataFile)[sapply(dataFile, class) == "integer"]
      updateSelectInput(session, "se_ETcountColRef",
                        choices = numCols)
    } else {
      shinyjs::disable("se_ETcountColRef")
      updateSelectInput(session, "se_ETcountColRef",
                        choices = "")
    }
  })
  
  # Not working
  output$hp_NumGOLinesRef <- renderUI({
    if (is.null(input$file_ETRef)) {
      helpText(" ", style = "font-weight: bold; color: #000000; visibility: hidden;")
    } else if (input$rb_ETrefStyle == "other") {
      tab <- ETReffile()
      helpText(paste("Number of lines:", nrow(tab)), style = "font-weight: bold; color: #000000;")
    } else {
      helpText(" ", style = "font-weight: bold; color: #000000; visibility: hidden;")
    }
  })
  
  ### Analysis code
  
  ETFisherOut <- eventReactive(input$bt_ETFisher, {
    
    NAs <- is.na(ETmainfile()[[input$se_ETcol]]) | is.na(ETmainfile()[[input$se_ETcountCol]])
    
    NAsRef <- is.na(ETReffile()[[input$se_ETcolRef]]) | is.na(ETReffile()[[input$se_ETcountColRef]])
    
    vec <- as.character(ETmainfile()[[input$se_ETcol]][!NAs])
    
    if (input$rb_ETsep == "mult") {
      vec <- unlist(strsplit(vec, split = input$tx_GOsep))
    } 
    
    vecRef <- if (input$rb_ETrefStyle == "same") {
      as.character(ETmainfile()[[input$se_ETcolRef]][!NAsRef])
    } else {
      as.character(ETReffile()[[input$se_ETcolRef]][!NAsRef])
    }
    
    if (input$rb_ETsepRef == "mult") {
      vecRef <- unlist(strsplit(vecRef, split = input$tx_GOsepRef))
    } 
    
    # Keep only the main IDs that appear on Reference
    foundOnRef <- vec %in% unique(vecRef) 
    vec <- vec[foundOnRef]
    
    foundOnRef_Ref <- vecRef %in% unique(vec)
    vecRef <- vecRef[foundOnRef_Ref]

    if (is.null(vec)) {
      message("Returning NULL")
      return(NULL)
    }

    nIDs <- length(unique(vec))
    
    prog.interval <- seq(0, nIDs, length.out = 100)
    
    tab <- if (input$rb_ETcountedData == "count") {
      unlist(table(vec))
    } else {
      dups <- duplicated(vec)
      setNames(ETmainfile()[[input$se_ETcountCol]][!NAs][foundOnRef],
               vec)[!dups]
    }
    
    tabRef <- if (input$rb_ETcountedDataRef == "count") {
      unlist(table(vecRef))
    } else {
      dups <- duplicated(vecRef)
      setNames(ETReffile()[[input$se_ETcountColRef]][!NAsRef][foundOnRef_Ref],
               vecRef)[!dups]
    }

    withProgress(message = "Performing tests", value = 0,
                 {
                   Fisher <- lapply(seq_along(unique(vec)), function(ID) {
                     target <- vec[ID]
                     
                     target_exp <- tab[target]
                     target_ref <- tabRef[target]
                     
                     rest_exp <- sum(tab) - target_exp
                     rest_ref <- sum(tabRef) - target_ref
                     
                     mat <- matrix(c(target_exp, target_ref, rest_exp, rest_ref), ncol = 2)
                     
                     f <- fisher.test(mat, alternative = "greater")
                     
                     setProgress(findInterval(ID, prog.interval)/100,
                                 detail = paste(ID, "out of", nIDs))
                     
                     data.frame(ID = target, p.value = f$p.value)
                   })
                 })
    
    out <- do.call(rbind, Fisher)
    out$q.value <- p.adjust(out$p.value, method = "fdr")
    out$TargetCount <- tab
    out$ReferenceCount <- tabRef
    out$TargetRatio <- out$TargetCount/sum(out$TargetCount)
    out$ReferenceRatio <- out$ReferenceCount/sum(out$ReferenceCount)
    out$ExpectedByChance <- out$ReferenceRatio*sum(out$TargetCount)
    out$OverRepresentation <- out$TargetCount/out$ExpectedByChance
    out
    
  })
  
  output$tb_ETstat <- DT::renderDataTable({
    fisher <- ETFisherOut()
    if (!is.null(fisher)) {
      dt <- DT::datatable(fisher, rownames = FALSE,
                          options = list(autoWidth = TRUE,
                                         columnDefs = list(list(width = '30%', targets = 0))))
      DT::formatSignif(dt, c("p.value", "q.value", "TargetRatio", "ReferenceRatio", 
                             "ExpectedByChance", "OverRepresentation"), 3)
    } else {
      NULL
    }
  })
  
  output$bt_writeETstat <-  downloadHandler(
    filename = function() {
      "Enrichment_Test.csv"
    },
    content = function(file) {
      if (input$se_ETfileType == "csv1") {
        write.csv(x = ETFisherOut(), file = file, row.names = FALSE)
      } else {
        write.csv2(x = ETFisherOut(), file = file, row.names = FALSE)
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
  minResTF <- if (!is.na(minRes) && minRes > 0) {
    nRes >= minRes
  } else TRUE
  maxResTF <- if (!is.na(maxRes) && maxRes > 0) {
    nRes <= maxRes
  } else TRUE
  headersTF <- if (!is.null(headers)) {
    headersText <- readLines(headers$datapath)
    names(sequences) %in% headersText
  } else TRUE
  filtered <- sequences[minResTF & maxResTF & headersTF]
  return(filtered)
}

extract.go <- function(GOvec, type = c("all", "C", "F", "P"),
                       removeCat = TRUE, removeDups = FALSE,
                       GOsep) {
  if (GOsep != "") {
    vec <- unlist(strsplit(as.character(na.omit(GOvec)), split = GOsep))
  } else {
    vec <- as.character(na.omit(GOvec))
  }
  
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




