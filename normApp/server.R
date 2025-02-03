library(flowDensity)
library(flowCore)
library(shiny)
library(FlowCIPHE)
library(flowPeaks)
library(openCyto)
library(colorspace)
library(ggplot2)
library(RColorBrewer)
library(parallel)
library(car)
library(mclust)
library(future.apply)

# Increase size limit of files for uploading  
options(expressions = 5e5, shiny.maxRequestSize = 100 * 1024^3)
# source("functions.R")
function(input, output, session) {
  gated_result <- reactiveVal()
  # Initialize objects
  listObject <- reactiveValues(
    listFCS = NULL,
    flow.frames = NULL,
    CSTFile=NULL,
    listCST=NULL,
    markerFCS=NULL,
    flow.frames.CST=NULL,
    flow.frames.PCA=NULL,
    flow.frames.clusters=NULL,
    data.table=NULL,
    csvFileTransfo=NULL,
    beads=NULL,
    flow.frames.beads=NULL,
    flow.frames.cells=NULL,
    gate=NULL,
    bead_gated=NULL,
    flow.frames.bead_peak=NULL,
    groups=NULL,
    concat=NULL,
    mapping_table=NULL,
    markerCurrentGroup=NULL,
    markerRefGroup=NULL)
  
  listObjects <- reactiveValues()
  
  
  ################################################################################
  ##################### UPLOAD DATA ##############################################
  
 
  output$file_inputs <- renderUI({
    req(input$selector)
    
    lapply(seq_len(input$selector), \(i) {
      fileInput(
        inputId = paste0("file_input_", i),
        label = paste("Upload group:", i),
        multiple=TRUE
      )
    })
  })
  
  
  output$beads_inputs<- renderUI({
    req(input$selector)
    
    lapply(seq_len(input$selector), \(i) {
      fileInput(
        inputId = paste0("beads_input_", i),
        label = paste("(Optional) Upload beads already extracted for group:", i),
        multiple=FALSE
      )
    })
  })
  
  
  # Upload cells (+beads) files
  plan(multisession, workers = parallel::detectCores() - 1)
  print(detectCores() - 1)
  observe({

    lapply(seq_len(input$selector), function(i) {
      observeEvent(input[[paste0("file_input_", i)]], {
        tryCatch({
          progress <- Progress$new()
          progress$set(message = paste("Processing group", i), value = 0.5)
          
          files <- input[[paste0("file_input_", i)]]
          newfilesnames <- files$datapath
          original_filenames <- basename(files$name)
          
          if (is.null(listObjects[[paste0("group_", i)]])) {
            listObjects[[paste0("group_", i)]] <- list(
              listFCS = NULL,
              flow.frames = NULL,
              markerFCS = NULL,
              data.table = NULL
            )
          }
          
          currentGroup <- listObjects[[paste0("group_", i)]]
          currentGroup$original_filenames <- c(currentGroup$original_filenames, original_filenames) 
          currentGroup$listFCS <- c(currentGroup$listFCS, newfilesnames)
          
          validated_files <- future_lapply(currentGroup$listFCS, function(file) {
            flow_data <- read.FCS(file, truncate_max_range = FALSE, ignore.text.offset = TRUE)
            if (!is.matrix(exprs(flow_data))) stop()
            flow_data
          })
          
          concatenatedFCS <- if (length(validated_files) == 1) {
            validated_files[[1]]
          } else {
            concatenate.FCS.CIPHE(validated_files,"FlagNormApp")
          }
          
          currentGroup$markerFCS <- concatenatedFCS@parameters@data$name
          names(currentGroup$markerFCS) <- NULL
          currentGroup$flow.frames <- concatenatedFCS
          
          currentGroup$data.table <- data.frame(Fluo = currentGroup$markerFCS, Arg = rep("none", length(currentGroup$markerFCS)))
          rownames(currentGroup$data.table) <- currentGroup$data.table$Fluo
          currentGroup$data.table$Fluo <- NULL
          
          listObjects[[paste0("group_", i)]] <- currentGroup
          
          updateSelectizeInput(session, "groupTransfo", choices = names(listObjects))
          updateSelectizeInput(session, "groupQC", choices = names(listObjects))
          updateSelectizeInput(session, "referenceGroup", choices = names(listObjects))
          
          progress$close()
          
        }, error = function(e) {
          progress$close()
        })
      })
    })
  })
  
  # Upload beads files
  observe({
    lapply(seq_len(input$selector), \(i) {
      observeEvent(input[[paste0("beads_input_", i)]], {
        
        files <- input[[paste0("beads_input_", i)]]
        
        newfilesnames <- files$datapath
        
        
        # Initialize the listObject for this group if not already done
        if (is.null(listObjects[[paste0("group_", i)]])) {
          listObjects[[paste0("group_", i)]] <- list(
            listFCS = NULL,
            flow.frames = NULL,
            markerFCS = NULL,
            data.table = NULL
          )
        }
        
        # Access the listObject for the current group
        currentGroup <- listObjects[[paste0("group_", i)]]
        flow_data <- read.FCS(newfilesnames, truncate_max_range = FALSE, ignore.text.offset = TRUE)
        
        currentGroup$flow.frames.beads<-flow_data
        
        
        if (is.null(currentGroup$markerFCS)){
          currentGroup$markerFCS<-flow_data@parameters@data$name
          names(currentGroup$markerFCS) <- NULL
        }
        # Update the group in listObjects
        listObjects[[paste0("group_", i)]] <- currentGroup
        
        updateSelectizeInput(session, "groupQC", choices = names(listObjects))
        
        updateSelectizeInput(session, "referenceGroup", choices = names(listObjects))
      })
    })
  })
  
  
  
  # Initialize markers 
  observeEvent(input$groupTransfo,{
    req(input$groupTransfo)  
    if (!is.null(listObjects)) {
      if (!is.null(listObjects[[input$groupTransfo]])) {
        currentGroup <- listObjects[[input$groupTransfo]]
        
        # Update UI elements for the group
        if ("FSC-A" %in% currentGroup$markerFCS && "SSC-A" %in% currentGroup$markerFCS) {
          updateSelectInput(session, "marker_x",  choices = currentGroup$markerFCS, selected = "FSC-A")
          updateSelectInput(session, "marker_y",  choices = currentGroup$markerFCS, selected = "SSC-A")
          updateSelectInput(session, "marker_x2", choices = currentGroup$markerFCS, selected = "FSC-A")
          updateSelectInput(session, "marker_y2", choices = currentGroup$markerFCS, selected = "SSC-A")
          
        } else {
          updateSelectInput(session, "marker_x",  choices = currentGroup$markerFCS)
          updateSelectInput(session, "marker_y",  choices = currentGroup$markerFCS)
          updateSelectInput(session, "marker_x2", choices = currentGroup$markerFCS)
          updateSelectInput(session, "marker_y2", choices = currentGroup$markerFCS)
          
        }
        
      }
      
    }
  })
  
  
  # Display different groups to normalize
  observe({
    
    groups<-names(listObjects)[! names(listObjects) %in% as.vector("group_1")]
    updateSelectizeInput(session, "groupNorm", choices = c("",groups))
    
  })
  
  
  # Display common Markers between ref and group
  observeEvent(input$groupNorm, {
    req(input$groupNorm)
    tryCatch({
   
      currentGroup <- listObjects[[input$groupNorm]]
      refGroup <- listObjects[["group_1"]]
      req(currentGroup$flow.frames)
      
      if (!is.null(currentGroup$flow.frames)) {
        markerCurrentGroup <- colnames(currentGroup$flow.frames)
      } else {
        markerCurrentGroup <- colnames(currentGroup$flow.frames.beads)
      }
      
      if (!is.null(refGroup$flow.frames)) {
        markerRefGroup <- colnames(refGroup$flow.frames)
      } else if (!is.null(refGroup$flow.frames.beads)) {
        markerRefGroup <- colnames(refGroup$flow.frames.beads)
      }
      
      listObject$markerRefGroup <- markerRefGroup
      listObject$markerCurrentGroup <- markerCurrentGroup
      
      commonMarkers <- intersect(markerCurrentGroup, markerRefGroup)
      
      
      if (!is.null(currentGroup$commonMarkers)) {
        selectedMarkers <- currentGroup$commonMarkers  
      } else {
        selectedMarkers <- commonMarkers  
        currentGroup$commonMarkers <- selectedMarkers  
        listObjects[[input$groupNorm]] <- currentGroup
      }
      
      updateSelectizeInput(session, "commonMarkers", 
                           selected = selectedMarkers,   
                           choices = commonMarkers)      
    }, error = function(e) {
      showModal(modalDialog(
        title = "Error",
        paste("An error occurred:", e$message),
        easyClose = TRUE
      ))
    })
  })
  
  # Display common Markers between ref and group
  observeEvent(input$commonMarkers, {
    req(input$groupNorm)
    currentGroup <- listObjects[[input$groupNorm]]
    req(input$commonMarkers)
    
    
    currentGroup$commonMarkers <- input$commonMarkers
    listObjects[[input$groupNorm]] <- currentGroup
  })
  
  
  # with unmatch column
  diff_colnames <- reactive({
    req(listObject$markerRefGroup)
    req(listObject$markerCurrentGroup)
    setdiff(listObject$markerRefGroup,listObject$markerCurrentGroup)  
  })
  

  output$dynamic_selectors <- renderUI({
    req(input$groupNorm)
    req(diff_colnames())
    
    selected_group_name <- input$groupNorm
    marker_unmatch <- setdiff(listObject$markerCurrentGroup, listObject$markerRefGroup)
    
    lapply(seq_along(diff_colnames()), function(i) {
      ref_marker <- diff_colnames()[i]
      
      selected_value <- if (!is.null(listObjects[[selected_group_name]]$mapping) &&
                            !is.null(listObjects[[selected_group_name]]$mapping[[ref_marker]])) {
        listObjects[[selected_group_name]]$mapping[[ref_marker]]
      } else {
        NULL
      }
      
      selectizeInput(paste0("slider_", i), 
                     label = paste("Ref", ref_marker, "match with :"), 
                     choices = c("", marker_unmatch), 
                     selected = selected_value,  
                     width = "150px")
    })
  })
  
  
  # 
  observe({
    req(input$groupNorm)
    req(diff_colnames())
    
    selected_group_name <- input$groupNorm
    
    if (is.null(listObjects[[selected_group_name]]$mapping)) {
      listObjects[[selected_group_name]]$mapping <- list()
    }
    if (is.null(listObjects[["group_1"]]$mapping)) {
      listObjects[["group_1"]]$mapping <- list()
    }
    
    mapping_group <- listObjects[[selected_group_name]]$mapping  
    mapping_ref <- if (!is.null(listObjects[["group_1"]]$mapping[[selected_group_name]])) {
      listObjects[["group_1"]]$mapping[[selected_group_name]]
    } else {
      list()
    }
    
    for (i in seq_along(diff_colnames())) {
      ref_marker <- diff_colnames()[i]
      selected_value <- input[[paste0("slider_", i)]]
      
      if (!is.null(selected_value) && selected_value != "") {
        mapping_group[[ref_marker]] <- selected_value  
        mapping_ref[[selected_value]] <- ref_marker  
      }
    }
    
    listObjects[[selected_group_name]]$mapping <- mapping_group
    listObjects[["group_1"]]$mapping[[selected_group_name]] <- mapping_ref
 
    print(paste("Mapping enregistré pour", selected_group_name, ":"))
    print(as.vector(unlist(listObjects[[selected_group_name]]$mapping)))
    print(names(unlist(listObjects[[selected_group_name]]$mapping)))
    output$result <- renderPrint({
      list(
        group = listObjects[[selected_group_name]]$mapping
      )
    })
  })
  
 
  
  # Render DataTable for the group
  output$table_group_transfo <- DT::renderDataTable({
    req(input$groupTransfo)
    currentGroup <- listObjects[[input$groupTransfo]]
    req(currentGroup$data.table)  
    
    DT::datatable(
      currentGroup$data.table,
      editable = list(target = 'cell', columns = which(names(currentGroup$data.table) == "Arg")),
      selection = "none",
      options = list(
        columnDefs = list(
          list(targets = which(names(currentGroup$data.table) == "Arg"),
               render = JS('function(data, type, row, meta) {
                      if (type === "display") {
                        return parseFloat(data).toFixed(2);
                      }
                      return data;
                    }'))
        )
      )
    )
  })
  
  # update group QC
  observeEvent(input$groupQC,{
    req(input$groupQC)  
    currentGroup <- listObjects[[input$groupQC]]
    updateSelectInput(session, "marker_x4", choices = currentGroup$markerFCS)
    updateSelectInput(session, "marker_y4", choices = currentGroup$markerFCS)
    updateSelectInput(session, "marker_x5", choices = currentGroup$markerFCS)
    updateSelectInput(session, "marker_y5", choices = currentGroup$markerFCS)
    
    if ("FSC-A" %in% currentGroup$markerFCS && "FSC-H" %in% currentGroup$markerFCS) {
      updateSelectInput(session, "marker_x3", choices = currentGroup$markerFCS, selected = "FSC-A")
      updateSelectInput(session, "marker_y3", choices = currentGroup$markerFCS, selected = "FSC-H")
    } else {
      updateSelectInput(session, "marker_x3", choices = currentGroup$markerFCS)
      updateSelectInput(session, "marker_y3", choices = currentGroup$markerFCS)
    }
  })
  
  observeEvent(input$HELP_FCS_UPLOAD,{
    
    showModal(modalDialog(
      title = "About FCS uploading",
      "Indicate the number of file groups you wish to standardize. One file group = one manipulation. If you already have a file containing only the beads, you can import it. You still need to perform a QC on them to clean them of doublets.  ",
      easyClose = TRUE,
      footer = NULL
    ))
    
  })
  observeEvent(input$HELP_SEPARATE_BEADS,{
    
    showModal(modalDialog(
      title = "About SEPARATE BEADS FROM CELLS ",
      "The aim here is to select the cluster(s) corresponding to the beads and then move on to the 'Beads QC' section. ",
      easyClose = TRUE,
      footer = NULL
    ))
    
  })
  
  observeEvent(input$HELP_PREPROCESSING,{
    
    showModal(modalDialog(
      title = "ABOUT PREPROCESSING",
      "It is not mandatory to transform before extracting the beads and normalizing. If you do not wish to transform, please select 'None' and press 'Transform'.",
      easyClose = TRUE,
      footer = NULL
    ))
    
  })
  
  observeEvent(input$HELP_KMEANS,{
    
    showModal(modalDialog(
      title = "ABOUT KMEANS",
      " You can leave the number of clusters at 5 in most cases. Clustering will be applied to the concatenate of files in the group selected at the top of the page.'The Mini-batch K-means clustering algorithm is a version of the standard K-means algorithm in machine learning. It uses small, random, fixed-size batches of data to store in memory, and then with each iteration, a random sample of the data is collected and used to update the clusters.' ",
      easyClose = TRUE,
      footer = NULL
    ))
    
  })
  
  
  ################################################################################
  ##################### BEAD EXTRACTION ##########################################
  
  
  ####### STEP 1 : PREPROCESSING
  
  observeEvent(input$uploadTransformation,{
    tryCatch({
      req(input$groupTransfo)
      currentGroup<-listObjects[[input$groupTransfo]]
      
      currentGroup$csvFileTransfo<-read.table(input$uploadTransformation$datapath, row.names = 1, header=TRUE, sep=",", quote="")
      currentGroup$data.table<-currentGroup$csvFileTransfo
      listObjects[[input$groupTransfo]]<-currentGroup
      
    }, error = function(e) {
      
      showModal(modalDialog(
        title = "Error",
        paste("An error occurred:", e$message),
        easyClose = TRUE
      ))
    })
    
    
  })
  
  
  # Export transformed FCS
  output$exportFCSTransfo <- downloadHandler(
    filename = function() {
      paste0(input$groupTransfo, "_", Sys.Date(), "_transf", ".zip")
    },
    content = function(file) {
      temp_dir <- tempdir()
      file_paths <- list()
      
      for (id in 1:length(listObjects[[input$groupTransfo]]$flow.frames.deconcatenate)) {
        name <- basename(listObjects[[input$groupTransfo]]$original_filenames[id])
        new_name <- paste0("transf_", name)
        file_path <- file.path(temp_dir, new_name)
        
        write.FCS(listObjects[[input$groupTransfo]]$flow.frames.deconcatenate[[id]], file_path)
        file_paths <- c(file_paths, file_path)
      }
      
      zip_path <- file.path(temp_dir, "custom.zip")
      zip(zip_path, files = file_paths, flags = "-j")
      
      file.copy(zip_path, file)
    }
  )
  

  groupData <- reactive({
    req(input$groupTransfo)
    currentGroup <- listObjects[[input$groupTransfo]]
    req(currentGroup)
    
    if (!is.null(currentGroup$flow.frames.transformed)) {
      return(currentGroup$flow.frames.transformed)  
    } else {
      return(currentGroup$flow.frames)  
    }
  })
  
  
  output$plot_barplot <- renderPlot({
    tryCatch({
      req(input$marker_x)
      data <- groupData()  
      req(data)
      
      
      channel_data <- data@exprs[, input$marker_x]
      
      
      hist(
        channel_data,
        col = "#BCEE68",
        border = "#BCEE68",
        xlab = paste("Fluo Expression :",input$marker_x),
        ylab = "nCells",
        main = input$groupTransfo,
        breaks = input$num_breaks
      )
    }, error = function(e) {
      plot.new()
      text(0.5, 0.5, "Plot no available", cex = 1.2)
    })
  })
  
  
  output$plot_density <- renderPlot({
    tryCatch({
      data <- groupData()  
      req(data)
      
      if (nrow(data@exprs) > 10000) {
        subset_indices <- sample(nrow(data@exprs), 10000)
        exprs(data) <- data@exprs[subset_indices, ]
      }
      
      plotDens(
        data,
        channels = c(input$marker_x, input$marker_y),
        main = input$groupTransfo  
      )
    },  error = function(e) {
      plot.new()
      text(0.5, 0.5, "Plot no available", cex = 1.2)
      
    })
  })
  
  
  # Transformation 
  observeEvent(input$applyTransformation, {
    tryCatch({
      req(input$groupTransfo)
      currentGroup <- listObjects[[input$groupTransfo]]
      req(currentGroup)
      
      if (input$compensation) {
        currentGroup$flow.frames <- FlowCIPHE::compensate.CIPHE(currentGroup$flow.frames)
      }
      
      if (input$transfo == "arcsinh") {
        arg <- c(currentGroup$data.table[, 1])
        numeric_rows <- grepl("^[0-9]+(\\.[0-9]+)?$", arg)
        filtered_table <- currentGroup$data.table[numeric_rows, , drop = FALSE]
        marker <- c(rownames(filtered_table))
        arg <- as.numeric(filtered_table[, 1])
        currentGroup$flow.frames.transformed <- arcsinhCIPHE(currentGroup$flow.frames, marker = marker, arg)
      } else if (input$transfo == "None") {
        currentGroup$flow.frames.transformed <- currentGroup$flow.frames
      }
      
      
      if ("FlagNormApp" %in% currentGroup$markerFCS){
        currentGroup$flow.frames.deconcatenate<-FlowCIPHE::deconcatenate.FCS.CIPHE(currentGroup$flow.frames.transformed, "FlagNormApp")
      }else{
        currentGroup$flow.frames.deconcatenate<-list(currentGroup$flow.frames.transformed)
      }
      
      
      
      listObjects[[input$groupTransfo]] <- currentGroup
      
    }, error = function(e) {
      showModal(modalDialog(
        title = "Error",
        paste("An error occurred:", e$message),
        easyClose = TRUE
      ))
    })
  })
  
  
  ####### STEP 2 : PCA
  
  # observeEvent(input$HELP_PCA,{
  #   
  #   showModal(modalDialog(
  #     title = "ABOUT PCA",
  #     "PCA will be applied to all selected markers. PCA will be applied to the concatenate of files in the group selected at the top of the page.",
  #     easyClose = TRUE,
  #     footer = NULL
  #   ))
  #   
  # })
  # groupPCA <- reactive({
  #   req(input$groupTransfo)
  #   currentGroup <- listObjects[[input$groupTransfo]]
  #   req(currentGroup)
  #   req(currentGroup$flow.frames.PCA)
  #   
  #   return(currentGroup$flow.frames.PCA)
  #   
  #   
  # })
  # 
  # output$plot_PCA <- renderPlot({
  #   data <- groupPCA()
  #   req(data)
  #   
  #   if (nrow(data@exprs) > 10000) {
  #     subset_indices <- sample(nrow(data@exprs), 10000)
  #     exprs(data) <- data@exprs[subset_indices, ]
  #   }
  #   
  #   plotDens(
  #     data,
  #     channels = c("PC1", "PC2"),
  #     main = paste("PCA for Group:", input$groupTransfo)
  #   )
  # })
  # 
  # observeEvent(input$RunPCA, {
  #   tryCatch({
  #     req(input$groupTransfo)
  #     currentGroup <- listObjects[[input$groupTransfo]]
  #     req(currentGroup)
  #     req(input$marker_PCA)
  #     
  #     if (length(input$marker_PCA) < 2) {
  #       showModal(modalDialog(
  #         title = "Error",
  #         "Please select at least two markers for PCA.",
  #         easyClose = TRUE
  #       ))
  #       return()
  #     }
  #     
  #     transformed_data <- currentGroup$flow.frames.transformed
  #     req(transformed_data)
  #     
  #     data <- transformed_data@exprs[, input$marker_PCA, drop = FALSE]
  #     if (ncol(data) < 2) stop("Not enough markers for PCA.")
  #     
  #     prc <- princomp(data)
  #     
  #     pc1 <- prc$scores[, 1]
  #     pc2 <- prc$scores[, 2]
  #     
  #     transformed_data <- FlowCIPHE::enrich.FCS.CIPHE(transformed_data, as.vector(pc1), "PC1")
  #     transformed_data <- FlowCIPHE::enrich.FCS.CIPHE(transformed_data, as.vector(pc2), "PC2")
  #     
  #     currentGroup$flow.frames.PCA <- transformed_data
  #     listObjects[[input$groupTransfo]] <- currentGroup
  #   }, error = function(e) {
  #     showModal(modalDialog(
  #       title = "Error",
  #       paste("An error occurred:", e$message),
  #       easyClose = TRUE
  #     ))
  #   })
  #   
  # })
  
  ####### STEP 2 : KMEANS
  
  groupClusters <- reactive({
    
    req(input$groupTransfo)
    currentGroup <- listObjects[[input$groupTransfo]]
    req(currentGroup)
    req(currentGroup$flow.frames.mnKmeans)
    
    return(currentGroup$flow.frames.mnKmeans)
  })
  
  # run mnkmeans
  observeEvent(input$RunmnKmeans, {
    tryCatch({
      progress <- Progress$new()
      progress$set(message = "Run mnKmneans", value = 0.5)
      req(input$groupTransfo)
      currentGroup <- listObjects[[input$groupTransfo]]
      req(currentGroup)
      req(currentGroup$flow.frames.transformed)
      
      mini_batch_result <- ClusterR::MiniBatchKmeans(
        currentGroup$flow.frames.transformed@exprs[, c("FSC-A", "SSC-A")],
        clusters = input$nbCluster,
        batch_size = 1000,
        max_iters = 100,
        initializer = "kmeans++"
      )
      
      mnKmeans <- ClusterR::predict_MBatchKMeans(
        data = currentGroup$flow.frames.transformed@exprs[, c("FSC-A", "SSC-A")],
        CENTROIDS = mini_batch_result$centroids
      )
      
      mnKmeans_result <- as.matrix(mnKmeans)
      colnames(mnKmeans_result) <- "mnKmeans"
      currentGroup$flow.frames.mnKmeans <- FlowCIPHE::enrich.FCS.CIPHE(currentGroup$flow.frames.transformed, mnKmeans_result, "mnKmeans")
      listObjects[[input$groupTransfo]] <- currentGroup
      
      updateSelectizeInput(session, "clusters",
                           choices = unique(currentGroup$flow.frames.mnKmeans@exprs[, "mnKmeans"]),
                           selected = unique(currentGroup$flow.frames.mnKmeans@exprs[, "mnKmeans"]))
      progress$close()
    }, error = function(e) {
      showModal(modalDialog(
        title = "Error",
        paste("An error occurred:", e$message),
        easyClose = TRUE
      ))
    })
  })
  
  ####### STEP 3 : SEPARATE BEADS FROM CELLS
  
  output$plot_histo_clusters <- renderPlot({
    tryCatch({
      data <- groupClusters()
      req(data)
      
      dataFCS <- Subset(data, data@exprs[, "mnKmeans"] == input$clusters)
      hist(
        dataFCS@exprs[, input$marker_x2],
        breaks = 100,
        col = "#BCEE68",
        border = "#BCEE68",
        xlab = input$marker_x2,
        main = paste0("Cluster(s): ", paste(input$clusters, collapse = ", "))
      )
    }, error = function(e) {
      plot.new()
      text(0.5, 0.5, "No histogram available.", cex = 1.2)
    })
  })
  
  # plot clusters
  output$plot_clusters <- renderPlot({
    tryCatch({
      data <- groupClusters()
      req(data)
      
      dataFCS <- Subset(data, data@exprs[, "mnKmeans"] == input$clusters)
      if (nrow(dataFCS@exprs) > 10000) {
        subset_indices <- sample(nrow(dataFCS@exprs), 10000)
        exprs(dataFCS) <- dataFCS@exprs[subset_indices, ]
      }
      
      if (input$colorPlot == "Density") {
        plotDens(
          dataFCS,
          channels = c(input$marker_x2, input$marker_y2),
          main = paste0("Cluster(s): ", paste(input$clusters, collapse = ", "))
        )
      } else {
        data <- as.data.frame(dataFCS@exprs)
        predefined_colors <- setNames(
          colorRampPalette(brewer.pal(9, "Set1"))(20),
          as.character(1:20)
        )
        
        unique_clusters <- unique(data$mnKmeans)
        cluster_colors <- predefined_colors[as.character(unique_clusters)]
        
        plot(
          data[[input$marker_x2]],
          data[[input$marker_y2]],
          col = cluster_colors[as.character(data$mnKmeans)],
          pch = 16,
          xlab = input$marker_x2,
          ylab = input$marker_y2,
          main = paste0("Cluster(s): ", paste(unique_clusters, collapse = ", "))
        )
        
        legend(
          "topright",
          legend = unique_clusters,
          col = cluster_colors,
          pch = 16,
          title = "Clusters"
        )
      }
    }, error = function(e) {
      plot.new()
      text(0.5, 0.5, "No clusters available.", cex = 1.2)
    })
  })
  
  
  observeEvent(input$table_group_transfo_cell_edit, {
    
    req(input$groupTransfo)
    currentGroup <- listObjects[[input$groupTransfo]]
    row<-c(input$table_group_transfo_cell_edit$row)
    
    value<-c(input$table_group_transfo_cell_edit$value)
    currentGroup$data.table[row,1]<-value
    listObjects[[input$groupTransfo]] <- currentGroup
    
  })
  
  
  
  # Attribute cluster to beads
  observeEvent(input$beadAttribution, {
    tryCatch({
      req(input$groupTransfo)
      
      
      currentGroup <- listObjects[[input$groupTransfo]]
      
      
      currentGroup$beads <- NULL
      currentGroup$flow.frames.beads <- NULL
      currentGroup$flow.frames.cells <- NULL
      
      
      currentGroup$beads <- input$clusters
      currentGroup$flow.frames.beads <- Subset(
        currentGroup$flow.frames.mnKmeans,
        currentGroup$flow.frames.mnKmeans@exprs[,"mnKmeans"] %in% input$clusters
      )
      currentGroup$flow.frames.cells <- Subset(
        currentGroup$flow.frames.mnKmeans,
        !(currentGroup$flow.frames.mnKmeans@exprs[, "mnKmeans"] %in% input$clusters)
      )
      
      
      listObjects[[input$groupTransfo]] <- currentGroup
      
      showNotification(
        paste0("You have assigned the cluster(s):", paste(input$clusters, collapse = ", "),
               " to beads for group: ", input$groupTransfo),
        type = "message",
        duration = 3
      )
    }, error = function(e) {
      showModal(modalDialog(
        title = "Error",
        paste("An error occurred:", e$message),
        easyClose = TRUE
      ))
    })
  })
  
  
  
  ########################################################################
  ##################### BEAD QC ##########################################
  
  output$plot_QC <- renderPlot({
    tryCatch({
      beads <- groupBeads()
      req(beads)
      bead_data <- exprs(beads)
      
      
      if (nrow(bead_data) > 10000) {
        sampled_indices <- sample(nrow(bead_data), 10000)
        bead_data <- bead_data[sampled_indices, ]
      }
      
      beads_sampled <- flowFrame(bead_data)
      
      plotDens(
        beads_sampled,
        channels = c(input$marker_x3, input$marker_y3),
        main = paste("Beads for Group:", input$groupQC),
        xlim=input$X_scale,
        ylim=input$Y_scale
      )
    }, error = function(e) {
      plot.new()
      text(0.5, 0.5, "No beads available ", cex = 1.2)
    })
  })
  
  # remove doublets
  observeEvent(input$removeDoublets, {
    tryCatch({
      req(input$groupQC)
      
      for (group_name in names(listObjects)) {
        currentGroup <- listObjects[[group_name]]
        req(currentGroup)
        req(currentGroup$flow.frames.beads)
        
        gate.singlet <- openCyto:::.singletGate(currentGroup$flow.frames.beads, channels = c("FSC-A", "FSC-H"))
        not.selected.singlet <- notSubFrame(
          obj = currentGroup$flow.frames.beads,
          channels = c("FSC-A", "FSC-H"),
          filter = gate.singlet@boundaries
        )
        exprs(currentGroup$flow.frames.beads) <- currentGroup$flow.frames.beads@exprs[-not.selected.singlet@index, ]
        
        gate <- flowDensity(
          currentGroup$flow.frames.beads,
          c("FSC-A", "FSC-H"),
          position = c(FALSE, FALSE),
          ellip.gate = TRUE,
          use.percentile = c(TRUE, TRUE),
          percentile = c(0.9, 0.9),
          scale = 0.99
        )
        currentGroup$flow.frames.beads <- getflowFrame(gate)
        
        
        listObjects[[group_name]] <- currentGroup
      }
      
    }, error = function(e) {
      showModal(modalDialog(
        title = "Error",
        paste("An error occurred:", e$message),
        easyClose = TRUE
      ))
    })
  })
  
  # Transform data in logicle to better visualization
  observeEvent(input$transformVizu, {
    tryCatch({
      if (!is.null(listObjects)) {
        for (group_name in names(listObjects)) {
          currentGroup <- listObjects[[group_name]]
          req(currentGroup)
          req(currentGroup$flow.frames.beads)
          
          
          currentGroup$flow.frames.beads <- FlowCIPHE::logicle.CIPHE(currentGroup$flow.frames.beads, value = 500)
          
          
          if (!is.null(currentGroup$flow.frames.beads.gated.regate)) {
            currentGroup$flow.frames.beads.gated.regate <- FlowCIPHE::logicle.CIPHE(currentGroup$flow.frames.beads.gated.regate, value = 500)
          }
          
          
          if (!is.null(currentGroup$flow.frames.beads.gated)) {
            currentGroup$flow.frames.beads.gated <- FlowCIPHE::logicle.CIPHE(currentGroup$flow.frames.beads.gated, value = 500)
          }
          
          listObjects[[group_name]] <- currentGroup
        }
      }
    }, error = function(e) {
      showModal(modalDialog(
        title = "Error",
        paste("An error occurred:", e$message),
        easyClose = TRUE
      ))
    })
  })
  
  
  
  groupBeads <- reactive({
    req(input$groupQC)
    currentGroup <- listObjects[[input$groupQC]]
    req(currentGroup)
    req(currentGroup$flow.frames.beads)
    
    return(currentGroup$flow.frames.beads)  
    
    
  })
  
  groupBeadsGated <- reactive({
    req(input$groupQC)
    currentGroup <- listObjects[[input$groupQC]]
    req(currentGroup)
    req(currentGroup$flow.frames.beads.gated)
    return(currentGroup$flow.frames.beads.gated)
    
  })
  
  
  output$plot_histogram_peak2<- renderPlot({
    tryCatch({
      req(input$groupQC)
      currentGroup <- listObjects[[input$groupQC]]
      req(currentGroup)
      req(currentGroup$flow.frames.beads)
      beads<-currentGroup$flow.frames.beads
      hist(
        beads@exprs[, input$marker_x4],
        breaks = 100,
        col = "#BCEE68",
        border ="#BCEE68",
        xlab = input$marker_x4,
        main = paste("Beads before gating for :", input$groupQC)
      )
    }, error = function(e) {
      plot.new()
      text(0.5, 0.5, "No available", cex = 1.2)
    })
  })
  
  
  # identify modes
  observeEvent(input$mclustMode, {
    tryCatch({
    for (group_name in names(listObjects)) {
      currentGroup <- listObjects[[group_name]]
      req(currentGroup$flow.frames.beads)
      
      if (!(input$marker_x4 %in% colnames(currentGroup$flow.frames.beads))){
        marker<-currentGroup$mapping[[input$marker_x4]]
        gmm_result <- Mclust(as.matrix(currentGroup$flow.frames.beads@exprs[, marker]), G = input$nbPeaks)
      
      } else{
        
        gmm_result <- Mclust(as.matrix(currentGroup$flow.frames.beads@exprs[, input$marker_x4]), G = input$nbPeaks)
      }
      
      
 
      
      if ("mclust" %in% colnames(exprs(currentGroup$flow.frames.beads))) {
        exprs(currentGroup$flow.frames.beads) <- exprs(currentGroup$flow.frames.beads)[, colnames(exprs(currentGroup$flow.frames.beads)) != "mclust"]
      }
      
      
      currentGroup$flow.frames.beads <- FlowCIPHE::enrich.FCS.CIPHE(currentGroup$flow.frames.beads, gmm_result$classification, "mclust")
      listObjects[[group_name]] <- currentGroup
    }
    
    
    
    
    
    output$peakChoice <- renderUI({
      
      
      radioButtons(
        inputId = "chosenPeak",
        label = "Select Peak", 
        choices = sort(unique(gmm_result$classification))
      )
    })
  }, error = function(e) {
    showModal(modalDialog(
      title = "Error",
      paste("An error occurred:", e$message),
      easyClose = TRUE
    ))
  })
  })
  
  
  output$plot_histogram_peak3 <- renderPlot({
    
    tryCatch({
      
      req(input$groupQC)
      currentGroup <- listObjects[[input$groupQC]]
      req(currentGroup)
      req(currentGroup$flow.frames.beads)
      
      df <- as.data.frame(currentGroup$flow.frames.beads@exprs)
      req(df$mclust)
      ggplot(df, aes(df[[input$marker_x4]], fill = as.factor(df$mclust))) +
        geom_histogram() +
        labs(x = input$marker_x4, fill = "mclust")
    }, error = function(e) {
      plot.new()
      text(0.5, 0.5, "No available", cex = 1.2)
    })
  })
  
  
  # Choose peak
  observeEvent(input$chosenPeak, {
    req(input$groupQC)
    req(input$chosenPeak)
    
    
    currentGroup <- listObjects[[input$groupQC]]
    
    
    currentGroup$flow.frames.beads.gated <- currentGroup$flow.frames.beads[
      exprs(currentGroup$flow.frames.beads)[, "mclust"] == input$chosenPeak
    ]
    
    
    listObjects[[input$groupQC]] <- currentGroup 
    
    
    output$plot_peak <- renderPlot({
      tryCatch({
        req(input$groupQC)  
        currentGroup <- listObjects[[input$groupQC]]
        
        req(input$chosenPeak) 
        
        
        plotDens(
          currentGroup$flow.frames.beads.gated,
          channels = c(input$marker_x5, input$marker_y4),
          main = paste("Beads peaks n°", input$chosenPeak, "for group :", input$groupQC),
          cex = 2.5
        )
        
        
        if (!is.null(currentGroup$flow.frames.beads.gated.regate)) {
          gated_data <- currentGroup$flow.frames.beads.gated.regate@exprs
          points(
            gated_data[, input$marker_x5],
            gated_data[, input$marker_y4],
            col = "green",
            pch = 20
          )
        }
      }, error = function(e) {
        
        plot.new()
        text(0.5, 0.5, "No available", cex = 1.2)
      })
    })
  })
  
  
  observeEvent(input$reset_gate, {
    req(input$groupQC)
    tryCatch({
    currentGroup <- listObjects[[input$groupQC]]
    
    
    currentGroup$gated_population <- NULL
    currentGroup$flow.frames.beads.gated.regate<-NULL
    listObjects[[input$groupQC]] <- currentGroup
    }, error = function(e) {
      plot.new()
      text(0.5, 0.5, "No available", cex = 1.2)
    })
  })
  
  
  output$plot_histogram_peak<- renderPlot({
    tryCatch({
      req(input$groupQC)
      currentGroup<-listObjects[[input$groupQC]]
      
      req(input$chosenPeak)
      req(req(input$chosenPeak))
      hist(
        currentGroup$flow.frames.beads.gated@exprs[, input$marker_x5],
        breaks = 100,
        col = "#BCEE68",
        border ="#BCEE68",
        xlab = input$marker_x5,
        main =  paste("Beads peaks n°",input$chosenPeak, "for group :" , input$groupQC)
      )
      if (!is.null(currentGroup$flow.frames.beads.gated.regate)) {
        hist(
          currentGroup$flow.frames.beads.gated.regate@exprs[, input$marker_x5],
          breaks = 100,
          col = "#BCEE68",
          border ="#BCEE68",
          xlab = input$marker_x5,
          main =  paste("RE-GATE Beads peaks n°",input$chosenPeak, "for group :" , input$groupQC)
        )
      }
    }, error = function(e) {
      plot.new()
      text(0.5, 0.5, "No available", cex = 1.2)
    })
  })
  
  
  # Gating beads
  observeEvent(input$plot_brush, {
    
      req(input$groupQC)  
    
      currentGroup <- listObjects[[input$groupQC]]
      
      
      brush <- input$plot_brush
      req(brush)
      tryCatch({
      
      x_min <- brush$xmin
      x_max <- brush$xmax
      y_min <- brush$ymin
      y_max <- brush$ymax
      
      
      channel_x <- currentGroup$flow.frames.beads.gated@exprs[, input$marker_x5]
      channel_y <- currentGroup$flow.frames.beads.gated@exprs[, input$marker_y4]
      
      selected_indices <- which(
        channel_x >= x_min & channel_x <= x_max &
          channel_y >= y_min & channel_y <= y_max
      )
      
      
      gated_data <- currentGroup$flow.frames.beads.gated[selected_indices, ]
      
      
      currentGroup$flow.frames.beads.gated.regate <- gated_data
      listObjects[[input$groupQC]] <- currentGroup

    }, error = function(e) {
      showModal(modalDialog(
        title = "Error",
        paste("An error occurred:", e$message),
        easyClose = TRUE
      ))
    })
  })
  
  # untransform 
  observeEvent(input$UnTransformVizu, {
    req(input$groupQC)
    tryCatch({

      for (group_name in names(listObjects)) {
        
        currentGroup <- listObjects[[group_name]]
        
        req(currentGroup)
        
        req(currentGroup$flow.frames.beads)
        
        currentGroup$flow.frames.beads <- invers.logicle.CIPHE(currentGroup$flow.frames.beads, 500)
        
        if (!is.null(currentGroup$flow.frames.beads.gated)) {
          currentGroup$flow.frames.beads.gated <- invers.logicle.CIPHE(currentGroup$flow.frames.beads.gated, 500)
        }
        if (!is.null(currentGroup$flow.frames.beads.gated.regate)) {
          currentGroup$flow.frames.beads.gated.regate <- invers.logicle.CIPHE(currentGroup$flow.frames.beads.gated.regate, 500)
        }
        
        listObjects[[group_name]] <- currentGroup
        
      }
    }, error = function(e) {
      plot.new()
      text(0.5, 0.5, "No  available ", cex = 1.2)
    })
  })
  
  # valide bead for normalization
  observeEvent(input$valideBeads,{
    req(input$groupQC)
    for (group_name in names(listObjects)) {
      currentGroup <- listObjects[[group_name]]
      if (!is.null(currentGroup$flow.frames.beads.gated.regate)) {
        currentGroup$flow.frames.beads.gated <- currentGroup$flow.frames.beads.gated.regate
      }
      
      listObjects[[group_name]] <- currentGroup
    }
    
  })
  
  
  
  
  ################################################################################
  ##################### NORMALIZATION  ##########################################
  
  
  
  # Normalization Run
  observeEvent(input$Normalize, {
    tryCatch({
      
      referenceGroup <- listObjects[["group_1"]]
      
      withProgress(message = "Normalization in progress...", value = 0, {
        
        
        # Normalization in each group
        for (group_name in names(listObjects)) {
          if (group_name != "group_1") {
            
            names<-c()
            
            currentGroup <- listObjects[[group_name]]
     
            channels<-currentGroup$commonMarkers
            
            
            # If there are unmatch columns between ref and group
            if (!is.null(currentGroup$mapping)){
              
              refMarker<-names(unlist(currentGroup$mapping))
              
              currentMarker<-unlist(currentGroup$mapping)
              
              for(i in length(refMarker)){
                
                if(!is.null(referenceGroup$flow.frames)){
                  
                  colnames(referenceGroup$flow.frames)<-gsub(refMarker[i], currentMarker[i], colnames(referenceGroup$flow.frames))
                }
                if(!is.null(referenceGroup$flow.frames.cells)){
                  
                  colnames(referenceGroup$flow.frames.cells)<-gsub(refMarker[i], currentMarker[i], colnames(referenceGroup$flow.frames.cells))
                }
                
              }
              
              channels<-c(channels,unlist(currentGroup$mapping))
              
              }
              
           
            
            if(!is.null(referenceGroup$flow.frames)){
              
              fcs<-referenceGroup$flow.frames
              
              if (nrow(referenceGroup$flow.frames@exprs) > 10000) {
                
              subset_indices <- sample(nrow(referenceGroup$flow.frames@exprs), 10000)
              
              exprs(fcs) <- fcs@exprs[subset_indices, ]
              
              }
              concat<-as.data.frame(fcs@exprs[,channels])
              
              concat$group<-rep(paste0("group_1_ref"),nrow(concat))
              
              if(!is.null(referenceGroup$flow.frames.cells)){
                
                fcs<-referenceGroup$flow.frames.cells
              
                if (nrow(referenceGroup$flow.frames.cells@exprs) > 10000) {
                  
                  subset_indices <- sample(nrow(referenceGroup$flow.frames.cells@exprs), 10000)
                  
                exprs(fcs) <- fcs@exprs[subset_indices, ]
                
                }
                
                concatOnlyCells<-as.data.frame(fcs@exprs[,channels])
                
                concatOnlyCells$group<-rep(paste0("group_1_ref_onlyCells"),nrow(concatOnlyCells))
                
              }
              
            }else{
              
              concat<-data.frame()
            }
              fcs<-currentGroup$flow.frames
              
              if (nrow(currentGroup$flow.frames@exprs) > 10000) {
                
                subset_indices <- sample(nrow(currentGroup$flow.frames@exprs), 10000)
           
                exprs(fcs) <- fcs@exprs[subset_indices, ]
              }
            dfCurrentGroupBefNorm<-as.data.frame(fcs@exprs[,channels])
            
            dfCurrentGroupBefNorm$group<-rep(paste0(group_name,"_before_norm"), nrow(dfCurrentGroupBefNorm))
            
            concat<-rbind(concat, dfCurrentGroupBefNorm)
            
            if(!is.null(currentGroup$flow.frames.cells)){
              fcs<-currentGroup$flow.frames.cells
              
              if (nrow(currentGroup$flow.frames.cells@exprs) > 10000) {
              subset_indices <- sample(nrow(currentGroup$flow.frames.cells@exprs), 10000)
            
              exprs(fcs) <- fcs@exprs[subset_indices, ]
              }
              dfCurrentGroupBefNormOnlyCells<-as.data.frame(fcs@exprs[,channels])
              
              dfCurrentGroupBefNormOnlyCells$group<-rep(paste0(group_name,"_before_norm_onlyCells"), nrow(dfCurrentGroupBefNormOnlyCells))
              
              concatOnlyCells<-rbind(concatOnlyCells, dfCurrentGroupBefNormOnlyCells)
              
            }
            
            harmonization_factors<-c()
            
            if(!is.null(currentGroup$flow.frames.mnKmeans)){
              currentGroup$flow.frames.normalized<-currentGroup$flow.frames.mnKmeans
            }else{
              currentGroup$flow.frames.normalized<-currentGroup$flow.frames
            }
            if(!is.null(currentGroup$flow.frames.cells)){
              currentGroup$flow.frames.cells.normalized<-currentGroup$flow.frames.cells
            }
            
            
            for (x in channels) {
              print(paste("Processing channel:", x, "for group:", group_name))
              
              if (!(x %in% colnames(referenceGroup$flow.frames.beads.gated))){
                print(names(currentGroup$mapping)[which(currentGroup$mapping == x)])
              ref_median <- median(exprs(referenceGroup$flow.frames.beads.gated)[,names(currentGroup$mapping)[which(currentGroup$mapping == x)]])

               
              } else{
              ref_median <- median(exprs(referenceGroup$flow.frames.beads.gated)[, x])
              
              }
              group_median <- median(exprs(currentGroup$flow.frames.beads.gated)[, x])
              harmonization_factor <- ref_median / group_median
              harmonization_factors<-c(harmonization_factors,round(harmonization_factor,4))
              
              print(paste("Harmonization factor for channel", x, ":", harmonization_factor))
              
              if(!is.null(currentGroup$flow.frames.mnKmeans)){
                exprs(currentGroup$flow.frames.normalized)[, x] <-
                  exprs(currentGroup$flow.frames.mnKmeans)[, x] * harmonization_factor
              }else{
                exprs(currentGroup$flow.frames.normalized)[, x] <-
                  exprs(currentGroup$flow.frames)[, x] * harmonization_factor
              }
              if(!is.null(currentGroup$flow.frames.cells)){
                exprs(currentGroup$flow.frames.cells.normalized)[, x] <-
                  exprs(currentGroup$flow.frames.cells)[, x] * harmonization_factor
              }
              
            }
            
 
            harmonization_table<-data.frame(channel=channels, Factor=harmonization_factors)
            rownames(harmonization_table) <- harmonization_table$channel
            harmonization_table$channel <- NULL
            
            if ("FlagNormApp" %in% currentGroup$markerFCS){
              
              currentGroup$flow.frames.normalized.deconcatenate<-FlowCIPHE::deconcatenate.FCS.CIPHE(currentGroup$flow.frames.normalized, "FlagNormApp")
              
            }else{
              
              currentGroup$flow.frames.normalized.deconcatenate<-list(currentGroup$flow.frames.normalized)
            }
            
            currentGroup$harmonization_table <- harmonization_table
            
            names<-c(names,group_name)
            fcs<-currentGroup$flow.frames.normalized
            if (nrow(currentGroup$flow.frames.normalized@exprs) > 10000) {
            subset_indices <- sample(nrow(currentGroup$flow.frames.normalized@exprs), 10000)
           
            exprs(fcs) <- fcs@exprs[subset_indices, ]
            }
            dfCurrentGroupAfterNorm<-as.data.frame(fcs@exprs[,channels])
            dfCurrentGroupAfterNorm$group<-rep(paste0(group_name,"_after_norm"), nrow(dfCurrentGroupAfterNorm))
            concat<-rbind(concat, dfCurrentGroupAfterNorm)
            
            currentGroup$concat<-concat
            
            
            if(!is.null(currentGroup$flow.frames.cells)){
              fcs<-currentGroup$flow.frames.cells.normalized
              if (nrow(currentGroup$flow.frames.cells.normalized@exprs) > 10000) {
              subset_indices <- sample(nrow(currentGroup$flow.frames.cells.normalized@exprs), 10000)

              exprs(fcs) <- fcs@exprs[subset_indices, ]
              }
              dfCurrentGroupAfterNormOnlyCells<-as.data.frame(fcs@exprs[,channels])
              
              dfCurrentGroupAfterNormOnlyCells$group<-rep(paste0(group_name,"_after_norm_onlyCells"), nrow(dfCurrentGroupAfterNormOnlyCells))
              concatOnlyCells<-rbind(concatOnlyCells, dfCurrentGroupAfterNormOnlyCells)
              currentGroup$concatOnlyCells<-concatOnlyCells 
            }
            
            currentGroup$channels<-channels
            listObjects[[group_name]]<-currentGroup
          }
        }
        
        })
  
        output$groupVizTable<-renderUI({
          groups<-names(listObjects)[! names(listObjects) %in% as.vector("group_1")]
          selectInput("groupVizTable", "Harmonization table for :", choices = groups, multiple=FALSE, width = 500)
        })

    
    }, error = function(e) {
      showModal(modalDialog(
        title = "Error",
        paste("An error occurred:", e$message),
        easyClose = TRUE
      ))
    })
  })
  
  
  # Display harmonization table
  output$harmTableVizu<-renderUI({
    
    req(input$groupVizTable)
    currentGroup<- listObjects[[input$groupVizTable]]
    req(currentGroup$harmonization_table)
    
    DT::datatable(
      currentGroup$harmonization_table
    ) %>%
      formatStyle(columns = "Factor", 
                  background = styleInterval(c(0, 0.8, 2)-1e-6, c("white", "lightblue", "white", "#F08080")))
    
    
  })
  
 
  observe({
    req(input$groupVizTable)
    req(listObjects[[input$groupVizTable]]$channels)
    updateSelectInput(session, "marker_x6",  choices = listObjects[[input$groupVizTable]]$channels)
    
    
  })
  
  observe({
    req(input$groupVizTable)
    currentGroup <- listObjects[[input$groupVizTable]]
    req(currentGroup$concat)
    req(input$cellsOrcellsAndBeads)
    if(input$cellsOrcellsAndBeads == "rawFCS"){
      updateSelectInput(session, "groupVerif", selected=unique(currentGroup$concat$group), choices = unique(currentGroup$concat$group))
    }else{
      updateSelectInput(session, "groupVerif", selected=unique(currentGroup$concatOnlyCells$group), choices = unique(currentGroup$concatOnlyCells$group))
    }
    
  })
  
  
  observe({
    req(input$groupVizTable)
    currentGroup <- listObjects[[input$groupVizTable]]
    req(currentGroup$concat)
    updateSelectInput(session, "cellsOrcellsAndBeads", selected="rawFCS" , choices = c("rawFCS"))
    
    req(currentGroup$concatOnlyCells)
    updateSelectInput(session, "cellsOrcellsAndBeads", selected="onlyCells" , choices = c("rawFCS", "onlyCells"))
  })
  
  
  # Plot after normalization
  output$vizAfterNorm<-renderPlot({
    
    req(input$groupVizTable)
    currentGroup <- listObjects[[input$groupVizTable]]
    req(input$cellsOrcellsAndBeads)
    req(currentGroup$concat)
    if(input$cellsOrcellsAndBeads == "rawFCS"){
      subset<-currentGroup$concat[which(currentGroup$concat$group %in% input$groupVerif),]
    }else{
      subset<-currentGroup$concatOnlyCells[which(currentGroup$concatOnlyCells$group %in% input$groupVerif),]
    }
    ggplot(subset, aes(x=subset[[input$marker_x6]], fill=group)) +
      geom_density(alpha=.25)+
      scale_x_log10() +
      xlab(input$marker_x6)
    
    
  })
  
  
  # Export harmonization table
  
  output$downloadHarmTable<- downloadHandler(

    filename = function() {
      paste("harmonization_table-", input$groupVizTable,"-",Sys.Date(), ".csv", sep="")
    },
    content = function(file) {
      write.csv(listObjects[[input$groupVizTable]]$harmonization_table, file, quote=FALSE)
    }
  )
  
  # Export transformation 
  
  output$exportTransformation <- downloadHandler(
    
    filename = function() {
      paste("transf_table-", input$groupTransfo,"-",Sys.Date(), ".csv", sep="")
    },
    content = function(file) {
      write.csv(listObjects[[input$groupTransfo]]$data.table, file, quote=FALSE)
    }
  )
  
  # Export normalized FCS
  
  output$downloadFCS <- downloadHandler(
    
    filename = function() {
      
      paste0(input$zipFileName, ".zip")
    },
    content = function(file) {
      
      temp_dir <- file.path(tempdir(), paste0("export_", Sys.time()))
      
      dir.create(temp_dir)
      
      file_paths <- list()
      
      for (group in names(listObjects)) {
        
        group_data <- listObjects[[group]]
        
        if (!is.null(group_data$flow.frames.normalized.deconcatenate)) {
          
          for (id in 1:length(group_data$flow.frames.normalized.deconcatenate)) {
            
            name <- basename(group_data$original_filenames[id])
            
            new_name <- paste0("Norm_", name)
            
            file_path <- file.path(temp_dir, new_name)
            
            
            flow_frame <- group_data$flow.frames.normalized.deconcatenate[[id]]
            
            
            if ("mnkmeans" %in% colnames(exprs(flow_frame))) {
              
              exprs(flow_frame) <- exprs(flow_frame)[, !(colnames(exprs(flow_frame)) == "mnkmeans")]
            }
            
            
            write.FCS(flow_frame, file_path)
            file_paths <- c(file_paths, file_path)
          }
        }
      }
      
      zip_path <- file.path(temp_dir, "custom.zip")
      zip(zip_path, files = file_paths, flags = "-j")
      
      file.copy(zip_path, file)
      
      unlink(temp_dir, recursive = TRUE)
    }
  )
  
}




























