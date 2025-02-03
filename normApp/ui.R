

library(shiny)
library(shinydashboard)
library(DT)
library(shinybusy)


ui <- dashboardPage(
  
  skin = "blue",
  dashboardHeader(title = 'Normalization App'),
  dashboardSidebar(
    sidebarMenu(
      id = "tab",
      menuItem(" Upload FCS", tabName = "upload", icon = icon("upload")),
      menuItem(" Beads extraction", tabName = "beadExtraction", icon = icon("scissors")),
      menuItem("Beads QC", tabName = "Beads_QC"),      
      
      menuItem("Normalization", tabName = "Normalization", icon = icon("balance-scale"))
      
      
    )
  ),
  dashboardBody(
    
    
    
    fluidRow(
      uiOutput("progressBoxes")
    ),
    
    
    
    ################################################################################
    ##################### UPLOAD DATA ##############################################
    
    conditionalPanel(
      
      condition = "input.tab == 'upload'",
      fluidPage(
        tabItems(
          
          tabItem(
            
            tabName = "upload",
            
            # Choose files to annotate
            box(width=12, title=p("Choose one or multiple FCS files",actionButton("HELP_FCS_UPLOAD", "", icon("question-circle"), class = "btn-xs")),
                
                status = "primary",
                
                solidHeader = TRUE,
                
                collapsible = TRUE,
                
                box(width=8,fluidRow(
                  HTML("<h4> <span style='background-color: #fdfd96;'> STEP 1 :</span>  Upload Files </h4>" ),
                  column(4,
                         selectInput(
                           inputId = "selector",
                           label = "How many groups?",
                           choices = 2:10
                         ),
                         HTML("<h4> <span style='color: green;font-weight: bold;'> GROUP 1 ->  </span><span style='font-weight: bold;'>REFERENCE GROUP </span></h4>" )
                         )),
                
                fluidRow(
                  column(5,
                         uiOutput(outputId = "file_inputs")
                  ),
                  column(5,  uiOutput(outputId = "beads_inputs")))),
                fluidRow(
                  box(width=8,
                      HTML("<h4> <span style='background-color: #fdfd96;'> STEP 2 :</span>  Select markers to normalize </h4>" ),
                                column(4,selectInput(
                                  "groupNorm",
                                  "Group",
                                  choices = c(""), multiple=FALSE,selected = NULL
                                )),
                                column(5, selectInput(
                                  "commonMarkers",
                                  "Common markers with Ref ",
                                  choices = c(""), multiple=TRUE,selected = NULL
                                )),
                                # column(2, tags$br(),actionButton("chooseMarkerNorm", label = "OK")),
                    
                      # selectInput(
                      #   "unNormMarkers", 
                      #   HTML("<b>Select markers <span style='background-color: #fdfd96;'>NOT</span> to be normalized</b>"), 
                      #   choices = c(""), 
                      #   selected = NULL, 
                      #   width = 500, 
                      #   multiple = TRUE
                      # )
                      
                  ),
                ),
                      fluidRow(box(width=8,HTML("<h5> <span style='background-color: #fdfd96;'> (OPTIONNAL)</span> Match markers between REF and others groups  : </h5>" ),
                               tags$br(),
                     
                               tags$br(),
                               tags$br(),
                               column(4,uiOutput("dynamic_selectors")),
                               
                               column(4,verbatimTextOutput("result"))))
             
              
                
            ),
            tags$br(),
            # Checkbox to determine if files are already clustered
            
          )
        ),
      )
    ),
    
    ################################################################################
    ##################### BEAD EXTRACTION ##########################################
    
    conditionalPanel(
      condition = "input.tab == 'beadExtraction'",
      
      uiOutput("fcsTransformed"),
      
      fluidPage(
        
        box(width = 12,HTML("<h4> Supervised Beads Extraction </h4>" ),HTML("<h5> This section allows you to extract beads from the selected group. If you have already imported extracted beads for the selected group, please ignore this section for this group and proceed to 'Bead QC'. Otherwise, follow the 4 steps below.</h5> "),
            column(6, selectizeInput("groupTransfo", "Select group", choices = c(""), multiple = FALSE, width="60%")),
            
        ),
        
        ####### STEP 1 : PREPROCESSING
        
        fluidRow(
          
          column(
            width = 12,
            box(
              width = 12,
              title = p("STEP 1 : Preprocessing",actionButton("HELP_PREPROCESSING", "", icon("question-circle"), class = "btn-xs")),
              status = "primary",
              solidHeader = TRUE,
              collapsible = TRUE,
              tabItems(
                tabItem(
                  tabName = "beadExtraction",
                  
                  # Option pour appliquer la matrice de compensation
                  checkboxInput("compensation", "Apply spillover matrix", value =FALSE),
                  uiOutput("compensation_message"),
                  
                  # Sélection de la transformation
                  fluidRow(
                    
                    column(
                      6, 
                      fileInput(
                        inputId = "uploadTransformation",
                        label = p("(Option) Upload CSV transformation file",actionButton("HELP_CSV_UPLOAD", "", icon("question-circle"), class = "btn-xs")),
                        buttonLabel = "Upload CSV",
                        placeholder = "No file selected",
                        multiple = FALSE,
                        accept = c(".csv"),
                        width = "100%",
                        
                      ),
                      selectInput("transfo", "Transformation", choices = list("arcsinh","None"), selected = "None", width = 200)
                    ), 
                    tags$br(),
                    
                  ),
                  
                  fluidRow(
                    column(
                      4,  
                      DTOutput("table_group_transfo"),
                      
                      tags$div(id = "status_message"),
                      actionButton(inputId = "applyTransformation", label = "Transform", style = "background-color: #BCEE68; border-width: 2px; border-color:  #BCEE68; margin-top: 10px;"),
                      downloadButton("exportTransformation", "Export CSV"),
                      downloadButton("exportFCSTransfo", "Export FCS")
                      
                    ),
                    
                    column(
                      8,
                      fluidRow(
                        column(
                          3,
                          selectInput("marker_x", "X-axis Marker", choices = c(""), selected = NULL, width = "100%")
                        ),
                        column(
                          3,
                          selectInput("marker_y", "Y-axis Marker", choices = c(""), selected = NULL, width = "100%")
                        ),
                        column(4, sliderInput(
                          inputId = "num_breaks",
                          label = "Number of Breaks:",
                          min = 5,
                          max = 100,
                          value = 30,
                          step = 5))
                      ),
                      fluidRow(
                        column(
                          6,
                          plotOutput("plot_density", height = "400px", width = "400px")
                        ),
                        column(
                          6,
                          plotOutput("plot_barplot", height = "400px", width = "400px")
                        )
                      )
                    )
                  )
                  
                )
              )
            )
            ,
            
            # ####### STEP 2 : PCA
            # 
            # box(width = 12,
            #     title = p("STEP 2 : PCA",actionButton("HELP_PCA", "", icon("question-circle"), class = "btn-xs")),
            #     status = "primary",
            #     solidHeader = TRUE,
            #     collapsible = TRUE,
            #     fluidRow(
            #       column(5,
            #              selectizeInput("marker_PCA", "Markers for PCA", choices = c(""), multiple = TRUE, width=700)),
            #       tags$br(),
            #       column(2, 
            #              # Apply transformation
            #              actionButton(inputId = "RunPCA", label = "Run PCA", style = " background-color: #BCEE68;border-width: 2px; border-color:  #BCEE68;")),
            #       column(10,
            #              
            #              plotOutput("plot_PCA", height = "500px", width = "500px"))
            #     ) ),
            box(width = 12,
                title = p("STEP 2 : Mini-batch K-means",actionButton("HELP_KMEANS", "", icon("question-circle"), class = "btn-xs")),
                status = "primary",
                solidHeader = TRUE,
                collapsible = TRUE,
                fluidRow(
                  column(6,
                         numericInput("nbCluster", "nb Clusters", min=2, max=50, value=5)),
                  tags$br(),
                  column(2,actionButton(inputId = "RunmnKmeans", label = "Run mini batch Kmeans", style = " background-color: #BCEE68;border-width: 2px; border-color: #BCEE68;"))
                  # column(2,numericInput("h", "h", value = 1.5, width=700))
                ) 
            ),
            box(width = 12,
                title = p("STEP 3 : Separate beads from cells ",actionButton("HELP_SEPARATE_BEADS", "", icon("question-circle"), class = "btn-xs")),
                status = "primary",
                solidHeader = TRUE,
                collapsible = TRUE,
                fluidRow(
                  column(6,
                         selectizeInput("clusters", "Select cluster(s) to visualize", choices = c(""), multiple = TRUE, width="100%")),
                  tags$br(),
                  
                  
                  
                ),
                tags$br(),
                fluidRow(
                  # Select markers for visualization
                  column(2,selectInput("marker_x2", "X-axis Marker", choices = c(""), selected = NULL, width = 200)),
                  column(2, selectInput("marker_y2", "Y-axis Marker", choices = c(""), selected = NULL, width = 200)),
                  column(4, selectInput("colorPlot", "Color by :", choices = c("Density","mnKmeansClusters"), selected = "mnKmeansClusters", width = 200))
                  
                ),
                fluidRow(
                  column(6,
                         # Display the plot
                         plotOutput("plot_clusters", height = "500px", width = "500px")),
                  column(6,
                         plotOutput("plot_histo_clusters", height = "500px", width = "500px")
                  ),
                  column(5,actionButton(inputId = "beadAttribution", label = "Assign this cluster to beads", style = " background-color: #BCEE68;border-width: 2px; border-color: #BCEE68;"))
                )
                
            ),
            
            
          )
        )
      )
    ),
    ########################################################################
    ##################### BEAD QC ##########################################
    
    conditionalPanel(
      condition = "input.tab == 'Beads_QC'",
      fluidPage(
        box(width = 12,HTML("<h4> Beads QC </h4>" ),HTML("<h5> This section allows you to clean the beads for each group and select the High peak.</h5> "),
            
            
        ),
        
        tabItems(
          tabItem(
            tabName = "Beads_QC",
            box(
              width = 12,
              
              title = "Beads Control Quality", 
              status = "primary",
              solidHeader = TRUE,
              collapsible = TRUE,
              
              fluidRow(
                box(width=12,
                    HTML("<h4><strong>    STEP 1 : Remove doublets </strong></h4>"),
                    column(5,
                           selectizeInput("groupQC", "Visualization group :", choices = c(""), multiple = FALSE, width="60%"),
                           plotOutput("plot_QC", height = "400px", width = "400px")),
                    column(
                      2,
                      selectInput(
                        "marker_x3", 
                        "X-axis Marker", 
                        choices = c(""), 
                        selected = NULL, 
                        width = "100%"
                      ),
                      selectInput(
                        "marker_y3", 
                        "Y-axis Marker", 
                        choices = c(""), 
                        selected = NULL, 
                        width = "100%"
                      ),
                      tags$br(),
                      
               
                      sliderInput("X_scale", "X scale",
                                  min = 1, max = 1000000,
                                  value = c(0,500000), step=100),
                      sliderInput("Y_scale", "Y scale",
                                  min = 1, max = 1000000,
                                  value = c(0,500000), step=100),
                      
                      actionButton(
                        inputId = "removeDoublets", 
                        label = "Remove doublets (All groups)", 
                        style = "background-color: #BCEE68; border-width: 2px; border-color: #BCEE68;", 
                        width = "60%"
                      ),
                      tags$br(),
                      tags$br(),
                      
                      
                    ),
                )
              ),
              
              box(width=12,
                  HTML("<h4><strong>    STEP 2 : Select mode </strong></h4>"),
                  tags$br(),
                  column(3, 
                         actionButton("transformVizu", "Transform data (All groups)", style = "background-color: #BCEE68; border-width: 2px; border-color: #BCEE68;"),
                         tags$br(),
                         tags$br(),
                         selectInput(
                           "marker_x4", 
                           "X-axis Marker", 
                           choices = c(""), 
                           selected = NULL, 
                           width = "100%"
                         ),
                         numericInput("nbPeaks", "nb of Peaks", min=1, max=20, value=2),
                         actionButton("mclustMode", "Launch mclust (All groups)", style = "background-color: #BCEE68; border-width: 2px; border-color: #BCEE68;"),
                         tags$br(),
                         uiOutput("peakChoice")),
                  
                  
                  column(4,  plotOutput("plot_histogram_peak2", height = "400px", width = "400px")),
                  column(4,  plotOutput("plot_histogram_peak3", height = "400px", width = "400px")),
                  
              ),
              
              box(width=12,
                  HTML("<h4><strong> (optionnal)    STEP 3 : Re-gate  </strong></h4>"),
                  tags$br(), 
                  
                  fluidRow(
                    
                    column(
                      3,
                      
                      selectInput(
                        "marker_x5", 
                        "X-axis Marker", 
                        choices = c(""), 
                        selected = NULL, 
                        width = "100%"
                      ),
                      # numericInput(
                      #   "gate_x", 
                      #   "Set Gate Threshold on X-axis", 
                      #   value = 0, 
                      #   min = 0, 
                      #   max = 100000, 
                      #   step = 0.1
                      # ),
                      # 
                      actionButton("reset_gate", "Reset Gate")
                      # tags$div(
                      #   
                      #   tags$span("←  More Stringent ", style = "color: red;"),
                      # 
                      # )
                    ),
                    column(
                      3,
                      
                      selectInput(
                        "marker_y4", 
                        "Y-axis Marker", 
                        choices = c(""), 
                        selected = NULL, 
                        width = "100%"
                      ),
                      # numericInput(
                      #   "gate_y", 
                      #   "Set Gate Threshold on Y-axis", 
                      #   value = 0, 
                      #   min = 0, 
                      #   max = 100000, 
                      #   step = 0.1
                      # ),
                      
                      tags$br(),
                      # column(
                      #   6,
                      #   
                      #   
                      #   actionButton(
                      #     inputId = "gatingCST", 
                      #     label = "Gating CST", 
                      #     style = "background-color: #BCEE68; border-width: 2px; border-color: #BCEE68;", 
                      #     width = "100%"
                      #   )
                      # ),
                      
                      
                      
                    ),
                    # column(6,box(
                    #   width = "60%", 
                    #   title = " Gating beads details", 
                    #   status = "primary", 
                    #   solidHeader = TRUE, 
                    #   collapsible = TRUE, 
                    #   tags$div(
                    #     style = "text-align: justify;",
                    #     tags$p("We aim to clean the beads using gates and select only one type of bead (e.g., CST Hi). The selected beads will be used for normalization."),
                    #     
                    #     tags$p(tags$b("PS: For CST beads in flow cytometry, beads are well-defined with X=G-PE-A and Y=BUV650-A.")),
                    #     
                    #     
                    #     tags$ol(
                    #       tags$li(tags$b("If your data is not transformed: "), "The first step is to check the box 'Transform data for better visualization'. Then, you can search for the marker pair where your beads are well-defined."),
                    #       tags$li(tags$b("Set Gate Threshold on X-axis "), "The beads you want to retrieve are above this threshold for X."),
                    #       tags$li(tags$b("Set Gate Threshold on Y-axis "), "The beads you want to retrieve are above this threshold for Y."),
                    #       tags$li(tags$b("Ellipse Gate: "), "Allows you to define the gate around the beads. Decrease the scale for higher stringency, and increase it for lower stringency.")
                    #     )
                    #   )
                    #   
                    #   
                    # ))
                  ),
                  
                  
                  
                  
                  tags$br(),
                  fluidRow(
                    
                    column(
                      4,
                      plotOutput("plot_peak", 
                                 brush = brushOpts(
                                   id = "plot_brush",       # L'ID du brush pour capturer la sélection
                                   resetOnNew = TRUE        # Réinitialise le brush lorsqu'un nouveau graphique est généré
                                 ),height = "400px", width = "400px"),
                      # plotOutput("plot_peak_gated", height = "400px", width = "400px")
                    ),
                    column(
                      4,
                      plotOutput("plot_histogram_peak", height = "400px", width = "400px"),
                      # plotOutput("plot_histogram_peak_gated", height = "400px", width = "400px")
                    ),
                    column(4,actionButton("UnTransformVizu", "unTransform data (All groups)", style = "background-color: #BCEE68; border-width: 2px; border-color: #BCEE68;"),
                           tags$br(),tags$br(),actionButton("valideBeads", "Valide beads (All groups)", style = "background-color: #BCEE68; border-width: 2px; border-color: #BCEE68;")),
                    
                    
                    
                    
                    
                    
                  ),
                  
                  
                 
                  
                  
              )
              
            )
          )
        )
      )
    )
    ,
    
    ################################################################################
    ##################### NORMALIZATION  ##########################################
    
    
    
    conditionalPanel(
      
      condition = "input.tab == 'Normalization'",
      fluidPage(
        tabItems(
          
          tabItem(
            
            tabName = "Normalization",
            
            # Choose files to annotate
            box(width=12,title="Normalization",
                
                status = "primary",
                
                solidHeader = TRUE,
                
                collapsible = TRUE,
                # fluidRow(
                #   
                #   box(width=8,
                #       HTML("<h4>  <span style='background-color: #fdfd96;'>STEP 1 :</span>  Select reference group </h4>" ),
                #       selectInput("referenceGroup", "", choices = c(""), multiple=FALSE,selected = NULL, width = 500)  )
                #   
                # ),
                # fluidRow(
                #   box(width=8,
                #       fluidRow( HTML("<h4> <span style='background-color: #fdfd96;'> STEP 2 :</span>  Select markers to normalize </h4>" ),
                #                 column(4,selectInput(
                #                   "groupNorm",
                #                   "Group",
                #                   choices = c(""), multiple=FALSE,selected = NULL
                #                 )),
                #                 column(5, selectInput(
                #                   "commonMarkers",
                #                   "Common markers ",
                #                   choices = c(""), multiple=TRUE,selected = NULL
                #                 )),
                #                 # column(2, tags$br(),actionButton("chooseMarkerNorm", label = "OK")),
                #                 uiOutput("slider_ui")),
                #       # selectInput(
                #       #   "unNormMarkers", 
                #       #   HTML("<b>Select markers <span style='background-color: #fdfd96;'>NOT</span> to be normalized</b>"), 
                #       #   choices = c(""), 
                #       #   selected = NULL, 
                #       #   width = 500, 
                #       #   multiple = TRUE
                #       # )
                #       fluidRow(HTML("<h5> <span style='background-color: #fdfd96;'> (OPTIONNAL)</span> Match markers between REF and others groups  : </h5>" ),
                #                tags$br(),
                #                actionButton("find_more", "Find more common columns"),
                #                tags$br(),
                #                tags$br(),
                #                DTOutput("mapping_table"))
                #   ),
                # ),
                # 
                # 
                # 
                
                
                column(2, actionButton(inputId = "Normalize", label = "Normalize", style = " background-color:  #BCEE68;border-width: 2px; border-color:  #BCEE68;"),
                       tags$br()),
                

                fluidRow(
                  tags$br(),
                  column(3,
                         tags$br(),
                         uiOutput("groupVizTable"),
                         uiOutput("harmTableVizu")),
                  
                  
                  column(2,tags$br(),
                         selectInput("cellsOrcellsAndBeads", "Choose", choices = c(""), selected = NULL, width = 200),
                         selectInput("marker_x6", "X-axis Marker", choices = c(""), selected = NULL, width = 200),
                         selectInput("groupVerif", "Group", choices = c(""), selected = NULL,multiple=TRUE, width = 200)),
                  column(4,plotOutput("vizAfterNorm",height = "500px", width = "500px"))
                ),
                fluidRow(
                  column(4,
                         textInput("zipFileName", "Enter ZIP file name (without extension):", value = "NormalizedFCS")),
                  column(4, downloadButton("downloadFCS", "Download normalized FCS"),tags$br(),downloadButton("downloadHarmTable", "Download harmonized Table")))
            ),
            
            
          )
        ),
        
      ),
      
    )
  )
  
  
  
)