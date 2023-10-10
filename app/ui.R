#if (interactive()) {
##################################################
# Define UI for application
ui <- fluidPage(
  
  # Capture user access information
  tags$head(
    tags$title("SIGNAL - Selection by Iterative pathway Group and Network Analysis Looping"),
    tags$script(src="getIP.js"),
    # sources below are for d3 network graph layout and interaction
    tags$script(src="mbostock_d3.js"),  # original file: https://mbostock.github.io/d3/talk/20111116/d3/d3.js
    tags$script(src="mbostock_d3_layout.js"),  # original file: https://mbostock.github.io/d3/talk/20111116/d3/d3.layout.js
    tags$script(src="custom_network.js"),
    tags$script(src="custom_network2.js")
  ),
  
  # Show a header using the header style
  headerPanel(includeHTML("header.html")),
  
  # style
  theme = "./css/signal.css",
  
  # use shinyjs
  useShinyjs(), br(),
  
  # Sidebar with a slider input for number of bins
  sidebarLayout(
    
    sidebarPanel(
      # Global site tag (gtag.js) - Google Analytics 
      tags$head(
        HTML("<!-- Global site tag (gtag.js) - Google Analytics -->"),
        tags$script(HTML("(function(w,d,s,l,i){w[l]=w[l]||[];w[l].push({'gtm.start':
                            new Date().getTime(),event:'gtm.js'});var f=d.getElementsByTagName(s)[0],
                             j=d.createElement(s),dl=l!='dataLayer'?'&l='+l:'';j.async=true;j.src=
                              'https://www.googletagmanager.com/gtm.js?id='+i+dl;f.parentNode.insertBefore(j,f);
                             })(window,document,'script','dataLayer','GTM-56JHMG');"
        ))
        
      ),
      tags$body(HTML('<noscript><iframe src="https://www.googletagmanager.com/ns.html?id=GTM-56JHMG"
                            height="0" width="0" style="display:none;visibility:hidden"></iframe></noscript>'
      )
      ),
      # Parameters to be selected
      selectInput(inputId = "organism",
                  label = "Select your organism:",
                  choices = c("Human", "Mouse")
      ),
      selectInput(inputId = "pathway",
                  label = "Select a Database for Enrichment Analysis:",
                  choices = c("KEGG: Biological Processes", "KEGG: Disease Pathways", "KEGG: All Pathways")
      ),
      pickerInput(inputId = 'pathway_alt',
                  label = "Select a Database for Enrichment Analysis:",
                  choices = list(
                    "KEGG:" = c('KEGG: Biological Processes', 
                                'KEGG: Disease Pathways', 
                                'KEGG: All Pathways'),
                    "GO:" = c('GO: Biological Processes', 
                              'GO: Molecular Functions', 
                              'GO: Cellular Components', 
                              'GO: All Pathways')
                  )
                  ),
      pickerInput(inputId = "network",
                  label = "Select Interactions for Network Analysis:",
                  choices = c("STRING: Experimental & Database", "Advanced Options")
      ),
      conditionalPanel(condition = "input.network == 'Advanced Options'",
                       checkboxGroupInput(inputId ="STRING_interaction_sources", 
                                          label = "Select One or More Interaction Sources:",
                                          choices = c(Neighborhood = "neighborhood",
                                                      Fusion = "fusion",
                                                      Cooccurence = "cooccurence",
                                                      Coexpression = "coexpression",
                                                      Experimental = "experimental",
                                                      Database = "database",
                                                      Textmining = "textmining"),
                                          inline = F
                                          , textOutput("txt"))
      ),
      selectInput(inputId = "interaction_confidence_cutoff",
                  label = "Interaction Confidence for Network Analysis:",
                  choices = c("Low (>0.15)" = 150, 
                              "Medium (>0.4)" = 400,
                              "High (>0.7)" = 700),
                  selected = 400
      ),
      fileInput(inputId= "file1",
                label = 'Choose an input file to upload',
                buttonLabel = "Browse...",
                # Restrict input file types to .txt and .csv files
                accept=c("txt/csv", "text/comma-separated-values,text/plain", ".csv")
      ),
      # cutoff values depending the cutoff method chosen
      uiOutput("cutoffTypes"),
      # textInput("cutoff_valueH", "High Confidence Cutoff Value", placeholder = "High-conf cutoff"),
      textInput("cutoff_valueH", "High Confidence Cutoff Value"),
      bsPopover("cutoff_valueH", "High confidence cutoff value:", "Please enter a value for high confience cutoff, use \"-\" sign for negative value", placement = "bottom", trigger = "hover", options = NULL),
      # textInput("cutoff_valueM", "Med-conf Cutoff Value", placeholder = "Med-conf cutoff"),
      textInput("cutoff_valueM", "Medium Confidence Cutoff Value"),
      bsPopover("cutoff_valueM", "Medium confidence cutoff value:", "Please enter a value for medium confience cutoff, use \"-\" sign for negative value", placement = "bottom", trigger = "hover", options = NULL),
      checkboxInput("secondaryCutoff", "Add an Additional Criteria"),
      uiOutput("secondaryColumn_choice"),
      conditionalPanel(
        condition = "input.secondaryCutoff == 1",
        fluidRow(
          column(3,
                 selectInput("secondaryDirection", "Direction:", c('\u2265','\u2264'), selected='\u2265')
          ),
          column(6,
                 textInput(inputId="secondaryValue", label="Value")
          )
        )
      ),
      bsPopover("secondaryCutoff", "Use the values of another column to set a cutoff that all medium and high confidence hits must meet", placement = "bottom", trigger = "hover", options = NULL),
      checkboxInput("includeBackground", "Add Genome Background"),
      bsPopover("includeBackground", "To include known coding genes that are not on your input gene list as background", placement = "bottom", trigger = "hover", options = NULL),          
      actionButton("goButton", "Analyze my data",
                   style="padding:4px; font-size:120%; color: #fff; background-color: rgb(1, 81, 154); border-color: #2e6da4"),
      actionButton("refresh", "Reset", icon("undo"),
                   style="padding:4px; font-size:120%; color: #fff; background-color: rgb(1, 81, 154); border-color: #2e6da4"),
      width = 3
    ),
    # Show a plot of the generated distribution
    mainPanel(
      tabsetPanel(id = "inTabset",
                  tabPanel(title = "Home", value = "landingPage",
                           htmlOutput("spacer1"),
                           h3("The SIGNAL platform is designed to facilitate robust hit selection from normalized high-throughput datasets.", align = "center"),
                           HTML('<center><img src="images/SIGNAL_schmematic_Home.png" width="1200"></center>'),
                           br(),
                           h4("For a detailed user guide and sample dataset click on the “User Guide” tab.", align = "center", style = "color:grey"),
                           br(),
                           br(),
                           br(),
                           HTML("<center><p>Our manuscript describing the design, development, and testing of SIGNAL can be found <a href='https://www.cell.com/cell-systems/fulltext/S2405-4712(21)00081-8'>here</a>.</p></center>"),
                           HTML("<center><p>To cite: Katz, S., Song, J., Webb, K. P., Lounsbury, N. W., Bryant, C. E., &amp; Fraser, I. D. C. (2021). SIGNAL: A web-based iterative analysis platform integrating pathway and network approaches optimizes hit selection from genome-scale assays. <em>Cell Systems, 12</em>(4), 338-352.e335. <a href='doi:https://doi.org/10.1016/j.cels.2021.03.001'>doi:https://doi.org/10.1016/j.cels.2021.03.001</a></p></center>")
                           #img(src = "images/SIGNAL_schmematic_Home.png", height="75%", width="75%", align = "center")
                  ),
                  tabPanel(title = "Input", value = "contents",
                           htmlOutput("spacer2"),
                           # dataTableOutput("contents"),
                           # br(),
                           # br(),
                           # fluidRow(div(style='display:inline-block; padding-left:0px; padding-right:0px; horizontal-align:top;',
                           #   column(width =  6, 
                           #          #div(style='display:inline-block; padding-left:0px; padding-right:0px; horizontal-align:top;',
                           #          textAreaInput(inputId="manualHits", label="Optional: Enter gene IDs that should be kept \n as high confidence hits throughout the analysis:", width="100%", height="200px", value="", placeholder = "7099, 4615, 51135")
                           #   ),
                           #   column(width = 3,
                           #          div(style=' padding-top:50px;',
                           #          radioButtons("manualIDtypes", h5("Gene ID type:"),
                           #                           choices = list("Entrez ID" = "EntrezID", "Gene Symbol" = "GeneSymbol"), selected = "EntrezID")
                           #          ))))
                           # ,
                           
                           tabsetPanel(id = 'contents',
                                       # Input genes table
                                       tabPanel(title = "Mapped Genes", value = "contents",
                                                dataTableOutput("contents"),
                                                br(),
                                                br(),
                                                fluidRow(div(style='display:inline-block; padding-left:0px; padding-right:0px; horizontal-align:top;',
                                                             column(width =  6,
                                                                    #div(style='display:inline-block; padding-left:0px; padding-right:0px; horizontal-align:top;',
                                                                    textAreaInput(inputId="manualHits", label="Optional: Enter gene IDs that should be kept \n as high confidence hits throughout the analysis:", width="100%", height="200px", value="", placeholder = "7099, 4615, 51135")
                                                             ),
                                                             column(width = 3,
                                                                    div(style=' padding-top:50px;',
                                                                        radioButtons("manualIDtypes", h5("Gene ID type:"),
                                                                                     choices = list("Entrez ID" = "EntrezID", "Gene Symbol" = "GeneSymbol"), selected = "EntrezID")
                                                                    ))))
                                       ),
                                       tabPanel(title = "Unmapped Genes", value = "contents2",
                                                dataTableOutput("contents2")
                                       )
                           )
                           
                           
                           
                           #textAreaInput(inputId="manualHits", label="Optional: Enter gene IDs that should be kept \n as high confidence hits throughout the analysis:", width="100%", height="200px", value="", placeholder = "7099, 4615, 51135"),
                           # div(style="display:grid",textAreaInput(inputId="manualHits", label="Optional: Enter gene IDs that should be kept \n as high confidence hits throughout the analysis:", width="100%", height="200px", value="", placeholder = "7099, 4615, 51135")),
                           # div(style="display:grid",radioButtons("manualIDtypes", h5("Gene ID type:"),
                           #                                               choices = list("Entrez ID" = "EntrezID", "Gene Symbol" = "GeneSymbol"), selected = "EntrezID")),
                           # #div(style="display: inline-block;vertical-align:top; width: 150px;",selectInput("ddllgra", "Function:",c('mean','median','sd','count','min','max'), selected='mean')),
                           # #div(style="display: inline-block;vertical-align:top; width: 150px;",textInput(inputId="xlimitsmax", label="x-max", value = 0.5))),
                  ),
                  tabPanel(title = "Enriched Pathways", value = "enrichedPathways",
                           htmlOutput("spacer3"),
                           dataTableOutput("enrichedPathways"),
                           htmlOutput("notes")
                  ),
                  tabPanel(title = "Gene Hits", value = "signalHits",
                           htmlOutput("spacer4"),
                           tabsetPanel(id = 'signalHits',
                                       #SIGNAL Hits table
                                       tabPanel(title = "SIGNAL Gene Hits", value = "signalHits",
                                                dataTableOutput("signalHits")),
                                       tabPanel(title = "Gene Hits By Iteration", value="geneList",
                                                dataTableOutput("geneList")),
                                       tabPanel(title="Graph: Gene Hits By Iteration", value="geneHitsByIteration",
                                                dataTableOutput("geneHitsTableByIteration"),
                                                plotOutput("geneHitsByIteration")),
                                       tabPanel(title="High Confidence Hits Not Selected by SIGNAL", value="nonSIGNALhits",
                                                dataTableOutput("nonSIGNALhitsTable")),
                                       tabPanel(title = "Pathway Enrichments", value = "pathwayEnrich.cond",
                                                dataTableOutput("pathwayEnrich.cond"))
                           )
                  ),
                  tabPanel(title = "Network", value = "myNetworkGraph",
                           h4('"Select up to 3 pathways for network graph analysis'), hr(),
                           div(style="display:inline-block",textInput(inputId="mySelection", label="Your selected pathway IDs", value = 0.0)),
                           div(style="display:inline-block",uiOutput("submitGraph")),
                           div(style="display:inline-block",uiOutput("link2Graph")),
                           dataTableOutput("myNetworkGraph")
                  ),
                  tabPanel(title = "Network Graph", value = "graphViews",
                           htmlOutput("spacer5"),
                           tabsetPanel(id = 'igraphViews',
                                       # Display in igraph
                                       tabPanel(title="1st Degree Network", value="graphView1",
                                                HTML("<div id='graphView1'></div>"),
                                                #htmlOutput("graphLegend1"),
                                                htmlOutput("graphView1i", width = "100%", height = "700px")
                                       ),
                                       tabPanel(title="2nd Degree Network", value="graphView2",
                                                HTML("<div id='graphView2'></div>"),
                                                #htmlOutput("graphLegend2"),
                                                htmlOutput("graphView2i", width = "100%", height = "700px")
                                       ),
                                       tabPanel(title = "Network Table", value = "PathNetTable",
                                                HTML("<div id='PathNetTable'></div>"),
                                                dataTableOutput("PathNetTable")
                                       ),
                                       tabPanel(title = "Clicked Pathways Table", value = "ClickedDataTable",
                                                HTML("<div id='ClickedTable'></div>"),
                                                dataTableOutput("ClickedDataTable")
                                       )
                           )
                  ),
                  tabPanel(title = "Download", value = "downloads",
                           htmlOutput("spacer7"),
                           htmlOutput("downloadFiles"),
                           downloadButton('downloadButton', 'Download all files')
                  ),
                  tabPanel(title = "User Guide", value = "userGuide",
                           htmlOutput("spacer8"),
                           navlistPanel(id = 'userGuideTab', widths = c(3, 6),
                                        tabPanel(title = "1. Sample dataset", value = "sampleDataset",
                                                 uiOutput("sampledataset_page")),
                                        tabPanel(title = "2. Preparing your data", value = "dataPrep",
                                                 uiOutput("dataPrep_page")),
                                        tabPanel(title = "3. Running an analysis", value = "runAnalysis",
                                                 uiOutput("runAnalysis_page")),
                                        tabPanel(title = "4. Reading SIGNAL results", value = "readResults",
                                                 uiOutput("readResults_page")),
                                        tabPanel(title = "5. Saving and securing your analysis", value = "saveAnalysis",
                                                 uiOutput("saveAnalysis_page")),
                                        tabPanel(title = "Appendix A: How to segment data into confidence tiers", value = "segmentData",
                                                 uiOutput("segmentData_page")),
                                        tabPanel(title = "Appendix B: Data security", value = "security",
                                                 uiOutput("dataSecurity_page")),
                                        tabPanel(title = "Appendix C: Running SIGNAL with alternative or bespoke databases", value = "SIGNALcode",
                                                 uiOutput("standaloneR"))
                           )
                  ),
                  
                  tabPanel(title = "Help", value = "helpUs",
                           htmlOutput("spacer9"),
                           tabsetPanel(id = 'helpTab',
                                       tabPanel(title = "Contact us", value = "contactUS",
                                                uiOutput("contactUS")),
                                       tabPanel(title = "Documentation", value = "readMe",
                                                uiOutput("documentation")),        
                                       tabPanel(title = "Updates", value = "changeLog",
                                                uiOutput("changeLog"))                          
                           )
                  )
      ),
      width = 9
    )
  ),
  
  # Show a footer using the header style
  headerPanel(includeHTML("footer.html"))
)

##################################################