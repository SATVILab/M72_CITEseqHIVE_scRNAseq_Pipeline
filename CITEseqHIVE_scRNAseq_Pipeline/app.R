source('global.R')

addResourcePath("docs", "_book")

ui <- dashboardPage(
  dashboardHeader(title = "M72 Immune Correlates scRNAseq Pilot"),
  dashboardSidebar(
    sidebarMenu(id='tab',
                useShinyjs(),
                menuItem("Home Page", tabName = "home", icon=icon("list")),
                menuItem("scRNAseq Analyzer", tabName = "input", icon = icon("edit")),
                conditionalPanel(condition = "input.tab == 'input'",
                                 div(
                                   fileInput("file", "Upload File", multiple = TRUE, accept = c('.rds', '.Rdata')),
                                   actionButton("reset", "Reset", icon = icon("undo"), style = "color: #fff; background-color: #dc3545; width: 87.25%"),
                                   actionButton("run", "Run", icon = icon("play"), style = "color: #fff; background-color: #28a745; width: 87.25%")
                                 )
                )
    )
  ),
  dashboardBody(
    tabItems(
      tabItem(tabName = "input",
              tabsetPanel(id = 'main_tabs',
                          tabPanel("Book",
                                   tags$iframe(src = "docs/intro.html", style = "width: 100%; height: 1600px; border: none;"))
              )
      ),
      tabItem(tabName = "home",
              tags$h1(HTML("Welcome to our app")))
    )
  )
)

server <- function (input, output, session){
  options(shiny.maxRequestSize=1000*1024^2)
  
  obj_D0_reactive <- reactiveVal(NULL)
  obj_D37_reactive <- reactiveVal(NULL)
  obj_M72_reactive <- reactiveVal(NULL)
  obj_UNSTIM_reactive <- reactiveVal(NULL)
  
  obj_D0_Lymphoid_reactive <- reactiveVal(NULL)
  obj_D37_Lymphoid_reactive <- reactiveVal(NULL)
  obj_M72_Lymphoid_reactive <- reactiveVal(NULL)
  obj_UNSTIM_Lymphoid_reactive <- reactiveVal(NULL)
  obj_Lymphoid_reactive <- reactiveVal(NULL)
  
  obj_D0_Myeloid_reactive <- reactiveVal(NULL)
  obj_D37_Myeloid_reactive <- reactiveVal(NULL)
  obj_M72_Myeloid_reactive <- reactiveVal(NULL)
  obj_UNSTIM_Myeloid_reactive <- reactiveVal(NULL)
  obj_Myeloid_reactive <- reactiveVal(NULL)
  
  obj_D0_T_reactive <- reactiveVal(NULL)
  obj_D37_T_reactive <- reactiveVal(NULL)
  obj_M72_T_reactive <- reactiveVal(NULL)
  obj_UNSTIM_T_reactive <- reactiveVal(NULL)
  obj_T_reactive <- reactiveVal(NULL)
  
  shinyjs::disable('run')
  
  observe({
    if (!is.null(input$file)){
      shinyjs::enable('run')
      print(input$file)
    } else {
      shinyjs::disable('run')
    }
  })
  
  observeEvent(input$reset, {
    shinyjs::reset('file')
    shinyjs::disable('run')
  })
  
  observeEvent(input$run, {
    shinyjs::disable('run')
    
    show_modal_spinner(text = "Preparing plots...")
    obj <- load_seurat_obj(input$file$datapath)
    
    if (is.vector(obj)){
      showModal(modalDialog(
        title = 'Error with file',
        HTML("<h5>There is an error with the file you uploaded. See below for more details.</h5><br>",
             paste(unlist(obj), collapse = "<br><br>"))
      ))
      shiny::enable('run')
    } else {
      
      # Generate subset objs
      obj_D0 <- subset(obj, subset = VISIT %in% c("D0"))
      obj_D37 <- subset(obj, subset = VISIT %in% c("D37"))
      obj_M72 <- subset(obj, subset = STIM %in% c("M72"))
      obj_UNSTIM <- subset(obj, subset = STIM %in% c("UNSTIM"))
      
      ## Progenitor cells and HSC/MPP cells left out
      obj_D0_Lymphoid <- subset(obj_D0, subset = sctype_classification %in% c("Pro-B cells", "Pre-B cells", "Naive B cells", "Memory B cells", "Plasma B cells", "Naive CD8+ T cells", "Naive CD4+ T cells", "Memory CD8+ T cells", "Memory CD4+ T cells", "Effector CD8+ T cells", "Effector CD4+ T cells", "CD8+ NKT-like cells", "CD4+ NKT-like cells", "Natural killer cells", "Antigen-specific CD8+ T cells", "Antigen-specific CD4+ T cells")) #add in gd-T cells but need to know how to put symbols
      obj_D37_Lymphoid <- subset(obj_D37, subset = sctype_classification %in% c("Pro-B cells", "Pre-B cells", "Naive B cells", "Memory B cells", "Plasma B cells", "Naive CD8+ T cells", "Naive CD4+ T cells", "Memory CD8+ T cells", "Memory CD4+ T cells", "Effector CD8+ T cells", "Effector CD4+ T cells", "CD8+ NKT-like cells", "CD4+ NKT-like cells", "Natural killer cells", "Antigen-specific CD8+ T cells", "Antigen-specific CD4+ T cells")) #add in gd-T cells but need to know how to put symbols
      obj_M72_Lymphoid <- subset(obj_M72, subset = sctype_classification %in% c("Pro-B cells", "Pre-B cells", "Naive B cells", "Memory B cells", "Plasma B cells", "Naive CD8+ T cells", "Naive CD4+ T cells", "Memory CD8+ T cells", "Memory CD4+ T cells", "Effector CD8+ T cells", "Effector CD4+ T cells", "CD8+ NKT-like cells", "CD4+ NKT-like cells", "Natural killer cells", "Antigen-specific CD8+ T cells", "Antigen-specific CD4+ T cells")) #add in gd-T cells but need to know how to put symbols
      obj_UNSTIM_Lymphoid <- subset(obj_UNSTIM, subset = sctype_classification %in% c("Pro-B cells", "Pre-B cells", "Naive B cells", "Memory B cells", "Plasma B cells", "Naive CD8+ T cells", "Naive CD4+ T cells", "Memory CD8+ T cells", "Memory CD4+ T cells", "Effector CD8+ T cells", "Effector CD4+ T cells", "CD8+ NKT-like cells", "CD4+ NKT-like cells", "Natural killer cells", "Antigen-specific CD8+ T cells", "Antigen-specific CD4+ T cells")) #add in gd-T cells but need to know how to put symbols
      obj_Lymphoid <- subset(obj, subset = sctype_classification %in% c("Pro-B cells", "Pre-B cells", "Naive B cells", "Memory B cells", "Plasma B cells", "Naive CD8+ T cells", "Naive CD4+ T cells", "Memory CD8+ T cells", "Memory CD4+ T cells", "Effector CD8+ T cells", "Effector CD4+ T cells", "CD8+ NKT-like cells", "CD4+ NKT-like cells", "Natural killer cells", "Antigen-specific CD8+ T cells", "Antigen-specific CD4+ T cells")) #add in gd-T cells but need to know how to put symbols
      
      obj_D0_Myeloid <- subset(obj_D0, subset = sctype_classification %in% c("Eosinophils", "Neutrophils", "Basophils", "Mast cells", "Classical Monocytes", "Non-classical monocytes", "Intermediate monocytes", "Macrophages", "Megakaryocyte", "Erythroid-like and erythroid precursor cells", "Myeloid Dendritic cells", "Plasmacytoid Dendritic cells", "Granulocytes", "ISG expressing immune cells"))
      obj_D37_Myeloid <- subset(obj_D37, subset = sctype_classification %in% c("Eosinophils", "Neutrophils", "Basophils", "Mast cells", "Classical Monocytes", "Non-classical monocytes", "Intermediate monocytes", "Macrophages", "Megakaryocyte", "Erythroid-like and erythroid precursor cells", "Myeloid Dendritic cells", "Plasmacytoid Dendritic cells", "Granulocytes", "ISG expressing immune cells"))
      obj_M72_Myeloid <- subset(obj_M72, subset = sctype_classification %in% c("Eosinophils", "Neutrophils", "Basophils", "Mast cells", "Classical Monocytes", "Non-classical monocytes", "Intermediate monocytes", "Macrophages", "Megakaryocyte", "Erythroid-like and erythroid precursor cells", "Myeloid Dendritic cells", "Plasmacytoid Dendritic cells", "Granulocytes", "ISG expressing immune cells"))
      obj_UNSTIM_Myeloid <- subset(obj_UNSTIM, subset = sctype_classification %in% c("Eosinophils", "Neutrophils", "Basophils", "Mast cells", "Classical Monocytes", "Non-classical monocytes", "Intermediate monocytes", "Macrophages", "Megakaryocyte", "Erythroid-like and erythroid precursor cells", "Myeloid Dendritic cells", "Plasmacytoid Dendritic cells", "Granulocytes", "ISG expressing immune cells"))
      obj_Myeloid <- subset(obj, subset = sctype_classification %in% c("Eosinophils", "Neutrophils", "Basophils", "Mast cells", "Classical Monocytes", "Non-classical monocytes", "Intermediate monocytes", "Macrophages", "Megakaryocyte", "Erythroid-like and erythroid precursor cells", "Myeloid Dendritic cells", "Plasmacytoid Dendritic cells", "Granulocytes", "ISG expressing immune cells"))
      
      obj_D0_T <- subset(obj_D0, subset = sctype_classification %in% c("Naive CD8+ T cells", "Naive CD4+ T cells", "Memory CD8+ T cells", "Memory CD4+ T cells", "Effector CD8+ T cells", "Effector CD4+ T cells", "CD8+ NKT-like cells", "CD4+ NKT-like cells", "Natural killer cells", "Antigen-specific CD8+ T cells", "Antigen-specific CD4+ T cells")) #add in gd-T cells but need to know how to put symbols
      obj_D37_T <- subset(obj_D37, subset = sctype_classification %in% c("Naive CD8+ T cells", "Naive CD4+ T cells", "Memory CD8+ T cells", "Memory CD4+ T cells", "Effector CD8+ T cells", "Effector CD4+ T cells", "CD8+ NKT-like cells", "CD4+ NKT-like cells", "Natural killer cells", "Antigen-specific CD8+ T cells", "Antigen-specific CD4+ T cells")) #add in gd-T cells but need to know how to put symbols
      obj_M72_T <- subset(obj_M72, subset = sctype_classification %in% c("Naive CD8+ T cells", "Naive CD4+ T cells", "Memory CD8+ T cells", "Memory CD4+ T cells", "Effector CD8+ T cells", "Effector CD4+ T cells", "CD8+ NKT-like cells", "CD4+ NKT-like cells", "Natural killer cells", "Antigen-specific CD8+ T cells", "Antigen-specific CD4+ T cells")) #add in gd-T cells but need to know how to put symbols
      obj_UNSTIM_T <- subset(obj_UNSTIM, subset = sctype_classification %in% c("Naive CD8+ T cells", "Naive CD4+ T cells", "Memory CD8+ T cells", "Memory CD4+ T cells", "Effector CD8+ T cells", "Effector CD4+ T cells", "CD8+ NKT-like cells", "CD4+ NKT-like cells", "Natural killer cells", "Antigen-specific CD8+ T cells", "Antigen-specific CD4+ T cells")) #add in gd-T cells but need to know how to put symbols
      obj_T <- subset(obj, subset = sctype_classification %in% c("Naive CD8+ T cells", "Naive CD4+ T cells", "Memory CD8+ T cells", "Memory CD4+ T cells", "Effector CD8+ T cells", "Effector CD4+ T cells", "CD8+ NKT-like cells", "CD4+ NKT-like cells", "Natural killer cells", "Antigen-specific CD8+ T cells", "Antigen-specific CD4+ T cells")) #add in gd-T cells but need to know how to put symbols
     
      # Make subset objs reactive
      obj_D0_reactive(obj_D0)
      obj_D37_reactive(obj_D37)
      obj_M72_reactive(obj_M72)
      obj_UNSTIM_reactive(obj_UNSTIM)
      
      obj_D0_Lymphoid_reactive(obj_D0_Lymphoid)
      obj_D37_Lymphoid_reactive(obj_D37_Lymphoid)
      obj_M72_Lymphoid_reactive(obj_M72_Lymphoid)
      obj_UNSTIM_Lymphoid_reactive(obj_UNSTIM_Lymphoid)
      obj_Lymphoid_reactive(obj_Lymphoid)
      
      obj_D0_Myeloid_reactive(obj_D0_Myeloid)
      obj_D37_Myeloid_reactive(obj_D37_Myeloid)
      obj_M72_Myeloid_reactive(obj_M72_Myeloid)
      obj_UNSTIM_Myeloid_reactive(obj_UNSTIM_Myeloid)
      obj_Myeloid_reactive(obj_Myeloid)
      
      obj_D0_T_reactive(obj_D0_T)
      obj_D37_T_reactive(obj_D37_T)
      obj_M72_T_reactive(obj_M72_T)
      obj_UNSTIM_T_reactive(obj_UNSTIM_T)
      obj_T_reactive(obj_T)
      
      insertTab(inputId = 'main_tabs',
                tabPanel("UMAP",
                         fluidRow(
                           column(
                             width = 8,
                             plotOutput(outputId = 'umap'),
                             downloadButton('download_umap', "Download UMAP")
                           ),
                           column(
                             width = 4,
                             selectizeInput(
                               'metadata_col',
                               "Metadata Column",
                               colnames(obj@meta.data)
                             )
                           )
                         ),
                         style = "height: 90%; width: 95%; padding-top: 5%;"),
                select = TRUE)
      insertTab(inputId = 'main_tabs',
                tabPanel("Compare UMAPs",
                         fluidRow(
                           column(
                             width = 8,
                             plotOutput(outputId = 'compare.umaps', height = "800px"),
                             downloadButton('download_compare.umaps', "Download UMAPs")
                           ),
                           column(
                             width = 4,
                             selectizeInput(
                               'group.by',
                               "Group by",
                               choices = c("HIVE" = "HIVE", "Participant ID" = "PID", "Visit" = "VISIT", "Stimulation" = "STIM", "Visit and Stimulation" = "VISIT_STIM", "Vaccination Status" = "group"),
                               selected = "STIM"
                             ),
                             selectizeInput(
                               'split.by',
                               "Split by",
                               choices = c("HIVE" = "HIVE", "Participant ID" = "PID", "Visit" = "VISIT", "Stimulation" = "STIM", "Visit and Stimulation" = "VISIT_STIM", "Vaccination Status" = "group"),
                               selected = "VISIT"
                             )
                           )
                         ),
                         style = "height: 90%; width: 95%; padding-top: 5%;"),
                select = TRUE)
      insertTab(inputId = 'main_tabs',
                tabPanel("Feature Plots",
                         fluidRow(
                           column(
                             width = 8,
                             plotOutput(outputId = 'featurePlot', height = "800px"),
                             downloadButton('download_featurePlot', "Download Feature Plot")
                           ),
                           column(
                             width = 4,
                             selectizeInput(
                               'gene',
                               "Genes",
                               rownames(obj)
                             )
                           )
                         ),
                         style = "height: 90%; width: 95%; padding-top: 5%;"),
                select = TRUE)
      insertTab(inputId = 'main_tabs',
                tabPanel("Violin Plots - Antigen Specific Genes",
                         fluidRow(
                           column(
                             width = 8,
                             plotOutput(outputId = 'AS.vlnplots_RNA', height = "1200px"),
                             downloadButton('download_AS.vlnplots_RNA', "Download Violin Plots")
                           ),
                           column(
                             width = 4,
                             selectizeInput(
                               'AS.gene',
                               "Genes",
                               AS.genes_RNA
                             ),
                             selectizeInput(
                               'layer.gene',
                               "Layer",
                               choices = c("Relative Expression" = "data", "Counts" = "counts"),
                               selected = "data"
                             )
                           )
                         ),
                         style = "height: 90%; width: 95%; padding-top: 5%;"),
                select = TRUE)
      insertTab(inputId = 'main_tabs',
                tabPanel("Violin Plots - Antigen Specific Proteins",
                         fluidRow(
                           column(
                             width = 8,
                             plotOutput(outputId = 'AS.vlnplots_ADT', height = "1200px"),
                             downloadButton('download_AS.vlnplots_ADT', "Download Violin Plots")
                           ),
                           column(
                             width = 4,
                             selectizeInput(
                               'AS.protein',
                               "Proteins",
                               AS.genes_ADT
                             ),
                             selectizeInput(
                               'layer.protein',
                               "Layer",
                               choices = c("Relative Expression" = "data", "Counts" = "counts"),
                               selected = "data"
                             )
                           )
                         ),
                         style = "height: 90%; width: 95%; padding-top: 5%;"),
                select = TRUE)
      insertTab(inputId = 'main_tabs',
                tabPanel("Feature Plots - Antigen Specific Genes",
                         fluidRow(
                           column(
                             width = 8,
                             plotOutput(outputId = 'AS.featureplots_RNA', height = "1200px"),
                             downloadButton('download_AS.featureplots_RNA', "Download Feature Plots")
                           ),
                           column(
                             width = 4,
                             selectizeInput(
                               'AS.gene',
                               "Genes",
                               AS.genes_RNA
                             )
                           )
                         ),
                         style = "height: 90%; width: 95%; padding-top: 5%;"),
                select = TRUE)
      insertTab(inputId = 'main_tabs',
                tabPanel("Feature Plots - Antigen Specific Proteins",
                         fluidRow(
                           column(
                             width = 8,
                             plotOutput(outputId = 'AS.featureplots_ADT', height = "1200px"),
                             downloadButton('download_AS.featureplots_ADT', "Download Feature Plots")
                           ),
                           column(
                             width = 4,
                             selectizeInput(
                               'AS.protein',
                               "Proteins",
                               AS.genes_ADT
                             )
                           )
                         ),
                         style = "height: 90%; width: 95%; padding-top: 5%;"),
                select = TRUE)
       
      remove_modal_spinner()
      shinyjs::enable('run')
      
      output$umap <- renderPlot({
        create_metadata_UMAP(obj, input$metadata_col)
      })
      output$compare.umaps <- renderPlot({
        compare.umaps(obj, input$group.by, input$split.by)
      })
      output$featurePlot <- renderPlot({
        create_feature_plot(obj, input$gene)
      })
      output$AS.vlnplots_RNA <- renderPlot({
        req(obj_T_reactive(), obj_M72_T_reactive(), obj_UNSTIM_T_reactive(), input$AS.gene, input$layer.gene)
        cat("Rendering plot for gene: ", input$AS.gene, " with layer: ", input$layer.gene, "\n")
        create_AS.vlnplots_RNA(obj_T_reactive(), obj_M72_T_reactive(), obj_UNSTIM_T_reactive(), input$AS.gene, input$layer.gene)
      })
      output$AS.vlnplots_ADT <- renderPlot({
        req(obj_T_reactive(), obj_M72_T_reactive(), obj_UNSTIM_T_reactive(), input$AS.protein, input$layer.protein)
        cat("Rendering plot for protein: ", input$AS.protein, " with layer: ", input$layer.protein, "\n")
        create_AS.vlnplots_ADT(obj_T_reactive(), obj_M72_T_reactive(), obj_UNSTIM_T_reactive(), input$AS.protein, input$layer.protein)
      })
      output$AS.featureplots_RNA <- renderPlot({
        req(obj_T_reactive(), input$AS.gene)
        cat("Rendering plot for gene: ", input$AS.gene, "\n")
        create_AS.featureplots_RNA(obj_T_reactive(), input$AS.gene)
      })
      output$AS.featureplots_ADT <- renderPlot({
        req(obj_T_reactive(), input$AS.protein)
        cat("Rendering plot for protein: ", input$AS.protein, "\n")
        create_AS.featureplots_ADT(obj_T_reactive(), input$AS.protein)
      })
    }
  })
}

shinyApp(ui, server)
