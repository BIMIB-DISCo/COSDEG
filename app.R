library(shiny)
library(shinydashboard)
library(Seurat)
library(shinyFiles)

source("chooser.R")



header <- dashboardHeader(
  title = "COSDEG"
)

body <- dashboardBody(
  fluidRow(
    column(width = 3,
           box(width = NULL, status = "warning",
               title = "Uploading Files",
               
               # Input: Select a file ----
               shinyDirButton('dataset_dir', label='Dir select', title='Please select a dataset', multiple=FALSE),
               textOutput("dataset_dir"),
               
               # Horizontal line ----
               tags$hr(),
               
               
               # Input: Select separator ----
               radioButtons("dataset_type", "Dataset type:",
                            choices = c("10X genomics mtx(.gz) matrix, barcodes and features" = "10X_mtx",
                              "10X hdf5" = "10X_h5",
                              "STARsolo" = "star",
                              "mtx file, barcodes file and features file " = "mtx"),
                            selected = "10X_mtx"
               ),

               # Horizontal line ----
               tags$hr(),
               
               textInput("dataset_name", "Dataset name", placeholder = "Choose a name"),
               
               ## Input: Select number of rows to display ----
               #radioButtons("disp", "Display", choices = c(Head = "head", All = "all"), selected = "head"),
               #p(
               # class = "text-muted",
                # paste("Note: a route number can have several different trips, each",
                 #      "with a different path. Only the most commonly-used path will",
                #       "be displayed on the map."
                # )
               #),
               actionButton("load_dataset", "Load dataset")
           )
    ),
    column(width = 9,
           box(width = NULL, solidHeader = FALSE, title = "Select datasets",
               status = "success",
             uiOutput("mydatasets_out")  
           ),
           box(width = NULL, collapsible = TRUE, collapsed = TRUE,
               verbatimTextOutput("selection")
           ),
           box(width = NULL, solidHeader = FALSE, title = "Project",
               status = "info",
               textInput("project_name", "Project name", placeholder = "Choose a name"),
               actionButton("create_project", "Create project")
           ),
           
           box(width = NULL,
               collapsible = TRUE, collapsed = TRUE,
               # Output: Data file ----
               verbatimTextOutput("contents")
           )
    )
  )
)


# check is a dataset folder
check_dataset_dir <- function(folder) {
  if(dir.exists(folder)) {
    mtx <- list.files(path = folder, 
               pattern = ".*\\.mtx|.*\\.mtx\\.gz|.*\\.h5", 
               full.names = TRUE, include.dirs = FALSE, 
               no.. = TRUE, recursive = TRUE
              )
    if(length(mtx)>1) {
      return(check_multidatasets_dir(folder))
    }
    else if (length(mtx) == 0) {
      return(NULL)
    }
    else {
      mtx_dir <- dirname(mtx[[1]])
      mtx_files <- list.files(path = folder, 
                        pattern = ".*\\.mtx|.*\\.mtx\\.gz|.*\\.h5|.*\\.tsv|.*\\.tsv\\.gz", 
                        full.names = TRUE, include.dirs = FALSE, 
                        no.. = TRUE, recursive = TRUE
      )
      mtx_files <- lapply(mtx_files, basename)
      
    }
      
    
  }
}

# check is multi-datasets folder

# Define UI for data upload app ----
ui <- fluidPage(
  dashboardPage(
    header,
    dashboardSidebar(disable = TRUE),
    body
  )
)

# Define server logic to read selected file ----
server <- function(input, output) {
  roots=c(wd='.', root=.Platform$OS.type)
  dataset_names <- reactiveVal(value = list())
  dataset_names_right <- reactiveVal(value = list())  
  
  seurat_objects <- reactiveVal(value = list())
  
  
  shinyDirChoose(input, 'dataset_dir', roots=roots, filetypes=c('', 'mtx'))
  output$dataset_dir <- renderText(
    tryCatch(
      {
        dir_name <- paste0(
          c(roots[input$dataset_dir[[2]]],input$dataset_dir[[1]][-1]), 
          collapse = .Platform$file.sep
        )
      },
      error = function(e) {
        # return a safeError if a parsing error occurs
        NULL
        #stop(safeError(e))
      }
    )
  )
  
  output$contents <- renderPrint({
    req(input$dataset_dir)
    seurat_objects()
  })
  
  observeEvent(input$load_dataset, {
    req(input$dataset_dir)
    
    tryCatch(
      {
        dir_name <- paste0(
          c(roots[input$dataset_dir[[2]]],input$dataset_dir[[1]][-1]), 
          collapse = .Platform$file.sep
        )
        print("reading...")
        expression_matrix <- Read10X(data.dir = dir_name)
        seurat_object <- CreateSeuratObject(counts = expression_matrix)
        
        dataset_name <- input$dataset_name
        if(stringr::str_length(dataset_name) == 0)
          dataset_name <- "NoName"
        
        Misc(object = seurat_object, slot = "name") <- dataset_name
        seurat_objects( c(seurat_objects(), seurat_object) )
        tmp_list <- seurat_objects()
        names(tmp_list)[[length(tmp_list)]] <- dataset_name
        seurat_objects(tmp_list)
        print(seurat_object)
        print(seurat_objects())
        print("reading complete")
        
      },
      error = function(e) {
        # return a safeError if a parsing error occurs
        #stop(safeError(e))
        print("error: the folder has no valid dataset")
        NULL
      }
    )
    
    
    tmp_list <- c(input$mydatasets$left, list(dataset_name))
    tmp_list <- as.list(make.unique( unlist( c(tmp_list, input$mydatasets$right) ) ))[1:length(tmp_list)]

    dataset_names_right(input$mydatasets$right)
    dataset_names(tmp_list)
  })
  
  output$mydatasets_out <- renderUI({
    #req(input$dataset_dir)
    #req(dataset_names())
    req(input$load_dataset)
    chooserInput("mydatasets", "Available datasets", "Selected datasets",
                 dataset_names(), dataset_names_right(), size = 10, multiple = TRUE
    )
  })
  
  output$selection <- renderPrint(
    input$mydatasets
  )
  
}

shinyApp(ui, server)


