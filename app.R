### app.R

library(shiny)
library(shinyjs)
library(shinydashboard)
library(Seurat)
library(shinyFiles)
library(ggplot2)
library(dplyr)
library(DT)

source("cosdeg.R")

source("chooser.R")

source("project_mgr.R")

header <- dashboardHeader(
  title = "COSDEG"
)


last_created_project_id <- NULL

sidebar <- dashboardSidebar(
  sidebarMenu(
    menuItem("Create Project", tabName = "Load Dataset"),
    shinySaveButton("save_project", label = "Save", title = "Save project", multiple = FALSE),
    shinyFilesButton("load_project", label = "Load", title = "Load project", multiple = FALSE)
  ), disable = FALSE, collapsed = TRUE
)

body <- dashboardBody(
  tags$script(HTML("$('body').addClass('fixed');")),
  fluidRow(
    column(width = 4,
      box(width = NULL, status = "warning",
        title = "Recently Open Projects",
        uiOutput("recent_projects_list", inline = FALSE),
        br()
      )
    ),
    column(width = 8,
      box(width = NULL, solidHeader = FALSE, title = "New Project",
          status = "info",
          textInput("project_name", "Project name", placeholder = "Choose a name"),
          br()
          #actionButton("create_project", "Create project")
      )
    )
  ),
  fluidRow(
    column(width = 4,
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
               
               actionButton("load_dataset", "Load dataset")
           )
    ),
    column(width = 8,
           box(width = NULL, solidHeader = FALSE, title = "Select datasets",
               status = "success",
               uiOutput("mydatasets_out")  
           ),
          # box(width = NULL, collapsible = TRUE, collapsed = TRUE,
          #     verbatimTextOutput("selection")
          # ),
    
           box(width = NULL,
               collapsible = TRUE, collapsed = TRUE,
               # Output: Data file ----
               verbatimTextOutput("contents")
           ),
           box(width = NULL, solidHeader = FALSE,
               status = "danger",
               div(class = "pull-right",
                actionButton("create_project", "Create project")
               )
           )
    )
  ),
  fluidRow(
    column(width = 7
    ),
    column(width = 3,
          # box(width = NULL, solidHeader = FALSE,
          #     status = "danger",
          #     actionButton("create_project", "Create project")
          # )
    )
  )
  
)


# check is a dataset folder
check_dataset_dir <- function(folder) {
  if (is.null(folder))
    return(NULL)
  if(dir.exists(folder)) {
    mtx <- list.files(path = folder, 
                      pattern = ".*\\.mtx|.*\\.mtx\\.gz|.*\\.h5", 
                      full.names = TRUE, include.dirs = FALSE, 
                      no.. = TRUE, recursive = TRUE
    )
    if(length(mtx)>1) {
      print("Found multiple datasets")
      mtx_files <- lapply(mtx, basename)
      test_fun <- function(x){stringi::stri_detect_regex(mtx_files, pattern=x)}
      test_string <- c(".h5")
      res_test <- sapply(test_string, test_fun)
      if (any(res_test)) { # if h5 add files and check rest
        #print("test")
        #print(mtx[res_test])
        #print("anti test")
        #print(mtx[!res_test])
        #browser()
        return(c(mtx[res_test], check_multidatasets_dir( dirname(mtx[!res_test])) ) )        
      }
      return(check_multidatasets_dir(folder))
    }
    else if (length(mtx) == 0) {
      print("No datasets found")
      return(NULL)
    }
    else {
      print("Found 1 dataset")
      mtx_dir <- dirname(mtx[[1]])
      mtx_files <- list.files(path = mtx_dir, 
                              pattern = ".*\\.mtx|.*\\.mtx\\.gz|.*\\.h5|.*\\.tsv|.*\\.tsv\\.gz", 
                              full.names = TRUE, include.dirs = FALSE, 
                              no.. = TRUE, recursive = TRUE
      )
      mtx_basefiles <- lapply(mtx_files, basename)
      test_fun <- function(x){any(stringi::stri_detect_regex(mtx_basefiles, pattern=x))}
      test_string <- c("barcodes.tsv","features.tsv","matrix.mtx")
      if (all(sapply(test_string, test_fun))) {
        return(mtx_dir)
      } 
      test_string <- c(".h5")
      if (all(sapply(test_string, test_fun))) {
        return(mtx_files)
        #return(mtx_dir)        
      }
      return(NULL)
    }
  }
  return(NULL)
}

# check is multi-datasets folder
check_multidatasets_dir <- function(folder) {
  if (is.null(folder) || length(folder) == 0)
    return(NULL)
  
  if(length(folder)>1)
    return(lapply(folder, check_multidatasets_dir))
  
  #print(folder)
  if(dir.exists(folder)) {
    mtx <- list.files(path = folder, 
                      pattern = ".*\\.mtx|.*\\.mtx\\.gz|.*\\.h5", 
                      full.names = TRUE, include.dirs = FALSE, 
                      no.. = TRUE, recursive = TRUE
    )
    
    if(length(mtx)==0)
      return(NULL)
    
    #print(mtx)
    dirs <- list()
    file_dirs <- dirname(mtx)
    file_dirs <- unique(file_dirs)
    for (d in file_dirs) {
      #print(d)
      dirs <- c(dirs, check_dataset_dir(d))
    }
    #browser()
    if(length(dirs)>0) {
      return(dirs)
    }
    return(NULL)
  }
}



project_tab <- function(project_id) {
  ns <- NS(project_id)
  tab <- tabPanel(title = projects[[project_id]][["tab_id"]],
    fluidRow(
      column(width = 3,
        box(width = NULL, status = "warning",
          title = "Quality Check",
          
          selectizeInput(
            inputId = ns("qc_seurat_obj"), label = "Select datasets", 
            choices = names(projects[[project_id]][["seurat_objs"]]),
            multiple = TRUE, options = list(maxItems = 10)
          ),
          
          
          sliderInput(ns("perc_zeros"), label = "Zero genes per cell cutoff", 0.80, 1, value=0.95, step = 0.01),
          sliderInput(ns("sigma_mito"), label = "Std mitochondrial content cutoff", 0.50, 2, value=1.0, step = 0.5),
          sliderInput(ns("multiplet_rate"), label = "Doublets rate formation", 0.02, 0.2, value=0.04, step = 0.01),
          actionButton(ns("do_qc"), "Quality check")
        )
      ),
      column(width = 9,
        box(width = NULL, status = "warning",
          title = "Emptylets",
          #plotOutput("emptylets")
          #tableOutput("summary_qc")
          DT::dataTableOutput(ns("status")),
          DT::dataTableOutput(ns("summary_qc")),
          plotOutput(ns("emptylets_plot"))
        )
      )
    )
  )
  return(tab)
}
  

create_proj_gui <- function (project_name, selected_seurat_objs = NULL) {
  #browser()
  project_id <- create_proj(project_name, selected_seurat_objs)
  last_created_project_id <<- project_id 
  appendTab("tabs", project_tab(project_id), select = TRUE)
}




meaningful_name_path <- function(dataset_dir_list) {
  names_ <- xfun::sans_ext(dataset_dir_list)
  
  names_ <- stringr::str_split(unlist(names_),pattern = "/")
  dup_names_ <- unique(unlist(names_)[duplicated(unlist(names_))])
  unique(unlist(stringr::str_split(unlist(dataset_dir_list),pattern = "/")))
  known <- "data|dataset|raw|filtered|barcode|feature|matrix"
  pat <- stringi::stri_replace_all_fixed(paste0(c(dup_names_,known),collapse = "|"), ".", "\\.")  
  probable_names <- lapply(names_,function(x){x[!stringi::stri_detect_regex(x, pattern = pat)]})
  return(probable_names)
}



ui <- fluidPage(
  shinyjs::useShinyjs(),
  dashboardPage(
    header,
    sidebar,
    dashboardBody(
      navbarPage(id = "tabs", position = "static-top", #, "fixed-top"
        "",
        tabPanel("Load Data",
          body
        )         
      )
    )
  )
)



qc_server <- function(id, input_id) {
  moduleServer(
    id,
    function(input, output, session) {
      
      project <- reactiveVal(projects[[input_id]])
      #output$summary_qc <- NULL
      #output$emptylets_plot <- NULL
      
      output$summary_qc <- DT::renderDataTable(DT::datatable(project()[["summary_qc"]], options = list(dom = ''), rownames = TRUE))
      
      if(is.data.frame(project()[["status"]]))
        output$status <- DT::renderDataTable(DT::datatable(t(project()[["status"]]), options = list(dom = ''), rownames = TRUE))
      output$status <- DT::renderDataTable(DT::datatable(project()[["status"]], options = list(dom = ''), rownames = TRUE))
      
      output$emptylets_plot <- renderPlot(project()[["p"]])
      
      observeEvent(input$do_qc, {
        browser()
        summary_qc <- data.frame()
        s_data_status <- data.frame()
        
        if (!is.null(input$qc_seurat_obj))
          s_data <- projects[[input_id]][["seurat_objs"]][input$qc_seurat_obj]
        else
          s_data <- projects[[input_id]][["seurat_objs"]]
        
        result <- pre_qc(s_data, summary_qc)
        summary_qc <- result$summary_qc
        #s_data_status <- result$status
        
        result <- filter_offgenes(s_data, summary_qc, s_data_status)
        data.filt <- result$data.filt
        #projects[[input_id]][["data.filt"]] <- data.filt
        projects[[input_id]][["selected_f0"]] <<- result$selected_f0
        summary_qc <- result$summary_qc
        s_data_status <- result$status
        
        print(paste("isolate(input$perc_zeros)", isolate(input$perc_zeros)))
        result <- filter_emptylets(data.filt, summary_qc, isolate(input$perc_zeros), s_data_status)
        data.filt <- result$data.filt
        #projects[[input_id]][["data.filt"]] <- data.filt
        projects[[input_id]][["selected_cE"]] <<- result$selected_cE
        summary_qc <- result$summary_qc
        s_data_status <- result$status
        
        
        projects[[input_id]][["data.filt"]] <<- data.filt
        projects[[input_id]][["summary_qc"]] <<- summary_qc
        projects[[input_id]][["status"]] <<- s_data_status
        print(projects[[input_id]][["summary_qc"]])
        
        browser()
        df <- list()
        for (n in names(s_data)) {
          df[[n]] <- data.frame(nFeature_RNA = s_data[[n]]$nFeature_RNA, type = c( rep(n, dim(s_data[[n]])[[2]]))) 
        }
        df <- bind_rows(df)
        
        browser()
        
        projects[[input_id]][["p"]] <<- df %>% 
          ggplot( aes(x=nFeature_RNA, color=type, fill=type)) +
          geom_histogram(alpha=0.6, bins = 80) +
          geom_vline(xintercept = 3, linetype="dotted", color = "black", size=1) +
          xlab("") +
          ylab("Frequency") +
          facet_wrap(~type)
        
        project(projects[[input_id]])
        #
        #print(project()) #### ERROR because it contains plot handles
        #browser()
        
      })
      
      #project_reactive(projects[[input_id]])
      
      #browser()
      #return(list(project_id = reactive(input_id), project = reactive(projects[[input_id]])))
      return(project) ##not needed
    }
  )
}


# Define server logic to read selected file ----
server <- function(input, output) {
  roots = c(wd='.', root=.Platform$OS.type)
  
  user_cosdeg_path <- getwd()
  

  dataset_names <- reactiveVal(value = list())
  dataset_names_right <- reactiveVal(value = list())  
  
  seurat_objects <- reactiveVal(value = list())
  
  recent_projects_file_conf_ <- reactiveVal(file.path(user_cosdeg_path,"recent_projects.txt"))
  
  recent_projects_ <- reactiveValues(renderd=data.frame(projectname = c(),filename = c()))
  rvs = reactiveValues(buttons = list(), observers = list()) 
  
  observeEvent(recent_projects_, {
    if (file.exists(recent_projects_file_conf_())) {
      recent_projects_files <- read.csv(recent_projects_file_conf_(), header = TRUE, sep = ",", strip.white = TRUE)
      #browser()
      renderd <- rbind(isolate(recent_projects_$renderd), recent_projects_files)
      rownames(renderd) <- renderd$filename
      if (length(renderd)>5)
        renderd <- renderd[1:5,]
      recent_projects_$renderd <- renderd
      
    }
  }, once = TRUE)
  
  observeEvent(recent_projects_$renderd, {
    #browser()
    if (length(recent_projects_$renderd)==0)
      output$recent_projects_list <- renderText("No Projects")
    else {
      #create the list of recent projects
      rvs$buttons <- apply(recent_projects_$renderd, MARGIN = 1, function(i){
        flowLayout(actionLink(inputId = paste0(i["projectname"]), label = paste0(i["projectname"])))
      })
      
      output$recent_projects_list <- renderUI({
        do.call(shiny::tagList,as.list(rvs$buttons))
      })
      
      browser()
      #remove former button events 
      lapply(rvs$observers, function(i) i$destroy())
      
      #add new events to the buttons
      rvs$observers = apply(
        isolate(recent_projects_$renderd), MARGIN = 1,
        function(i) {
          observeEvent(input[[i["projectname"]]], {
            print(paste("Loading:", isolate(i[["filename"]])))
            browser()
            project_id <- load_project(isolate(i[["filename"]]), projects)
            appendTab("tabs", project_tab(project_id), select = TRUE)
            
            prj <- list(
              projectname = projects[[project_id]][["project_name"]], 
              filename =  normalizePath(projects[[project_id]][["file_name"]]))
            
            renderd <- isolate(recent_projects_$renderd)[which(rownames(isolate(recent_projects_$renderd))!=prj$filename),]
            renderd <- rbind(prj, renderd)
            rownames(renderd) <- renderd[["filename"]]
            
            #if (file.exists(recent_projects_file_conf_())) {
              write.table(renderd, file = isolate(recent_projects_file_conf_()), col.names =  TRUE, row.names = FALSE, sep = ",", quote = FALSE)
            #}
            
            print(project_id)
            browser()
            #aa <- qc_server(project_id, project_id)
            qc_server(project_id, project_id)
            #projects[[project_id]] <- aa()
            
            
          }, ignoreInit = TRUE) #end oberveevent
        } #end function
      )
    }
    #browser()
  })
  
  
  
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
    
    seurat_objects_ <- list()
    
    switch (input$dataset_type,
            "10X_mtx" = read_func <- function(x) {Read10X(data.dir = x)},
            "10X_h5" = read_func <- function(x) {Read10X_h5(filename = x)},
            "star" = read_func <- function(x) {ReadSTARsolo(data.dir = x)},
            read_func <- function(x) {ReadMtx(x)}
    )
    
    tryCatch(
      {
        dir_name <- paste0(
          c(roots[input$dataset_dir[[2]]],input$dataset_dir[[1]][-1]), 
          collapse = .Platform$file.sep
        )
        
        multidataset_dirs <- check_multidatasets_dir(dir_name)
        multidataset_names <- meaningful_name_path(multidataset_dirs)
        #browser()
        already_read <- unlist(lapply(seurat_objects(), function(x) x@misc$dir_name))
        
        for (i in 1:length(multidataset_dirs)) {
          skip <- stringi::stri_detect_regex(already_read, pattern=multidataset_dirs[[i]])
          if( length(skip)>0)
            if (any(skip)) {
              print("dataset already loaded ... continue")
              next
            }
          
          print("reading...")
          
          expression_matrix <- tryCatch(
            expression_matrix_ <- read_func(multidataset_dirs[[i]]),
            error = function(e) {
              print(paste("error using:", deparse1(read_func)))
              print("trying different read method...")
              
              if (dir.exists(multidataset_dirs[[i]])) {
                expression_matrix_ <- Read10X(data.dir = multidataset_dirs[[i]])
              }
              else {
                expression_matrix_ <- Read10X_h5(filename = multidataset_dirs[[i]])
              }
              expression_matrix_
            }
          )
          
          seurat_object <- CreateSeuratObject(counts = expression_matrix)
          
          dataset_name <- input$dataset_name
          if(stringr::str_length(dataset_name) == 0)
            dataset_name <- "NoName"
          if (length(multidataset_dirs) > 1) {
            if(length(multidataset_names[[i]]) > 0)
              dataset_name <- multidataset_names[[i]][1]
          }   
          #set unique name to dataset
          unique_names <- 
            make.unique( 
              unlist( 
                c(names(seurat_objects()), names(seurat_objects_), dataset_name) 
              ) 
            )
        
          dataset_name <- unique_names[length(unique_names)]
          
          Misc(object = seurat_object, slot = "name") <- dataset_name
          Misc(object = seurat_object, slot = "dir_name") <- multidataset_dirs[[i]]
          
          seurat_objects_ <- c(seurat_objects_, seurat_object)
          names(seurat_objects_)[[length(seurat_objects_)]] <- dataset_name
          
          #browser()
          print(seurat_object@misc$name)
          print(seurat_object)
          #print(seurat_objects())
          print("reading completed")
          
          ####
          
        }
        

        seurat_objects(c(seurat_objects(), seurat_objects_))
        
        print("seurat_objects")
        print(seurat_objects())
        
        #update mydatasets widget
        tmp_list <- c(input$mydatasets$left, as.list(names(seurat_objects_)))
        dataset_names_right(input$mydatasets$right)
        dataset_names(tmp_list)
        
        
        
      },
      error = function(e) {
        # return a safeError if a parsing error occurs
        #stop(safeError(e))
        print("error: the folder has no valid dataset")
        NULL
      }
    )
    
  })
  
  observeEvent(input$create_project, {
    browser()
    req(input$mydatasets)
    req(input$project_name)
    s_obj <- seurat_objects()[input$mydatasets$right]
    create_proj_gui(input$project_name, selected_seurat_objs = s_obj)
    print(last_created_project_id)
    #aa <- qc_server(last_created_project_id, input$tabs)
    qc_server(last_created_project_id, last_created_project_id)
    #projects[[project_id]] <- aa()
  })
  
  output$mydatasets_out <- renderUI({
    #req(input$dataset_dir)
    #req(dataset_names())
    req(input$load_dataset)
    chooserInput("mydatasets", "Available datasets", "Selected datasets",
                 dataset_names(), dataset_names_right(), size = 10, multiple = TRUE
    )
  })
  
  shinyFileSave(input, "save_project", roots = roots)
  shinyFileChoose(input, "load_project", roots = roots)
  
  observeEvent(input$tabs, {
    if (input$tabs != "Load Data") {
      shinyjs::enable("save_project")
    } else {
      shinyjs::disable("save_project")
    }
    
  })
  
  observeEvent(input$load_project,{
    browser()
    req(!is.integer(input$load_project))
    print(paste("Loading:", parseFilePaths(roots, input$load_project)$datapath))
    
    project_id <- load_project(parseFilePaths(roots, input$load_project)$datapath, projects)
    appendTab("tabs", project_tab(project_id), select = TRUE)
    
    prj <- list(
      projectname = projects[[project_id]][["project_name"]], 
      filename =  normalizePath(projects[[project_id]][["file_name"]]))
    
    renderd <- isolate(recent_projects_$renderd)
    renderd <- renderd[which(rownames(renderd)!=prj$filename),]
    renderd <- rbind(prj, renderd)
    rownames(renderd) <- renderd[["filename"]]
    recent_projects_$renderd <- renderd
    
    #if (file.exists(recent_projects_file_conf_())) {
      write.table(renderd, file = isolate(recent_projects_file_conf_()), col.names =  TRUE, row.names = FALSE, sep = ",", quote = FALSE)
    #}
    
    print(project_id)
    browser()
     #aa <- qc_server(project_id, project_id)
    qc_server(project_id, project_id)
    #projects[[project_id]] <- aa()
    
  })
  
  observeEvent(input$save_project,{
    req(input$tabs != "Load Data")
    req(!is.integer(input$save_project))
    #browser()
    print(paste("Saving:", input$tabs, "in", parseSavePath(roots, input$save_project)))
    project_id <- save_project(input$tabs, parseSavePath(roots, input$save_project))
    
    
    browser()
    prj <- list(
      projectname = projects[[project_id]][["project_name"]], 
      filename =  normalizePath(projects[[project_id]][["file_name"]]))
    
    recent_projects_$renderd[prj$filename,] <- prj
    
    #if (file.exists(recent_projects_file_conf_())) {
      write.table(isolate(recent_projects_$renderd), file = recent_projects_file_conf_(), col.names =  TRUE, row.names = FALSE, sep = ",", quote = FALSE)
    #}
    
  })
  
  
  #observeEvent(input$do_qc, {
  #  browser()
  #  quality_check(projects[[isolate(input$tabs)]][["seurat_objs"]], input$perc_zeros, input$sigma_mito, input$multiplet_rate)
  #})
  #
  output$selection <- renderPrint(
    input$mydatasets
  )
  
  
  
  
} #end server







shinyApp(ui, server)

### end of file -- app.R

