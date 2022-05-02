### project_mgr.R

cosdeg_version <- "v0.8.0"

projects <- list()


###

normalize_proj_name <- function(project_name) {
  proj_names <- names(projects)
  tmp_project_name_ <- stringi::stri_replace_all_charclass(project_name, '[,.;:.|!?$\\ /\\|*%\\&()@#]', '_', merge=TRUE)
  tmp_project_name <- tmp_project_name_
  for (counter in  1:10) {
    if (any(stringi::stri_detect_fixed(proj_names, pattern = tmp_project_name))) {
      tmp_project_name <- paste0(tmp_project_name_, "_", counter)
    }
    else
      return(tmp_project_name)
  }
  if (counter == 10)
    return(NULL)
}


normalize_tab_name <- function(project_name) {
  tab_names <- c()
  for (proj in projects) {
    tab_names <- c(tab_names, proj[["tab_name"]])
  }
  
  tmp_project_name_ <- project_name
  tmp_project_name <- tmp_project_name_
  for (counter in  1:10) {
    if (any(stringi::stri_detect_fixed(tab_names, pattern = tmp_project_name))) {
      tmp_project_name <- paste0(tmp_project_name_, " (", counter,")")
    }
    else
      return(tmp_project_name)
  }
  if (counter == 10)
    return(NULL)
}



normalize_tab_id <- function(project_name) {
  tab_ids <- c()
  for (proj in projects) {
    tab_ids <- c(tab_ids, proj[["tab_id"]])
  }
  
  tmp_project_name_ <- stringi::stri_replace_all_charclass(project_name, '[,.;:.|!?$\\ /\\|*%\\&()@#]', '_', merge=TRUE)
  tmp_project_name <- tmp_project_name_
  
  for (counter in  1:10) {
    if (any(stringi::stri_detect_fixed(tab_ids, pattern = tmp_project_name))) {
      tmp_project_name <- paste0(tmp_project_name_, "_", counter)
    }
    else
      return(tmp_project_name)
  }
  if (counter == 10)
    return(NULL)
}



create_proj <- function(project_name, selected_seurat_objs = NULL ) {
  #browser()
  tmp_ <- normalize_proj_name(project_name)
  if (!is.null(tmp_))
    project_id <- tmp_
  else {
    # rise warning
  }
  
  project <- list()
  
  project[["COSDEG"]] <- paste("COSDEG", cosdeg_version)
  project[["project_name"]] <- project_name
  
  project[["file_name"]] <- NULL
  
  tmp_ <- normalize_tab_name(project_name)
  if (!is.null(tmp_))
    project[["tab_name"]] <- tmp_
  else {
    # rise warning
  }
  
  
  tmp_ <- normalize_tab_id(project_id)
  if (!is.null(tmp_))
    project[["tab_id"]] <- tmp_
  else {
    # rise warning
  }
  
  project[["modified"]] <- TRUE
  
  project[["seurat_objs"]] <- selected_seurat_objs  #copy
  
  projects[[project_id]] <<- project
  
  return(project_id)
}



save_project <- function(project_id, file_name) {
  browser()
  file_name <- file_name$datapath
  
  projects[[project_id]][["file_name"]] <- file_name
  projects[[project_id]][["project_id"]] <- project_id
  
  
  
  project <- projects[[project_id]]
  COSDEG <- project[["COSDEG"]]
  
  save(project, COSDEG , precheck = TRUE, file = file_name)
  return(project_id)
}


loadRData <- function(filename, varname=NULL){
  #loads an RData file, and returns it
  load(filename)
  if (is.null(varname))
    return(ls()[(ls() != "filename") & (ls() != "varname") ])
  else
    return(get(ls()[ls() == varname]))
}


load_project <- function(file_name, projects) {
  browser()
  #file_name <- file_name$datapath
  load(file_name)
  #browser()
  if (!any(ls() == "COSDEG"))
    return(FALSE)
  
  if (!any(ls() == "project"))
    return(FALSE)
  
  tmp_ <- normalize_proj_name(project[["project_name"]])
  if (!is.null(tmp_))
    project_id <- tmp_
  else {
    # rise warning
  }
  
  #project_ <- list()
  
  #project[["COSDEG"]] <- COSDEG
  #project[["project_name"]] <- project_name
  
  project[["file_name"]] <- file_name
  
  tmp_ <- normalize_tab_name(project[["project_name"]])
  if (!is.null(tmp_))
    project[["tab_name"]] <- tmp_
  else {
    # rise warning
  }
  
  
  tmp_ <- normalize_tab_id(project[["project_id"]])
  if (!is.null(tmp_))
    project[["tab_id"]] <- tmp_
  else {
    # rise warning
  }
  
  project[["modified"]] <- FALSE
  
  #project_[["seurat_objs"]] <- seurat_objs  #copy
  
  projects[[project_id]] <<- project
  
  return(project_id)
}


######





#project_name <- "dati 1"
#create_proj(project_name = project_name)
#create_proj(project_name = project_name)

#project_id <- names(projects)[[2]]
#save_project(project_id = project_id, file_name = list(datapath="./test.Rdata"))
#load_project(file_name = list(datapath="./test.Rdata"), projects = projects)


### end of file -- project_mgr.R