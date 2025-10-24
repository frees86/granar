#' @title Generate a root cross-section
#'
#' Functions to generate root cross section anatomy, based on global parameters, such as the mean size of cells and number of cell layers.
#' @param path The path to the input file
#' @param parameters the input parameter table. Can be obtain by read_param_xml()
#' @param verbatim TRUE = Generate text to follow the simulation process/ FALSE = no text
#' @param maturity_x TRUE = meta-xylem are labeled as stele cell /FALSE = no change in cell labeling (default)
#' @param paraview TRUE = cell wall data is set to make 3D object/ FALSE = cell wall data is not compatible for 3D object
#' @keywords root
#' @import xml2
#' @import purrr
#' @import dplyr
#' @import tidyverse
#' @import deldir
#' @import sf
#' @import packcircles
#' @export
#' @examples
#' # Load input
#' params <- read_param_xml(path = system.file("extdata", "root_monocot.xml", package = "granar"))
#' # Generate anatomy
#' root = create_anatomy(parameters = params)
#' # Visualize the simulation output
#' plot_anatomy(root)
#' # Write the simulation output
#' write_anatomy_xml(sim = root, path = system.file("extdata", "current_root.xml", package = "granar"))
#'

create_anatomy <- function(path = NULL,  # path to xml file
                           parameters = NULL,
                           verbatim = F,
                           maturity_x = F,
                           paraview = T){

  # CHECKING & SETTING PARAMETERS:
  ################################

  # Return NULL is no parameters are specified
  if( is.null(path) & is.null(parameters)){
    warning("Please specify a parameter set for the simulation")
    return(NULL)
  }
  if(!is.null(path)){
    params <- read_param_xml(path)
  }
  if(!(is.null(parameters))){
    params <- parameters
    if(nrow(params[params$name == "planttype",]) == 0){
      warning(paste0("Could not find the 'planttype' information in the parameter input"))
      return(NULL)
    }
    # Quality control
    to_find <- c("secondarygrowth", "planttype", "randomness", "xylem", "phloem", "stele", "endodermis", "exodermis", "epidermis", "aerenchyma", "pericycle", "cortex")
    for(tf in to_find){
      if (nrow(params[params$name == tf,]) == 0){
        warning(paste0("Could not find the '",tf,"' information in the parameter input"))
      }
    }
    cols_to_find <- c("name", "type", "value")
    for(ctf in cols_to_find){
      if (is.null(params[[ctf]])){
        warning(paste0("Could not find the '",ctf,"' column in the parameter input"))
        return(NULL)
      }
    }
  }

  # If the parameters do not contain information about the radius of the stele, we add it:
  if (length(params$value[params$name == "stele" & params$type == "radius"])==0) {
    if (length(params$value[params$name == "central_xylem" & params$type == "cell_diameter"])==1) {
      central_xylem_radius <- params$value[params$name == "central_xylem" & params$type == "cell_diameter"]/2
    } else {
      central_xylem_radius <- 0
    }
    stele_cell_radius <- params$value[params$name == "stele" & params$type == "cell_diameter"]
    stele_layers <- params$value[params$name == "stele" & params$type == "n_layers"]
    stele_radius <- central_xylem_radius + (stele_cell_radius*2 + 1) * stele_layers
    params <- rbind(params, data.frame(name = "stele",
                                      type = "radius",
                                      value = stele_radius))
  }

  # We set initial time:
  t_1 <- Sys.time()
  t1 <- proc.time()
  # We set the random factor:
  random_fact <- params$value[params$name == "randomness" & params$type == "intensity"] / 10 * params$value[params$name == "stele" & params$type == "cell_diameter"]
  # We set the random seed:
  seed <- params$value[params$name == "randomness" & params$type == "seed"]
  # We set the proportion of aerenchyma:
  proportion_aerenchyma <- params$value[params$name == "aerenchyma" & params$type == "proportion"]

  # CREATING THE CENTERS OF THE CELLS:
  ####################################

  # We create tables containing the information of each cell layer to be drawn:
  data_list <- cell_layer(params)
  layers <- data_list$layers # layers: cell_type, diameter, n_layers, order
  all_layers <- data_list$all_layers # expand layers, with e.g. the number of cells per layer and the radius

  # The x-coordinate and y-coordinate of the center of the cross section will be identical.
  # They are defined as the highest radius of the layers calculated above:
  center <- max(all_layers$radius)

  # We get the time used for creating layers:
  t2 <- proc.time()

  # We set the center of each cell:
  message("Creating the center of cells...")
  all_cells <- create_cells(all_layers, random_fact, random_seed=seed)
  # We get a summary of cells:
  summary_cells <- plyr::ddply(all_cells, plyr::.(type), summarise, n_cells = length(angle))
  # We relabel all cortex-related cells as "cortex":
  all_cells$type[grepl("cortex", all_cells$type)]<- "cortex"
  # We initialize the id_group variable:
  all_cells$id_group <- 0

  # Inclusion of the pith in the stele
  if(length(params$value[params$name == "pith" & params$type == "layer_diameter"]) > 0){
    all_cells <- make_pith(all_cells, params, center)
  }

  # # Addition of intercellular space and reshape cortex layers
  # # (Sub optimal process, may take a while)
  # if(length(params$value[params$name =="inter_cellular_space"]) > 0){
  #   all_cells <- rondy_cortex(params, all_cells, center)
  # }

  # We get the vascular system inside the stele:
  message("Creating vascular elements...")
  if(verbatim) message("Add vascular elements")
  # Case 1: No secondary growth
  if(length(params$value[params$name == "secondarygrowth"]) ==  0){
    all_cells <- vascular(all_cells, params, layers, center)
  } else if(params$value[params$name == "secondarygrowth"] == 0){
    all_cells <- vascular(all_cells, params, layers, center)
  # Case 2: Secondary growth
  } else if (params$value[params$name == "secondarygrowth"] == 1){
    # In case of secondary growth, we do circle packing:
    packing<-pack_xylem(all_cells, params, center)
    rm_stele <- all_cells%>%
      filter(type != "stele")
    new_cells <- rbind(packing, rm_stele)
    new_cells$id_cell <- 1:nrow(new_cells)
    all_cells <- new_cells
  }

  # Change parenchyma into stele, otherwise it will have serious problems later in MECHA
  all_cells$type[all_cells$type == "parenchyma"] = "stele"

  # PERFORMING THE TESSELATION TO GET THE CELL WALLS:
  ###################################################

  # # OPTION WITH WEIGHTED VORONOI (Apollonius package):
  # # We create a numeric matrix with only the x and y coordinates of the cell centers:
  # sites = NULL
  # for (i in seq(1,length(all_cells$x))) {
  #   sites <- rbind(sites, c(all_cells$x[i], all_cells$y[i]))
  # }
  # # We create a numeric vector containing only the radius of the corresponding cells:
  # radii <- all_cells %>%
  #   mutate(cell_radius = case_when(type=="cortex" ~ 1.0,
  #                                  .default = 1.0)) %>%
  #   pull(cell_radius)
  #
  # # We perform the Apollonius diagram calculation:
  # apo <- Apollonius(sites, radii)
  # # plotApolloniusGraph(apo, xlab = "x", ylab = "y")
  #
  # sites  <- apo[["diagram"]][["sites"]]
  # faces  <- apo[["diagram"]][["faces"]]
  # nsites <- nrow(sites)
  # radii  <- sites[, "weight"]
  # single_points <- apo[["graph"]][["sites"]]
  # edges  <- apo[["graph"]][["edges"]]
  # hsegments <- edges[["segments"]]
  # hrays     <- edges[["rays"]]
  #
  # sites <- as.data.frame(sites)
  # # sites <- all_cells %>% select(x,y)
  #
  # # View(apo)
  # # View(sites)
  # # View(hsegments)
  #
  # # Defining limits for the plot:
  # x <- extendrange(sites[, "x"]) # This extends the range by a small fraction
  # y <- extendrange(sites[, "y"]) # This extends the range by a small fraction
  # limits <- c(min(x[1L], y[1L]), max(x[2L], y[2L]))
  # # Plotting:
  # plot(NULL, xlim = limits, ylim = limits, asp = 1)
  # for(i in seq_along(hsegments)) {
  #   # print(hsegments[[i]])
  #   lines(hsegments[[i]], col="black", lwd = 2)
  #   points(sites$x[i], sites$y[i], pch=19, cex=0.2)
  # }
  #
  # # We create a new dataframe "segments12" that will contain four columns x1, y1, x2, y2:
  # segments12 = data.frame(matrix(nrow = length(hsegments), ncol = 5))
  # col_names= c("x1","y1","x2","y2","id_site")
  # colnames(segments12) = col_names
  # # We create a new dataframe "segments" that will contain columns x and y only:
  # segments = data.frame(matrix(nrow = length(hsegments)*2, ncol = 3))
  # col_names= c("x","y","id_site")
  # colnames(segments) = col_names
  #
  # # Initialization:
  # j <- 1
  # k <- 1
  # id_site <- 1
  # # For each line of hsegments:
  # for (i in seq(1,length(hsegments))) {
  #   # Implementing the four-column table:
  #   segments12$x1[i] = hsegments[[i]][1][1]
  #   segments12$y1[i] = hsegments[[i]][3][1]
  #   segments12$x2[i] = hsegments[[i]][2][1]
  #   segments12$y2[i] = hsegments[[i]][4][1]
  #   # segments12$id_site[i] = id_site
  #   # Implementing the two-column table:
  #   segments$x[j] = hsegments[[i]][1][1]
  #   segments$y[j] = hsegments[[i]][3][1]
  #   segments$x[j+1] = hsegments[[i]][2][1]
  #   segments$y[j+1] = hsegments[[i]][4][1]
  #   # segments$id_site[j] = id_site
  #   # segments$id_site[j+1] = id_site
  #   j <- j+2
  #   # k <- k+1
  #   # if (k>3) {
  #   #   k <-1
  #   #   id_site <- id_site + 1
  #   # }
  # }
  #
  # angle <- function(dir1_x1, dir1_y1, dir1_x2, dir1_y2,
  #                   dir2_x1, dir2_y1, dir2_x2, dir2_y2) {
  #
  #   print("Coucou!")
  #
  # }
  #
  #
  # # single_points = segments12 %>% mutate(x=x1, y=y1) %>% select(x,y) %>% slice(1)
  # # j <-2
  # # # For each line of hsegments:
  # # for (i in seq(1,length(hsegments))) {
  # #   starting_point <- segments12 %>% mutate(x=x1, y=y1) %>% select(x,y) %>% slice(i)
  # #   ending_point <- segments12 %>% mutate(x=x2, y=y2) %>% select(x,y) %>% slice(i)
  # #   starting_is_in_list <- paste(starting_point, collapse = ' ') %in% paste(single_points$x, single_points$y)
  # #   ending_is_in_list <- paste(ending_point, collapse = ' ') %in% paste(single_points$x, single_points$y)
  # #   if (!starting_is_in_list) {
  # #     print("Trying...")
  # #     print(starting_point$x[1])
  # #     single_points <- rbind(single_points, starting_point)
  # #     j <- j+1
  # #   }
  # #   if (!ending_is_in_list) {
  # #     print("Trying...")
  # #     print(ending_point$x[1])
  # #     single_points <- rbind(single_points, ending_point)
  # #     j <- j+1
  # #   }
  # # }
  # single_points <- as.data.frame(single_points)
  # single_points <- single_points %>% select(-3)
  # colnames(single_points) <- c("x", "y")
  # View(single_points)

  # PERFORMING TESSELATION WITH TRADITIONAL VORONOI DIAGRAM:
  ##########################################################

  # We make sure that the table do not contain any NA in x or y:
  initial_length = length(all_cells$x)
  NA_cells <- all_cells %>%
    filter(is.nan(x) | is.nan(y))
  all_cells <- all_cells %>%
    filter(!is.nan(x), !is.nan(y))
  filtered_length = length(all_cells$x)
  if (filtered_length < initial_length) {
    message("WATCH OUT: the table 'all_cells' before tesselation contained NA!")
    message("Here is the initial lines with NA:")
    print(NA_cells)
  }

  message("Performing tesselation...")
  vtess <- deldir(all_cells$x, all_cells$y, digits = 8)
  if(is.null(vtess)){return(NULL)}
  vorono_list <- cell_voro(all_cells, vtess, center)
  all_cells <- vorono_list$all_cells
  rs2 <- vorono_list$rs2

  cell_apo <- function(all_cells, vtess, center){

    # Get the size of the cells
    ids <- all_cells$id_cell
    all_cells$area <- NA
    all_cells$dist <- sqrt((all_cells$x - center)^2 + (all_cells$y - center)^2 )

    rs <- vtess$dirsgs[vtess$dirsgs$ind1 %in% ids |
                         vtess$dirsgs$ind2 %in% ids,]

    # Get the coordinates for every nodes in the voronoi
    rs <- rs %>% arrange(ind1)
    rs2 <- data.frame(x = rs$x1, y=rs$y1, id_cell = rs$ind1)
    rs2 <- rbind(rs2, data.frame(x = rs$x2, y=rs$y2, id_cell = rs$ind1))
    rs2 <- rbind(rs2, data.frame(x = rs$x2, y=rs$y2, id_cell = rs$ind2))
    rs2 <- rbind(rs2, data.frame(x = rs$x1, y=rs$y1, id_cell = rs$ind2))

    rs2 <- merge(rs2, all_cells[,c("id_cell", "type", "area", "dist", "angle", "radius", "id_layer", "id_group")], by="id_cell")

    rs2 <- rs2%>%filter(type != "outside")

    return(list(all_cells = all_cells, rs2 = rs2))

  }

  message("   > Tesselation has been done!")

  rs1 <- rs2 %>%
    dplyr::group_by(id_cell) %>%
    dplyr::mutate(my = mean(y),
                  mx = mean(x),
                  atan = atan2(y-my, x - mx)) %>%
    dplyr::arrange(id_cell, atan)%>%
    ungroup()

  rs1$id_point <- paste0(rs1$x,";",rs1$y)

  # Uniform cell by id_group
  if (max(rs1$id_group)>0) {
    if(verbatim) message("Smooth edge of large cells")
    rs1 <- smoothy_cells(rs1)
  }

  if(verbatim) message("Merging inter cellular space")
  rs1 <- fuzze_inter(rs1)

  ini_cortex_area <- sum(all_cells$area[all_cells$type %in% c( "cortex" ,"exodermis" , # ,"endodermis"
                                                               "epidermis", "inter_cellular_space")])

  if(proportion_aerenchyma > 0){
    rs1 = clear_nodes(rs1)
    if(verbatim) message("remove cells for aerenchyma")
    rs1 <- aerenchyma(params, rs1)
    # simplify septa
    if(verbatim) message("simplify septa between aerenchyma lacuna")
    rs1 <- septa(rs1)
  }else{
    cortex_area <- ini_cortex_area
  }

  # hairy epidermis # add-on 27-02-2020
  #-----------------------------------------

  if(length(params$value[params$name == "hair"] )!= 0){
    if(params$value[params$name == "hair" & params$type == "n_files"] > 0 ){
      if(verbatim) message("Add root hair")
      rs1 <- root_hair(rs1, params, center)
    }
  }
  if(verbatim) message("root hair, done")

  tt <- proc.time()
  # outputing the inputs
  output <- data.frame(io = "input", name = params$name, type = params$type,
                       value = params$value)

  if(length(which(is.na(rs1$x)))>0 & verbatim){
    print("NA in cell coordinate ... ")
  }
  rs1 = clear_nodes(rs1)

  # Reset the ids of the cells to be continuous
  ids <- data.frame(id_cell = unique(rs1$id_cell))
  ids$new <- c(1:nrow(ids))
  rs1 <- merge(rs1, ids, by="id_cell")
  rs1$id_cell <- rs1$new

  if(proportion_aerenchyma > 0){
    if(verbatim) message("create id_aerenchyma vector")
    id_aerenchyma <- unique(rs1$id_cell[rs1$aer == "aer"])
  }else{id_aerenchyma <- NA}

  all_cells <- merge(all_cells, ids, by="id_cell")
  all_cells$id_cell <- all_cells$new

  #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  # WATCH OUT: WE DO NOT DISTINGUISH ANY MORE METAXYLEM AND XYLEM
  # mX <- mean(rs1$area[rs1$type == "xylem"])
  # if(params$value[params$name == "planttype"] == 1){
  #   if(verbatim) message("for monocot, if xylem is above average, it is labeled as metaxylem")
  #   rs1$type[rs1$type == "xylem" & rs1$area > mX] <- "metaxylem"
  # }
  #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  if(length(params$value[params$name == "epidermis" & params$type == "remove"])>0){
    if(params$value[params$name == "epidermis" & params$type == "remove"]){
      all_cells = all_cells%>%filter(type != "epidermis")
      rs1 = rs1%>%filter(type != "epidermis")
    }
  }

  one_cells <- rs1%>%
    filter(!duplicated(id_cell)) #!duplicated(type)), !duplicated(id_group), !duplicated(area)

  all_cells <- merge(all_cells, one_cells, by = "id_cell")

  # adding the outputs by cell layers
  out <- plyr::ddply(all_cells, plyr::.(type.y), summarise, n_cells=length(type.y),
                     layer_area = sum(area.y),
                     cell_area = mean(area.y)) %>%
    mutate(name = type.y) %>%
    dplyr::select(-type.y) %>%
    tidyr::gather(key = "type", value = "value", n_cells, layer_area, cell_area) %>%
    mutate(io = "output")%>%
    dplyr::select(io, everything())
  output <- rbind(output, out)

  time <- as.numeric(Sys.time()-t_1)

  # finaly we add the outputs for the whole section
  output <-rbind(output, data.frame(io="output", name="all", type="n_cells", value = nrow(all_cells)))

  output <-rbind(output, data.frame(io="output", name="stelar", type="layer_area",
                                    value = sum(all_cells$area.y[all_cells$order < 4])))
  TCA <- sum(all_cells$area.y[all_cells$order > 3])
  output <-rbind(output, data.frame(io="output", name="cortex_alive_to_epidermis", type="layer_area",
                                    value = TCA))

  output <-rbind(output, data.frame(io="output", name="all", type="layer_area", value = sum(all_cells$area.y)))
  output <-rbind(output, data.frame(io="output", name="simulation", type="time", value = time))

  rs1$sorting <- c(1:nrow(rs1))

  nodes <- vertex(rs1)
  nodes <- nodes[!is.na(nodes$x), ]

  # In the MECHA python script, to have unmature metaxylem vessels
  # Metaxylem elements are turned into stele cell type
  if(maturity_x){
    tmp_m <- mean(nodes$area[nodes$type == "cortex"])
    nodes$type[nodes$type == "metaxylem" & nodes$area > tmp_m ] <- "stele"
  }

  # comment
  if(paraview){
    walls <- pv_ready(rs1)
    wall_length <- walls%>%select(-x, -y, -xx, -yy)%>% # ends_with(as.character(c(0:9)))
      select(starts_with("x"), starts_with("y"))%>%
      colnames()
    if(verbatim){
      print(wall_length)
    }

    wally <- walls[!duplicated(walls[,wall_length]),] %>%
      dplyr::select(wall_length)
    wally$id_wall <- c(1:nrow(wally))
    walls <- merge(walls, wally, by= wall_length)
    walls <- walls %>%
      # filter(!duplicated(id_wall))%>% # 30/09/2020
      arrange(sorting)

  }else{
    wally <- nodes[!duplicated(nodes[,c('x1', 'x2', 'y1', 'y2')]),] %>%
      dplyr::select(c(x1, x2, y1, y2))

    wally$id_wall <- c(1:nrow(wally))
    walls <- wally

    nodes <- merge(nodes, walls, by=c("x1", "x2", "y1", "y2"))
    nodes <- nodes %>%
      arrange(sorting)
  }

  id_aerenchyma <- unique(nodes$id_cell[nodes$type %in% c("aerenchyma", "inter_cellular_space")])
  id_aerenchyma <- id_aerenchyma-1

  print(Sys.time()-t_1)

  return(list(nodes = nodes,
              walls_nodes = walls,
              walls = wally,
              cells=all_cells,
              output = output,
              id_aerenchyma = id_aerenchyma))

}

`%!in%` <- compose(`!`, `%in%`)
