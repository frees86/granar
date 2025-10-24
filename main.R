# # To install the required external packages (if needed):
# install.packages("xml2")
# install.packages("htmltools")
# install.packages("retistruct")

# # To install weighted Voronoi diagrams:
# #######################################
# library(devtools)
# install_github("stla/gyro") # Note: if the installation fails, it may be necessary to install prior to gyro the packages cxhull and Polychrome
# install_github("stla/Apollonius")
#
# # Other possible source (required the latest versions of R and RStudio):
# library(devtools)
# install_github("andresgmejiar/lbmech")
# library(lbmech)

# # To load the required external packages:
library(xml2)
library(tidyverse)
library(plyr)
library(ggplot2)
library(dplyr)
library(deldir)
library(alphahull)
library(viridis)
library(Hmisc)
library(sf)
library(purrr)
library(packcircles)
library(retistruct)
library(smoothr)
library(geometry)
library(rgl)
library(stats)

# library(Apollonius)
library(splancs)

# To define the main path where all GRANAR files are stored locally:
main_path = "C:/Users/frees/granar/"
# To define the input path where all parameter files are stored:
input_path = paste0(main_path,"inputs/")

# # To load all internal R functions of the original package GRANAR:
# r_files <- list.files(paste0(main_path,"R/"), pattern = "\\.R$", full.names = TRUE)
# sapply(r_files, source)

# To load files from the original package GRANAR one by one:
source(paste0(main_path,'R/aer_in_geom_xml.R'))
source(paste0(main_path,'R/aerenchyma.R'))
source(paste0(main_path,'R/cell_layer.R'))
source(paste0(main_path,'R/cell_voro.R'))
source(paste0(main_path,'R/clear_nodes.R'))
source(paste0(main_path,'R/concavety.R'))
source(paste0(main_path,'R/create_anatomy.R'))
source(paste0(main_path,'R/create_cells.R'))
source(paste0(main_path,'R/fuzze_inter.R'))
source(paste0(main_path,'R/get_root_section.R'))
source(paste0(main_path,'R/granar_metadata.R'))
source(paste0(main_path,'R/layer_info.R'))
source(paste0(main_path,'R/make_it3d.R'))
source(paste0(main_path,'R/make_pith.R'))
source(paste0(main_path,'R/mecha_ApoSymp.R'))
source(paste0(main_path,'R/pack_xylem.R'))
source(paste0(main_path,'R/plot_3d_root.R'))
source(paste0(main_path,'R/plot_anatomy.R'))
source(paste0(main_path,'R/prep_geo.R'))
source(paste0(main_path,'R/pv_ready.R'))
source(paste0(main_path,'R/read_param_xml.R'))
source(paste0(main_path,'R/rondy_cortex.R'))
source(paste0(main_path,'R/root_hair.R'))
source(paste0(main_path,'R/septa.R'))
source(paste0(main_path,'R/smoothy_cells.R'))
source(paste0(main_path,'R/sum_area.R'))
source(paste0(main_path,'R/vascular.R'))
source(paste0(main_path,'R/vertex.R'))
source(paste0(main_path,'R/write_anatomy_obj.R'))
source(paste0(main_path,'R/write_anatomy_xml.R'))
source(paste0(main_path,'R/write_geo.R'))
source(paste0(main_path,'R/write_geo.R'))

# DEFINING NEW FUNCTIONS:
#########################

# FUNCTION FOR SHRINKING ONE CELL ACCORDING TO THE THICKENING OF THE CELL WALL:
shrinking_one_cell <- function(.data) {

  df <- .data
  if (length(df$x)==0 | length(df$y)==0) {
    message("PROBLEM: no cell can be shrinked!")
  }

  # We cover each coordinate of the wall of the current cell:
  for (i in seq(1,length(df$x))) {

    x = df$x[i]
    y=df$y[i]
    mx=df$mx[i]
    my=df$my[i]
    wall_thickness = df$wall_thickness[i]

    if (is.na(x) | is.na(y) | is.na(mx) | is.na(my)) {
      message("PROBLEM: this cell cannot be shrinked!")
      print(df)
    } else {
      if (x>mx & y>my) {
        angle = atan((y-my)/(x-mx))
        df$x[i] = x - cos(angle)*wall_thickness/2
        df$y[i] = df$y[i] - sin(angle)*wall_thickness/2
      } else if (x<mx & y>my) {
        angle = atan((y-my)/(mx-x))
        df$x[i] = x + cos(angle)*wall_thickness/2
        df$y[i] = df$y[i] - sin(angle)*wall_thickness/2
      } else if (x>mx & y<my) {
        angle = atan((my-y)/(x-mx))
        df$x[i] = x - cos(angle)*wall_thickness/2
        df$y[i] = df$y[i] + sin(angle)*wall_thickness/2
      } else if (x<mx & y<my) {
        angle = atan((my-y)/(mx-x))
        df$x[i] = x + cos(angle)*wall_thickness/2
        df$y[i] = df$y[i] + sin(angle)*wall_thickness/2
      }
    }
  }
  return(df)
}

# # Example of a simple dataframe where coordinates have to get closer to the middle:
# x <- c(0,1,1,0)
# y <- c(0,0,1,1)
# mx <- c(0.5,0.5,0.5,0.5)
# my <- c(0.5,0.5,0.5,0.5)
# m <- cbind(x,y,mx,my)
# df <- as.data.frame(m)
# df$x_new <- df$x
# df$y_new <- df$y
# df$angle <- df$x
# df$wall_thickness <- c(0.2, 0.4, 0.2, 0.2)
#
# df <- shrinking_one_cell(df)
# plot(df$x, df$y)
# plot(df$x_new, df$y_new)

# FUNCTION FOR SMOOTHING THE EDGES OF A SET OF CELLS:
####################################################

shrinking_all_cells <- function(cells) {

  # We apply the function 'smoothing_one_cell' to all cells in the simulation:
  new_cells <- cells %>%
    group_by(id_cell) %>%
    do(data.frame(shrinking_one_cell(.))) %>%
    drop_na()

  return(new_cells)
}

# # Example:
# sim1$nodes <- sim1$nodes %>%
#   mutate(wall_thickness=0.002)
#
# sim1$nodes <- shrinking_all_cells(sim1$nodes)
# plot_anatomy(sim1)

# FUNCTION FOR SMOOTHING THE EDGES OF A GIVEN CELL:
###################################################

smoothing_one_cell <- function(.data) {

  # # Example of a simple matrix to be smoothen:
  # x <- c(0,1,1,0)
  # y <- c(0,0,1,1)
  # m <- cbind(x,y)

  # We create a matrix with (x,y) coordinates of the different points of the matrix:
  m <- cbind(.data$x,.data$y)
  # As it is necessary to close the polygon, we add one last row corresponding to the first row:
  m <- rbind(m,m[1,])
  # We convert the matrix into a sf polygon:
  polygon <- st_polygon(list(m))

  # We get the smoothness factor to apply, if any:
  if (is.null(.data$smoothness)) {
    smoothness = 0
  } else {
    smoothness = .data$smoothness[1]
  }

  # SMOOTHING:
  #-----------

  # If the data is to be smoothened:
  if (as.numeric(smoothness) > 0) {

    # We smoothen the polygon (possible methods: 'chaikin', 'ksmooth', 'spline'):
    # smoothened_cell <- smoothr::smooth(polygon, method="chakin")
    smoothened_cell <- smoothr::smooth(polygon, method="ksmooth", smoothness=smoothness)
    # smoothened_cell <- smoothr::smooth(polygon, method="spline")

    # We create a new dataframe containing the (x,y) coordinates of the smoothened polygon:
    final_cell <- as.data.frame(st_coordinates(smoothened_cell)) %>%
      # We rename the first columns X and Y:
      dplyr::rename(x=X, y=Y) %>%
      # We remove unnecessary columns L1 and L2:
      select(-L1,-L2)

    } else {
    final_cell <- data.frame(cbind(.data$x,.data$y))
    colnames(final_cell)=c("x","y")
    }

  # RECREATING (x1,y1,x2,y2) COORDINATES AND WALL LENGTH:
  #------------------------------------------------------

  # At this stage, the table only contains x and y coordinates,
  # with the last line corresponding to a repetition of the first point to close the polygon.

  # We create a first table containing x1 and y1 coordinates:
  final_cell1 <- final_cell %>%
    # We define x1 and y1 as the first coordinates of an edge:
    dplyr::rename(x1=x, y1=y) %>%
    # We add an id column that will be used for merging this table with the second one:
    dplyr::mutate(id = row_number()) %>%
    # We remove the last line:
    dplyr::filter(dplyr::row_number() <= n()-1)

  # We create a second table containing only x2 and y2 coordinates:
  final_cell2 <- final_cell %>%
    # We remove the first line of the table:
    dplyr::slice(-1) %>%
    # We define x2 and y2 as the second coordinates of an edge:
    dplyr::rename(x2=x, y2=y) %>%
    # We add an id column that will be used for merging this table with the first one:
    dplyr::mutate(id = row_number())

  # We now combine the two tables so that each line has x1, y1, x2, y2, and x, y corresponding to the first coordinates:
  final_cell12 <- suppressWarnings(dplyr::left_join(final_cell1, final_cell2, by="id"))

  final_cell <- full_join(final_cell %>% dplyr::mutate(id = row_number()), final_cell12, by="id") %>%
  # final_cell <- reduce(list(final_cell1, final_cell2), left_join) %>%
    # # We remove the last line:
    # dplyr::filter(dplyr::row_number() <= n()-1) %>%
    # We remove the id column:
    select(-id)

  # Finally, we compute the distance "wall_length" according to x1, y1, x2, y2 coordinates:
  final_cell <- final_cell %>%
    dplyr::mutate(wall_length = ((x2-x1)^2 + (y2-y1)^2)^0.5) %>%
    dplyr::mutate(wall_length = case_when(is.na(wall_length) ~ 0,
                                          .default = wall_length)) %>%
    select(-c(x1,y1,x2,y2))

  return(final_cell)
}

# FUNCTION FOR SMOOTHING THE EDGES OF A SET OF CELLS:
#####################################################

smoothing_all_cells <- function(sim) {

  message("Smoothing the edges of the cells...")
  # We apply the function 'smoothing_one_cell' to all cells in the simulation:
  smoothened_cells <- sim$nodes %>%
    group_by(id_cell) %>%
    do(data.frame(smoothing_one_cell(.)))

    # The table smoothened_cells only contains id_cell and the x and y coordinates of each new points in the cell walls.
  # To combine this with the information from the original table, we extract the useful information from the original table:
  summary <- sim1$nodes %>%
    group_by(type, id_cell) %>%
    dplyr::summarize(radius=mean(radius),
                     mx=mean(mx),
                     my=mean(my))
  # And we combine this summary with the new table containing the x,y coordinates of the points of smoothened cells:
  new_cells <- right_join(summary, smoothened_cells)

  # # We make sure that the smoothened new cells do not contain additional types as the original simulated cells:
  # original_types <- sim$nodes %>% select(type, smoothness)
  # new_cells <- left_join(original_types, new_cells)

  return(new_cells)
}

# # Example of smoothing:
# smoothing=1
# if (smoothing > 0) {
#   # We create a table containing the coordinates of the edges of the new smoothened cells
#   new_cells <- smoothing_all_cells(sim1, smoothness = smoothing)
#   # We create a new simulation outputs by replacing the original nodes with the new nodes after smoothing:
#   new_sim1 <- sim1
#   new_sim1$nodes <- new_cells
#   # We plot the new simulation outputs:
#   plot_anatomy(new_sim1)
# }


# FUNCTION FOR CREATING A SUMMARY FROM THE SIMULATION:
#####################################################

computing_BRIDGES_outputs <- function(sim, params) {

  # Reference table for the summary:
  #---------------------------------

  # We select the few interesting variables from the simulation results and summarize them for each cell:
  cells_table <- sim$nodes %>%
    select(id_cell, type, wall_length, radius) %>%
    group_by(type, id_cell) %>%
    dplyr::summarize(d_tissue_in_mm_per_section = 2*mean(radius),
                     total_membrane_length_in_mm_per_cell = sum(wall_length))

  # We get the cell diameters from the parameters:
  filtered_diameters <- params %>%
    filter(type == "cell_diameter") %>%
    select(-type) %>%
    dplyr::rename(type = name, cell_diameter_in_mm = value)
  # We get the cell heights from the parameters:
  filtered_cell_heights <- params %>%
    filter(type == "cell_height") %>%
    select(-type) %>%
    dplyr::rename(type = name, cell_height_in_mm = value)
  # We get the cell wall thickness from the parameters:
  filtered_cell_wall <- params %>%
    filter(type == "wall_thickness") %>%
    select(-type) %>%
    dplyr::rename(type = name, cell_wall_thickness_in_mm = value)
  # We make one table with the 3 parameters:
  filtered_parameters <- full_join(filtered_diameters,filtered_cell_heights)
  filtered_parameters <- full_join(filtered_parameters, filtered_cell_wall)

  # We couple these parameters to the summary at the scale of cells:
  cells_table <- right_join(filtered_parameters, cells_table)

  # Computing the cross section of each actual cell:
  #-------------------------------------------------

  # We create a function that calculates the area of a polygon from the coordinates x and y of a dataframe:
  polygon_area = function(.data){
    m <- cbind(.data$x, .data$y)
    area = areapl(m)
    return(area)
  }

  # We now create a summary of the simulation results, calculating the cross section area of each cell:
  summary_cells <- sim$nodes %>%
    group_by(type, id_cell) %>%
    drop_na() %>%
    do(data.frame(cell_area=polygon_area(.)))

  # Computing total membrane surface and volume for each cell:
  #-----------------------------------------------------------

  # We bind the summary and the table together:
  summary_cells <- right_join(cells_table, summary_cells)

  # Computing the equivalent diameter of actual cells and their volume from their cross section:
  summary_cells <- summary_cells %>%
    mutate(equivalent_cell_diameter_in_mm = (cell_area / pi)^0.5*2,
           cell_volume_in_mm3_per_cell = cell_area * cell_height_in_mm)

  # Computing variables for each tissue:
  #-------------------------------------

  # We summarize the properties of cells for each tissue:
  summary_tissue <- summary_cells %>%
    group_by(type) %>%
    dplyr::summarize(
      # Number of cells per tissue:
      number_of_cells = n(),
      # Reference cell diameter from the parameters:
      cell_diameter_in_mm = mean(cell_diameter_in_mm),
      # Calculated mean equivalent cell diameter computed from the actual section of the cells:
      mean_equivalent_cell_diameter_in_mm = mean(equivalent_cell_diameter_in_mm),
      # Reference cell height from the parameters:
      cell_height_in_mm = mean(cell_height_in_mm),
      # Reference cell wall thickness from the parameters:
      cell_wall_thickness_in_mm = mean(cell_wall_thickness_in_mm),
      # Calculated internal and external radius of the tissue layer:
      internal_tissue_radius_in_mm = max(0, min(d_tissue_in_mm_per_section/2) - cell_diameter_in_mm/2),
      external_tissue_radius_in_mm = max(d_tissue_in_mm_per_section/2) + cell_diameter_in_mm/2,
      # Cumulated membrane length of the cells within the cross section:
      total_membrane_length_in_mm_per_section = sum(total_membrane_length_in_mm_per_cell),
      # Total area of the cells within the cross section:
      total_cells_area_in_mm2_per_section = sum(cell_area),
      # Total volume of the cells within one vertical cell layer:
      total_cells_volume_in_mm3_per_cell_layer = sum(cell_volume_in_mm3_per_cell))

  # We reorganize the tissue type in a specific order:
  summary_tissue <- summary_tissue %>%
    # We create a variable for sorting the type of cells in a specific order:
    mutate(type_order = case_when(type == "central_xylem" ~ 1,
                                  type == "xylem" ~ 2,
                                  type == "protoxylem" ~ 3,
                                  type == "phloem" ~ 4,
                                  type == "companion_cell" ~ 5,
                                  type == "stele" ~ 6,
                                  type == "pericycle" ~ 7,
                                  type == "endodermis" ~ 8,
                                  type == "cortex" ~ 9,
                                  type == "exodermis" ~ 10,
                                  type == "epidermis" ~ 11,
                                  .default = 100)) %>%
    # We reorganize the rows according to this specific order:
    arrange(type_order)

  # We compute the total cell membrane surface including the cross section at regular intervals (in mm2 per mm):
  summary_tissue <- summary_tissue %>%
    mutate(full_membranar_surface_area_in_mm2_per_mm
           = (total_membrane_length_in_mm_per_section * cell_height_in_mm) * 1 / (cell_height_in_mm + cell_wall_thickness_in_mm)
           + total_cells_area_in_mm2_per_section * (1 + 1 / (cell_height_in_mm + cell_wall_thickness_in_mm)))

  # We compute the total cell volume density, excluding the intercellular space (in mm2 per mm):
  summary_tissue <- summary_tissue %>%
    mutate(full_internal_cell_volume_in_mm3_per_mm = total_cells_volume_in_mm3_per_cell_layer * 1 / (cell_height_in_mm + cell_wall_thickness_in_mm))

  # We create another, special summary in order to count one edge of cell only one time
  # in case of neighbouring cells:
  table_restricted <- sim$nodes %>%
    # dplyr::select(id_cell, type, x, y, x1, x2, y1, y2, wall_length) %>% # For unknown reasons, this line lead to an error for smoothened_sim...
    # We now keep only the first line where a new combination of (x1,y1,x2,y2) appears:
    dplyr::distinct(x,y, .keep_all = TRUE)

  summary_restricted <- table_restricted %>%
    group_by(type) %>%
    dplyr::summarize(restricted_wall_length_in_mm_per_section = sum(wall_length))

  results <- full_join(summary_tissue, summary_restricted)

  return(results)
}


#######################################################################################################
# To run a simple anatomy simulation by GRANAR:
#######################################################################################################

# PERFORMING THE CLASSICAL SIMULATION WITHOUT SMOOTHING OR CELL SHRINKING:
#-------------------------------------------------------------------------

# We define the input path where the parameter file is stored:
input_path = paste0(main_path,"inputs/inputs.xml")
# We load the input parameters:
param1 <- read_param_xml(path = input_path)
# We create the corresponding root anatomy:
sim1  <- suppressWarnings(create_anatomy(parameters = param1))
# We plot the corresponding root section with default parameters:
suppressWarnings(plot_anatomy(sim1, xmin=0, ymin=0, xmax=0.8, ymax=0.8, xc=0.35, yc=0.35))

# PERFORMING A NEW SIMULATION WITH SMOOTHING AND CELL SHRINKING:
#---------------------------------------------------------------

# We get the wall_thickness and the required smoothness from the parameters:
filtered_wall_thickness_parameters <- param1 %>%
  filter(type == "wall_thickness") %>%
  select(-type) %>%
  dplyr::rename(type = name, wall_thickness = value)
# We couple both tables:
sim1$nodes <- left_join(sim1$nodes, filtered_wall_thickness_parameters)

# We get the required smoothness from the parameters:
filtered_smoothness_parameters <- param1 %>%
  filter(type == "smoothness") %>%
  select(-type) %>%
  dplyr::rename(type = name, smoothness = value)
# We couple both tables:
sim1$nodes <- left_join(sim1$nodes, filtered_smoothness_parameters)

# We make sure that values of wall_thickness and smoothness are also available for the protoxylem:
protoxylem_wall_thickness <- filtered_wall_thickness_parameters %>%
  filter(type=="pericycle")
protoxylem_smoothness <- filtered_smoothness_parameters %>%
  filter(type=="pericycle")
sim1$nodes <- sim1$nodes %>%
  mutate(wall_thickness = case_when(type=="protoxylem" ~ protoxylem_wall_thickness[1,2],
                                    .default=wall_thickness),
         smoothness = case_when(type=="protoxylem" ~ protoxylem_smoothness[1,2],
                                .default=smoothness))

# We now increase the thickness of the cell walls and reduce the size of the cells:
shrinked_sim1 <- sim1
shrinked_sim1$nodes <- shrinking_all_cells(sim1$nodes %>% drop_na())
plot_anatomy(shrinked_sim1, plotting_cell_centers = TRUE,
             xmin=0, ymin=0, xmax=0.8, ymax=0.8, xc=0.35, yc=0.35)

# We then smooth the new cells:
new_cells <- smoothing_all_cells(shrinked_sim1)
# We create a new simulation output by replacing the original nodes with the new nodes after smoothing:
smoothened_sim1 <- sim1
smoothened_sim1$nodes <- new_cells

# We plot the new simulation outputs:
plot_anatomy(smoothened_sim1, xmin=0, ymin=0, xmax=0.8, ymax=0.8, xc=0.35, yc=0.35)
# # Or we plot the section without specific cell layers:
# plot_anatomy(smoothened_sim1,
#              xmin=0, ymin=0, xmax=0.8, ymax=0.8,
#              hidden_cell_layers=c("exodermis", "epidermis"))

# We record the plot:
tiff(filename=paste0(main_path,"outputs/Plots/plot.tiff"), height = 12, width = 16, units = 'cm', compression = "lzw", res = 300)
plot_anatomy(smoothened_sim1, xmin=0, ymin=0, xmax=0.8, ymax=0.8, xc=0.35, yc=0.35)
# # Or we plot the section without specific cell layers:
# plot_anatomy(smoothened_sim1,
#              xmin=0, ymin=0, xmax=0.8, ymax=0.8, xc=0.35, yc=0.35,
#              hidden_cell_layers=c("exodermis", "epidermis"))
dev.off()

# COMPUTING INTERESTING PROPERTIES:
#----------------------------------

# We compute results for the original section BEFORE shrinking and smoothing:
original_results <- computing_BRIDGES_outputs(sim=sim1, params=param1)
# We compute results for the original section AFTER shrinking and smoothing:
smoothened_results <- computing_BRIDGES_outputs(sim=smoothened_sim1, params=param1)
# We extract the few interesting columns from the original section:
selected_original_results <- original_results %>%
  dplyr::select(type,
                internal_tissue_radius_in_mm,
                external_tissue_radius_in_mm,
                total_cells_area_in_mm2_per_section,
                total_cells_volume_in_mm3_per_cell_layer,
                restricted_wall_length_in_mm_per_section) %>%
  dplyr::rename(original_total_cells_area_in_mm2_per_section = total_cells_area_in_mm2_per_section,
                original_total_cells_volume_in_mm3_per_cell_layer = total_cells_volume_in_mm3_per_cell_layer,
                wall_length_in_mm_per_section = restricted_wall_length_in_mm_per_section)
# We remove the few useless columns from the section after shrinking and smoothing:
selected_smoothened_results <- smoothened_results %>%
  dplyr::select(-c(internal_tissue_radius_in_mm,
                   external_tissue_radius_in_mm,
                   restricted_wall_length_in_mm_per_section))
# We combine both tables and compute additional variables (i.e. cell wall area and cell wall volume):
final_results <- full_join(selected_smoothened_results, selected_original_results) %>%
  dplyr::mutate(
    # The cell wall area per section is the difference of area between original and shrinked & smoothened cells:
    total_cell_wall_area_in_mm2_per_section
    = original_total_cells_area_in_mm2_per_section - total_cells_area_in_mm2_per_section,
    total_cell_wall_volume_in_mm3_per_mm
    = (original_total_cells_volume_in_mm3_per_cell_layer
       - total_cells_volume_in_mm3_per_cell_layer)
    * 1 / (cell_height_in_mm + cell_wall_thickness_in_mm)
    + original_total_cells_area_in_mm2_per_section * cell_wall_thickness_in_mm
    * (1 + 1 / (cell_height_in_mm + cell_wall_thickness_in_mm))
  ) %>%
  # We remove the restricted wall length as it does not seem to work properly:
  select(-c(original_total_cells_volume_in_mm3_per_cell_layer,
            wall_length_in_mm_per_section,
            type_order))
View(final_results)
write.csv(final_results, paste0(main_path,"outputs/GRANAR_results.csv"), row.names=FALSE, quote=FALSE)


# # We can also plot the section showing apoplastic barriers
# # [apo_bar: displays apolastic barrier when col = "segment".
# # 1 endodermal casparian strip,
# # 2 fully suberized endodermis,
# # 3 fully suberized endodermis and an exodermal casparian strip,
# # 4 exodermis and endodermis are fully suberized]
# plot_anatomy(sim1, col="segment", apo_bar=2)
# dev.copy(png,filename=paste0(main_path,"outputs/Plots/plot1.png"))
# dev.off ()

# # We can save the results in xml format:
# write_anatomy_xml(sim = sim1, path = paste0(main_path,"outputs/Anatomies/Wheat.xml"))

# # We can reload a previously-stored root anatomy (plotting it does not work, though!):
# sim2 <- get_root_section(path = paste0(main_path,"outputs/Anatomies/Wheat.xml"))
