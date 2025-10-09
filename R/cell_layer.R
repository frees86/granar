#' @title Place one line of cell center
#'
#' Create a list
#' @param params The input dataframe
#' @keywords root layer
#' @export
#' @examples
#' # data_list <- cell_layer(params)
#'

cell_layer <- function(params){

  # We create a dataframe containing the main properties of layer groups, i.e. the cell diameter,
  # the number of layers per group, and the order by which the groups are ranked:
  layers <- params %>%
    filter(type %in% c("cell_diameter","n_layers","order")) %>%
    tidyr::spread(type, value) %>%
    # We remove any layer if n_layers is nil:
    filter(!is.na(n_layers) & n_layers > 0) %>%
    # We sort the different layers according to 'order':
    arrange(order)
  # NOTE: Here, vascular tissues have been removed, unless n_layers was >0!

  # We get the cell diameter of the outer layers (it is usually the epidermis):
  outer_layers <- layers %>%
    filter(order==max(order))
  outer_layer_diameter = outer_layers$cell_diameter
  # We create two additional "outside" layers to serve as boundary for the tesselation algorithm:
  layers <- rbind(layers, data.frame(name="outside",
                                     n_layers=2,
                                     cell_diameter=outer_layer_diameter,
                                     order = max(layers$order)+1))

  # We get the number of layers of stele cells (excluding pericycle and endodermis):
  if (length(params$value[params$name == "stele" & params$type == "n_layers"])==1) {
    layers$n_layers[layers$name == "stele"] <- params$value[params$name == "stele" & params$type == "n_layers"]
  } else {
    stele_radius <- params$value[params$name == "stele" & params$type == "radius"]
    layers$n_layers[layers$name == "stele"] <- round(stele_radius / layers$cell_diameter[layers$name == "stele"])
  }

  # We get one row per actual cell layer:
  all_layers <- NULL
  for(i in c(1:nrow(layers))){
    for(j in c(1:layers$n_layers[i])){
      all_layers <- rbind(all_layers, layers[i,])
    }
  }

  # We now calculate the number of cells per layer, the radius and perimeter of the circle
  # corresponding to each layer, and the angle between each cell within one layer:
  all_layers <- layer_info(all_layers)

  # We return a list containing the detailed or general information about the layers:
  return(list(all_layers = all_layers , layers = layers))
}
