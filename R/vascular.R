#' @title Add the vascular element in the stele
#'
#' remove stele cells and place vascular elements
#' add round shaped boundaries for xylem elements
#' @param params The input dataframe
#' @param layers the layer dataframe
#' @param center The cross-section center
#' @keywords root
#' @export
#' @examples
#' # all_cells <- vascular(all_cells, params, layers, center)
#'

vascular <- function(all_cells, params, layers, center){

  # READING GENERAL PARAMETERS:
  #############################

  # We get the number of central xylem files from the parameters:
  if(length(params$value[params$name == "central_xylem" & params$type == "cell_diameter"]) > 0){
    n_central_xylem_files <- params$value[params$name == "central_xylem" & params$type == "n_files"]
  } else {
    n_central_xylem_files = 0
  }
  # We get the number of metaxylem files from the parameters:
  if (length(params$value[params$name == "xylem" & params$type == "n_files"]) == 1) {
    n_xylem_files = params$value[params$name == "xylem" & params$type == "n_files"]
  } else {
    n_xylem_files = 0
  }

  # We get the number of protoxylem cells to create for each metaxylem cell:
  if (length(params$value[params$name == "protoxylem" & params$type == "n_cells_per_file"]) == 1) {
    n_protoxylem_cells = params$value[params$name == "protoxylem" & params$type == "n_cells_per_file"]
  } else {
    n_protoxylem_cells <- 0
  }

  # We get the plant type from the parameters:
  plant_type <- params$value[params$name == "planttype"]

  # We define k_max_cortex as the maximal id_cell from cortex cells:
  if(length(all_cells$id_group[all_cells$type == "cortex"])> 0) {
    k_max_cortex <- max(all_cells$id_group[all_cells$type == "cortex"])
  } else {
    k_max_cortex  <- 0
  }
  k_max <- max(all_cells$id_group)

  # CASE 1: MONOCOT PLANT
  #######################

  if(plant_type == 1){

    # We initialize the id_group:
    k <- 1

    # DEFINING THE CIRCULAR FRONTIERS OF CENTRAL XYLEM:
    #--------------------------------------------------

    # In case the central xylem has not been considered as one layer:
    if (length(params$value[params$name == "central_xylem"]) > 0) {

      # We get the expected diameter of the central xylem cell:
      central_xyl_d <- params$value[params$name == "central_xylem" & params$type=="cell_diameter"]
      stele_d <- params$value[params$name == "stele" & params$type=="cell_diameter"]

      # x_cir <- seq(-0.95,0.95,0.95/4)
      x_cir <- seq(-0.95,0.95,0.95/10)
      y_p <- sqrt(1-x_cir^2)
      y_m <- -sqrt(1-x_cir^2)
      cir <- tibble(x = rep(x_cir,2), y = c(y_p,y_m))
      cir <- cir*abs((central_xyl_d-stele_d)*0.95)/2
      cir$x <- cir$x + center
      cir$y <- cir$y + center
      xyl_central_frontier <- data.frame(angle = 0,
                                         radius = 0,
                                         x = cir$x,
                                         y = cir$y,
                                         id_layer = 201,
                                         id_cell = 1,
                                         type = "central_xylem",
                                         order = 1.1,
                                         id_group = k)
      k <- k+1

      # ADDING THE FRONTIERS OF THE CENTRAL XYLEM:
      all_cells <- all_cells[all_cells$type != "central_xylem",]
      all_cells <- rbind(all_cells, xyl_central_frontier)
      # all_cells$id_group[all_cells$type == "central_xylem"] <- all_cells$id_group[all_cells$type == "central_xylem" & all_cells$id_group != 0] + k_max_cortex
      all_cells$id_group[all_cells$type == "central_xylem"] <- all_cells$id_group[all_cells$type == "central_xylem" & all_cells$id_group != 0] + k_max
    }

    # Finally, we remove the stele cells inside the perimeter of the central xylem:
    all_cells <- all_cells %>%
      filter(!((x-center)^2 + (y - center)^2 < (central_xyl_d/2-stele_d/2)^2
               & type == "stele"))

    # DEFINING THE POSITION OF METAXYLEM VESSELS:
    #--------------------------------------------

    # We get the radius of stele cells from the parameters:
    stele_d <- params$value[params$name == "stele" & params$type=="cell_diameter"]

    # We get the radius of xylem cells from the parameters:
    if (length(params$value[params$name == "xylem" & params$type == "cell_diameter"]) == 1) {
      radius_xylem_cell = params$value[params$name == "xylem" & params$type == "cell_diameter"]/2
    } else {
      radius_xylem_cell = 0.012
    }
    # If there is only one xylem vessel, it will be located in the middle of the root section:
    if(n_xylem_files == 1){
      xyl <- data.frame(r = 0,
                        d = radius_xylem_cell*2)
    } else if (n_xylem_files > 0) {
      xyl <- data.frame(r = params$value[params$name == "xylem" & params$type == "radius"],
                        d = radius_xylem_cell*2)
    }

    # We initialize the dataframe containing (x,y) coordinates of metaxylem vessels:
    all_xylem <- NULL
    # We define the range of angles corresponding to the orientation of each metaxylem:
    if (n_xylem_files > 0) {
      angle_seq <- seq(from = 0, to = (2*pi)-(2 * pi) / n_xylem_files, by = (2 * pi) / n_xylem_files)

      # For each metaxylem vessel:
      for(angle in angle_seq){
        # We define the (x,y) coordinate of the center of the current metaxylem according to the angle:
        x <- center + (xyl$r[1] * cos(angle))
        y <- center + (xyl$r[1] * sin(angle))
        # We add this new metaxylem vessel to the general dataframe of metaxylem vessels:
        all_xylem <- rbind(all_xylem, data.frame(x = x,
                                                 y = y,
                                                 d = xyl$d[1],
                                                 angle = angle,
                                                 id_group = k))
        # We add this new metaxylem vessel to the general dataframe of all cells:
        all_cells <- rbind(all_cells, data.frame(
          angle = angle,
          radius = xyl$r[1],
          x = x,
          y = y,
          id_layer = 20,
          id_cell = 1,
          order = 1.5,
          type = "xylem",
          id_group = k))
        # We increment the id_group for the next vessel:
        k <- k+1
      }

      # Finally, we remove the stele cells inside the perimeter of metaxylem vessels:
      if (!is.null(all_xylem)) {
        for(i in c(1:nrow(all_xylem))){
          all_cells <- all_cells %>%
            filter(!((x-all_xylem$x[i])^2 + (y - all_xylem$y[i])^2 < (all_xylem$d[i]/2-stele_d/2)^2
                     & type == "stele"))
        }
      }
    }

    # DEFINING THE CIRCULAR FRONTIERS OF METAXYLEM VESSELS:
    #------------------------------------------------------

    if (n_xylem_files > 0) {
      x_cir <- seq(-0.95,0.95,0.95/10)
      y_p <- sqrt(1-x_cir^2)
      y_m <- -sqrt(1-x_cir^2)
      xyl_frontier <- NULL

      for (i_xyl in 1:nrow(all_xylem)) {
        tmp <- all_xylem[i_xyl,]
        cir <- tibble(x = rep(x_cir,2), y = c(y_p,y_m))
        cir <- cir*abs(tmp$d - stele_d)/2
        cir$x <- cir$x + tmp$x
        cir$y <- cir$y + tmp$y
        xyl_frontier <- rbind(xyl_frontier, data.frame(angle = tmp$angle,
                                                       radius = xyl$r[1],
                                                       x = cir$x,
                                                       y = cir$y,
                                                       id_layer = 20,
                                                       id_cell = 1,
                                                       type = "xylem",
                                                       order = 1.5,
                                                       id_group = k))
        # We increment the id_group:
        k <- k+1
      }

      # We combine the metaxylem frontiers with the other, non-xylem cells:
      all_cells <- all_cells[all_cells$type != "xylem",]
      all_cells <- rbind(all_cells, xyl_frontier)
      # We reset the id_group of metaxylem vessels:
      # all_cells$id_group[all_cells$type == "xylem"] <- all_cells$id_group[all_cells$type == "xylem" & all_cells$id_group != 0] + k_max_cortex
      all_cells$id_group[all_cells$type == "xylem"] <- all_cells$id_group[all_cells$type == "xylem" & all_cells$id_group != 0] + k_max
      }

    # CREATION OF PROTOXYLEM VESSELS:
    #--------------------------------
    if (n_protoxylem_cells> 0) {

      # We get the radius of metaxylem cells from the parameters:
      if (length(params$value[params$name == "xylem" & params$type == "cell_diameter"]) == 1) {
        radius_xylem_cell = params$value[params$name == "xylem" & params$type == "cell_diameter"]/2
      } else {
        radius_xylem_cell = 0.012
      }
      # We get the general radius of the circle around which metaxylem cells are located:
      xylem_radius = params$value[params$name == "xylem" & params$type == "radius"]
      # We get all angles corresponding to the position of the metaxylem files:
      angle_seq_proto <- seq(from = 0, to = (2*pi)-(2 * pi) / (n_xylem_files), by = (2 * pi) / (n_xylem_files))
      # For each metaxylem file:
      for(angle in angle_seq_proto){
        # We get the coordinate of the metaxylem vessel:
        x1 <- center + (xylem_radius * cos(angle))
        y1 <- center + (xylem_radius * sin(angle))
        # For each protoxylem cell to create:
        for (step in seq(1,n_protoxylem_cells)) {
          # We find the closest pericycle cell and assign it to a protoxylem cell:
          all_cells <- all_cells %>%
            mutate(type = as.character(type)) %>%
            mutate(dist_protoxyl = sqrt((x-x1)^2 + (y-y1)^2)) %>%
            mutate(dist_protoxyl = ifelse(type == "pericycle"
                                          & sqrt((x-center)^2 + (y-center)^2) > sqrt((x1-center)^2 + (y1-center)^2) + 0.5*radius_xylem_cell,
                                          dist_protoxyl,
                                          100)) %>%
            mutate(type = ifelse(dist_protoxyl == min(dist_protoxyl), "protoxylem", type))
        }

      }

      # Eventually, we remove the intermediary variable dist_protoxyl:
      if("dist_protoxyl" %in% colnames(all_cells)) {
        all_cells <- all_cells %>% select(-dist_protoxyl)
      }
    }

    # PHLOEM CREATION & POSITIONNING:
    #--------------------------------

    # We get the number of phloem sieve tubes from the parameters:
    if (length(params$value[params$name == "phloem" & params$type == "n_files"]) == 1) {
      n_phloem_files = params$value[params$name == "phloem" & params$type == "n_files"]
    } else {
      n_phloem_files = n_xylem_files
    }

    # We get the radius of phloem cells from the parameters:
    if (length(params$value[params$name == "phloem" & params$type == "cell_diameter"]) == 1) {
      radius_phloem_cell = params$value[params$name == "phloem" & params$type == "cell_diameter"]/2
    } else {
      radius_phloem_cell = 0.012
    }
    # We get the radius of companion cells from the parameters:
    if (length(params$value[params$name == "companion_cell" & params$type == "cell_diameter"]) == 1) {
      radius_companion_cell = params$value[params$name == "companion_cell" & params$type == "cell_diameter"]/2
    } else {
      radius_companion_cell = 0.006
    }
    # We get the number of companion cells to create for each phloem cell:
    if (length(params$value[params$name == "companion_cell" & params$type == "n_cells_per_file"]) == 1) {
      n_companion_cells = params$value[params$name == "companion_cell" & params$type == "n_cells_per_file"]
    } else {
      n_companion_cells = 0
    }

    if (n_phloem_files > 0) {

      # We define the positions of the center of a phloem cell within the stele:
      phl <- data.frame(r = params$value[params$name == "phloem" & params$type == "radius"],
                        d = params$value[params$name == "phloem" & params$type == "cell_diameter"])

      # We will cover the position of each phloem cell:
      angle_seq_ph <- seq(from = ((2 * pi) / n_phloem_files) /2,
                          to = (2 * pi),
                          by = (2 * pi) / n_phloem_files)

      # For each new phloem cell:
      for(angle in angle_seq_ph){

        x1 <- center + (phl$r[1] * cos(angle))
        y1 <- center + (phl$r[1] * sin(angle))

        # PHLOEM TRANSFORMATION:
        # We find the stele cell closest to the target position and assign it to a phloem vessel:
        all_cells <- all_cells %>%
          mutate(type = as.character(type)) %>%
          # We compute the distance of each cell from the theoretical phloem cell center:
          mutate(dist_phl = sqrt((x-x1)^2 + (y-y1)^2)) %>%
          # We attribute a ridiculous value of 100 mm to dist_phl for any cell that is not a stele cell:
          mutate(dist_phl = ifelse(type == "stele",
                                   dist_phl, # If True, we keep the original value of dist_phl
                                   100)) %>% # If False, we replace dist_phl by 100
          # We replace the stele type if the distance to the phloem circle is the lowest one:
          mutate(type = ifelse(dist_phl == min(dist_phl),
                               "phloem", # If True, we set the cell type as "phloem"
                               type)) # If False, we keep the original cell type

        # We get the actual new coordinate of the phloem cell that has just been created:
        phloem_cell = all_cells %>%
          filter(type=="phloem" & dist_phl==min(dist_phl))
        x_ph = phloem_cell$x[1]
        y_ph = phloem_cell$y[1]

        # COMPANION CELL TRANSFORMATION:
        # We replace the nearest stele cell by a companion cell,
        # and we repeat this for as many time as necessary for the current phloem cell:
        for (n in seq(1,n_companion_cells)){
          all_cells <- all_cells %>%
            mutate(type = as.character(type)) %>%
            # We now attribute a ridiculous value of 100 mm to dist_phl for any non-stele cell,
            # for any stele cell located closer to the root center than the phloem cell center:
            mutate(dist_phl = ifelse(type == "stele"
                                     & sqrt((x-center)^2 + (y-center)^2) > sqrt((x_ph-center)^2 + (y_ph-center)^2) + 0.5*radius_phloem_cell,
                                     dist_phl, # If True, we keep the original value of dist_phl
                                     100)) %>% # If False, we replace dist_phl by 100
            # We replace the stele type if the distance to the phloem circle is the lowest one:
            mutate(type = ifelse(dist_phl == min(dist_phl) & min(dist_phl) != 100,
                                 "companion_cell", # If True, we set the cell type as "companion_cell"
                                 type)) # If False, we keep the original cell type
            # # We move the (x,y) coordinates of the companion cell to have it closer to the phloem cell:
            # mutate(x = ifelse(type == "companion_cell",
            #                   x + radius*sin(angle)*cos(angle) - radius_companion_cell ,
            #                   x)) %>%
            # mutate(y = ifelse(type == "companion_cell",
            #                   y - radius*(sin(angle))^2 - radius_companion_cell,
            #                   y))
        }

        # We select the companion cells:
        companion_cells <- all_cells %>%
          select(type, id_cell, x, y, radius) %>%
          filter(type=="companion_cell")

        if (nrow(companion_cells)>0) {
          # We get the direction coefficient for the axis formed by the section center and the phloem cell:
          m_ph = (y_ph - center) / (x_ph - center)
          for (i in seq(1,nrow(companion_cells))) {
            m_cc = (companion_cells$y[i] - y_ph)/(companion_cells$x[i] - x_ph)
            alpha = atan(abs((m_cc-m_ph)/(1+m_cc*m_ph)))
            x = companion_cells$x[i] + companion_cells$radius[i]*sin(alpha)*cos(alpha)
            y = companion_cells$y[i] - companion_cells$radius[i]*(sin(alpha))^2
            companion_cells$x_new[i] = x
            companion_cells$y_new[i] = y
          }
        }
      }
      # all_cells <- full_join(all_cells, companion_cells)
      # View(all_cells)
      #
      # all_cells <- all_cells %>%
      #   mutate(x = ifelse(type=="companion_cell",
      #                     x_new,
      #                     x),
      #          y = ifelse(type=="companion_cell",
      #                     y_new,
      #                     y))
      # # Eventually, we remove the intermediary variables:
      # all_cells <- all_cells %>% select(-c(dist_phl, x_new, y_new))

      all_cells <- all_cells %>% select(-c(dist_phl))

      # ADJUSTING THE RADIUS OF PHLOEM CELLS:
      # We adjust the radius of the phloem cells and companion cells:
      all_cells <- all_cells %>%
        mutate(radius = ifelse(type == "phloem", radius_phloem_cell, radius))
      all_cells <- all_cells %>%
        mutate(radius = ifelse(type == "companion_cell", radius_companion_cell, radius))

      # REMOVING STELE CELLS INSIDE THE PHLOEM CELLS:
      # We select the group of phloem cells:
      all_phloem <- all_cells %>% filter(type == "phloem")
      # For each phloem cell:
      for(i in c(1:nrow(all_phloem))){
        # We refine the general table of all cells:
        all_cells <- all_cells %>%
          # We remove the stele cells, in case the distance between their center
          # and the center of the current phloem cell is lower than the radius of the phloem cell:
          filter(!((x-all_phloem$x[i])^2 + (y - all_phloem$y[i])^2 < (all_phloem$radius[i] - stele_d/2)^2 & type == "stele"))
      }

      # MAKING CIRCULAR FRONTIER FOR PHLOEM CELLS:
      # We create a table defining the frontiers for each phloem cells:
      x_cir <- seq(-0.95,0.95,0.95/4)
      y_p <- sqrt(1-x_cir^2)
      y_m <- -sqrt(1-x_cir^2)
      phloem_frontier <- NULL
      for (i_phloem in 1:nrow(all_phloem)) {
        tmp <- all_phloem[i_phloem,]
        cir <- tibble(x = rep(x_cir,2), y = c(y_p,y_m))
        cir <- cir*abs(tmp$radius)/2
        cir$x <- cir$x + tmp$x
        cir$y <- cir$y + tmp$y
        phloem_frontier <- rbind(phloem_frontier, data.frame(angle = tmp$angle,
                                                             radius = phl$r[1],
                                                             x = cir$x,
                                                             y = cir$y,
                                                             id_layer = 202,
                                                             id_cell = 1,
                                                             type = "phloem",
                                                             order = 1.5,
                                                             id_group = k))
        # We increment the id_group:
        k <- k + 1
      }

      # add phl frontier
      all_cells <- all_cells[all_cells$type != "phloem",]
      all_cells <- rbind(all_cells, phloem_frontier)
      # all_cells$id_group[all_cells$type == "phloem"] <- all_cells$id_group[all_cells$type == "phloem" & all_cells$id_group != 0] + k_max_cortex
      all_cells$id_group[all_cells$type == "phloem"] <- all_cells$id_group[all_cells$type == "phloem" & all_cells$id_group != 0] + k_max
    }

  # CASE 2: DICOT PLANT
  #####################

  } else if(plant_type == 2){ # DICOT
    xyl <- data.frame(r=numeric(2), d=numeric(2))
    xyl$r <- c(0, max(all_cells$radius[all_cells$type == "stele"]))
    xyl$d <- c(params$value[params$type == "max_size" & params$name == "xylem"], layers$cell_diameter[layers$name == "stele"])

    # Get the cells in between
    fit <- stats::lm(d ~ r, data=xyl)$coefficients
    rnew <- xyl$r[1]
    i <- 1
    rmax <- xyl$r[2]
    dmin <- xyl$d[2]
    keep_going <- T
    while(keep_going){
      xyl <- xyl %>% arrange(r)

      rnew <- xyl$r[i] + xyl$d[i] #+ (xyl$r[2]/10)
      dnew <- fit[1] + rnew*fit[2]
      while(rnew+(dnew/2) > rmax-(dmin/2)){
        rnew <- rnew - 0.05
        dnew <- dnew - 0.05
        keep_going = F
      }
      xyl <- rbind(xyl, data.frame(r = rnew,d = dnew))
      i <- i+1
    }
    xyl <- xyl %>% arrange(r)%>%
      filter(d > 0)
    while(xyl$d[nrow(xyl)] >= xyl$d[nrow(xyl)-1]){
      xyl$d[nrow(xyl)] <- xyl$d[nrow(xyl)] - 0.04
      xyl$r[nrow(xyl)] <- xyl$r[nrow(xyl)] - 0.02
    }
    all_xylem <- NULL
    i <- 1
    angle_seq <- seq(from = 0, to = (2*pi), by = (2 * pi) / n_xylem_files)
    x <- center + (xyl$r[1] * cos(angle_seq[1]))
    y <- center + (xyl$r[1] * sin(angle_seq[1]))
    all_xylem <- rbind(all_xylem, data.frame(x = x,
                                             y = y,
                                             d = xyl$d[1],
                                             angle = angle_seq[1],
                                             id_group = k))
    all_cells <- rbind(all_cells, data.frame(
      angle = angle_seq[1],
      radius = xyl$r[1],
      x = x,
      y = y,
      id_layer = 20,
      id_cell = 1,
      type = "xylem",
      order = 1.5,
      id_group = k
    ))

    k <- k+1
    for(angle in angle_seq){
      x <- center + (xyl$r[-1] * cos(angle))
      y <- center + (xyl$r[-1] * sin(angle))
      all_xylem <- rbind(all_xylem, data.frame(x = x,
                                               y = y,
                                               d = xyl$d[-1],
                                               angle = angle,
                                               id_group = k))
      all_cells <- rbind(all_cells, data.frame(
        angle = angle,
        radius = xyl$r[-1],
        x = x,
        y = y,
        id_layer = 20,
        id_cell = 1,
        type = "xylem",
        order = 1.5,
        id_group = k
      )
      )
      i <- k+1
    }
    # Phloem vessels are built between xylem ones
    phl <- data.frame(r = max(all_cells$radius[all_cells$type == "stele"]) - (params$value[params$type == "cell_diameter" & params$name == "stele"])/2,
                      d = params$value[params$type == "cell_diameter" & params$name == "stele"])
    angle_seq_ph <- seq(from = ((2 * pi) / n_xylem_files ) /2, to = (2*pi), by = (2 * pi) / n_xylem_files)
    for(angle in angle_seq_ph){
      x1 <- center + (phl$r[1] * cos(angle))
      y1 <- center + (phl$r[1] * sin(angle))
      #Find the closest stele cell and assign it as a phloem vessel
      all_cells <- all_cells %>%
        mutate(type = as.character(type)) %>%
        mutate(dist_phl = sqrt((x-x1)^2 + (y-y1)^2)) %>%
        mutate(dist_phl = ifelse(type == "stele", dist_phl, 100)) %>%
        mutate(type = ifelse(dist_phl == min(dist_phl), "phloem", type))

      # Get the compagnion cells
      for (compa in 1:2) {
        all_cells <- all_cells %>%
          mutate(type = as.character(type)) %>%
          mutate(dist_phl = sqrt((x-x1)^2 + (y-y1)^2)) %>%
          mutate(dist_phl = ifelse(type == "stele", dist_phl, 100)) %>%
          mutate(type = ifelse(dist_phl == min(dist_phl), "companion_cell", type))
      }

      for (compa in 1:6) {
        all_cells <- all_cells %>%
          mutate(type = as.character(type)) %>%
          mutate(dist_phl = sqrt((x-x1)^2 + (y-y1)^2)) %>%
          mutate(dist_phl = ifelse(type == "stele", dist_phl, 100)) %>%
          mutate(type = ifelse(dist_phl == min(dist_phl), "cambium", type))
      }
      all_cells <- all_cells%>% select(-dist_phl)

    }

    # We change the identity of stele cells to be replaced by xylem cells:
    for(i in c(1:nrow(all_xylem))){
      all_cells <- all_cells %>%
        filter(!((x-all_xylem$x[i])^2 + (y - all_xylem$y[i])^2 < (all_xylem$d[i]/1.5)^2 & type == "stele")) # find the cells inside the xylem poles and remove them
    }
    all_xylem <- all_cells[all_cells$type == "xylem",]%>%
      mutate(radius = round(radius,2),
             d = c(xyl$d[1], rep(xyl$d[2:length(xyl$d)],
                                 max(unique(all_cells$id_group[all_cells$type == "xylem"]))-1)))# %>%
    # Remove xylem cell that are to close to each other
    # filter(!duplicated(angle), !duplicated(radius))
    all_cells <- all_cells[all_cells$type != "xylem",]

    # Make circular frontier for metaxylem:
    x_cir <- seq(-0.95,0.95,0.95/4)
    y_p <- sqrt(1-x_cir^2)
    y_m <- -sqrt(1-x_cir^2)
    xyl_frontier <- NULL
    for (i_xyl in 1:nrow(all_xylem)) {
      tmp <- all_xylem[i_xyl,]
      cir <- tibble(x = rep(x_cir,2), y = c(y_p,y_m))
      cir <- cir*abs(tmp$d*0.8)/2
      cir$x <- cir$x + tmp$x
      cir$y <- cir$y + tmp$y
      xyl_frontier <- rbind(xyl_frontier, data.frame(angle = tmp$angle,
                                                     radius = tmp$radius,
                                                     x = cir$x,
                                                     y = cir$y,
                                                     id_layer = 20,
                                                     id_cell = 1,
                                                     type = "xylem",
                                                     order = 1.5,
                                                     id_group = k))

      k <- k + 1
    }
    all_cells <- rbind(all_cells, xyl_frontier)
    # all_cells$id_group[all_cells$type == "xylem"] <- all_cells$id_group[all_cells$type == "xylem" & all_cells$id_group != 0] + k_max_cortex
    all_cells$id_group[all_cells$type == "xylem"] <- all_cells$id_group[all_cells$type == "xylem" & all_cells$id_group != 0] + k_max

    if(length(params$value[params$name == "pith" & params$type == "layer_diameter"]) > 0){
      all_cells <- make_pith(all_cells, params, center)
    }
  }

  #/////////////////////////////////////////////////////////////////////////////////////////////////////

  # Eventually, we reset the cell id of all cells:
  all_cells$id_cell <- c(1:nrow(all_cells))

  return(all_cells)
}
