#'
#'
#'
#' @export
getPathwayImage <- function(wp, val, sbgn.dir = NULL, as.img = FALSE) {
  img <- NULL
  if (grepl("WP", wp)) img <- wikipathview(wp, val = val)
  if (grepl("SMP", wp)) img <- pathbankview(wp, val = val)
  if (grepl("R-HSA", wp)) img <- getReactomeSVG(wp, val = val)
  ##  if(grepl("R-HSA",wp)) img <- getReactomeSVG.SBGN(wp, val=val, sbgn.dir=sbgn.dir)

  if (is.null(img)) {
    return(NULL)
  }

  if (as.img) {
    img <- list(
      src = normalizePath(img),
      contentType = "image/svg+xml",
      width = "100%", height = "100%", ## actual size: 1040x800
      alt = paste(wp, "pathway (SVG image)")
    )
  }
  return(img)
}


#' @export
getReactomeSVG <- function(wp, val = NULL, as.img = FALSE) {
  require(xml2)

  ## wp="R-HSA-449147"
  url <- paste0("https://reactome.org/ContentService/exporter/diagram/", wp, ".svg")
  destfile <- tempfile(fileext = ".svg")
  down <- tryCatch(
    {
      download.file(url, destfile)
    },
    error = function(w) {
      return(NULL)
    }
  )

  if (as.img) {
    destfile <- list(
      src = normalizePath(destfile),
      contentType = "image/svg+xml",
      width = "100%", height = "100%", ## actual size: 1040x800
      alt = paste("Reactome pathway downloaded from", url)
    )
  }

  return(destfile)
}


#' Download and color nodes in PathBank SVGs
#'
#' Downloads a PathBank SVG file, modifies its XML structure to color pathway nodes
#' according to input values, and returns the path to the modified file.
#'
#' @param pb A character string specifying the PathBank pathway ID.
#' @param val A named numeric vector where names correspond to pathway element IDs,
#'            and values represent the level of regulation (e.g., fold-change or intensity).
#'            Positive values indicate up-regulation, while negative values indicate down-regulation.
#'
#' @return A character string representing the path to the modified SVG file,
#'         or `NULL` if the download fails.
#'
#' @details
#' The function downloads the SVG representation of a PathBank pathway, modifies its
#' XML content to color pathway nodes based on the provided `val` parameter:
#' - **Positive values**: Nodes are colored red, with the intensity increasing as the value grows.
#' - **Negative values**: Nodes are colored blue, with the intensity increasing as the absolute value grows.
#'
#' Opacity is determined by the square root of the absolute value, scaled and capped at a maximum of 1.
#'
#' Nodes with IDs matching the names in `val` are colored accordingly. Nodes without matching IDs
#' or `val` values are left unchanged.
#'
#' @examples
#' \dontrun{
#' pb <- "SMP000001"
#' val <- c("metabolite1" = 2.0, "metabolite2" = -1.5)
#' modified_svg <- pathbankview(pb, val)
#' if (!is.null(modified_svg)) {
#'   print(paste("Modified SVG saved at:", modified_svg))
#' }
#' }
#'
#' @export
pathbankview <- function(pb, val, as.img = FALSE, large_font = TRUE) {
  require(xml2)

  if (large_font) {
    url <- paste0("https://www.pathbank.org/view/", pb, "/download?type=simple_large_font_vector_image")
  } else {
    url <- paste0("https://www.pathbank.org/view/", pb, "/download?type=simple_vector_image")
  }
  destfile <- tempfile(fileext = ".svg")
  down <- tryCatch(
    {
      download.file(url, destfile)
    },
    error = function(w) {
      return(NULL)
    }
  )

  # Read the file line by line
  lines <- readLines(destfile)

  # Use gsub to replace the line
  lines <- gsub(
    'xmlns="http://www.w3.org/2000/svg"',
    'xmlns:svg="http://www.w3.org/2000/svg"',
    lines
  )

  # Write the lines back to the file
  writeLines(lines, destfile)

  # Load the SVG
  doc <- tryCatch(
    {
      read_xml(destfile)
    },
    error = function(w) {
      NULL
    }
  )
  if (is.null(doc)) {
    return(NULL)
  }

  label_nodes <- xml_find_all(doc, ".//text")
  labels <- xml2::xml_text(label_nodes)

  # Find nodes
  g_nodes <- xml_find_all(doc, ".//g")
  g_nodes_labels <- xml_attr(g_nodes, "data-element-id")
  g_nodes_types <- xml_attr(g_nodes, "data-element-type")

  # Find the 'rect' and 'circle' childs nodes of the g_nodes
  circle_nodes <- xml_find_first(g_nodes, ".//circle")
  tspan_nodes <- xml_find_first(g_nodes, ".//tspan")[which(g_nodes_types == "protein")]
  tspan_label <- xml_text(tspan_nodes)
  rect_nodes <- xml_find_first(xml2::xml_parent(xml2::xml_parent(tspan_nodes)), ".//rect")

  if (!is.null(val)) {
    if (sum(names(val) %in% g_nodes_labels) > 0) {
      # Circles (metabolites)
      found_indexes <- which(g_nodes_labels %in% names(val))
      g_nodes_labels <- g_nodes_labels[found_indexes]
      circle_nodes <- circle_nodes[found_indexes]
      val_circle <- val[g_nodes_labels]
      rr <- as.character(round(66 * pmin(1, abs(val_circle / 2.0))**0.5))
      rr <- stringr::str_pad(rr, width = 2, pad = "0")
      colors <- ifelse(val_circle > 0, paste0("#ff0000", rr), paste0("#0055ff", rr))
      xml_attr(circle_nodes, "fill") <- colors
      # Rects (proteins)
      found_indexes <- which(paste0("px:", tspan_label) %in% names(val))
      tspan_label <- tspan_label[found_indexes]
      rect_nodes <- rect_nodes[found_indexes]
      val_rect <- val[paste0("px:", tspan_label)]
      rr <- as.character(round(66 * pmin(1, abs(val_rect / 2.0))**0.5))
      rr <- stringr::str_pad(rr, width = 2, pad = "0")
      colors <- ifelse(val_rect > 0, paste0("#ff0000", rr), paste0("#0055ff", rr))
      xml_attr(rect_nodes, "fill") <- colors
    }
  }

  # rect_nodes <- xml_find_first(g_nodes, ".//circle")

  write_xml(doc, destfile)

  # Read the file line by line
  lines <- readLines(destfile)

  # Use gsub to replace the line
  lines <- gsub(
    'xmlns:svg="http://www.w3.org/2000/svg"',
    'xmlns="http://www.w3.org/2000/svg" xmlns:svg="http://www.w3.org/2000/svg"',
    lines
  )

  # Write the lines back to the file
  writeLines(lines, destfile)

  if (as.img) {
    destfile <- list(
      src = normalizePath(destfile),
      contentType = "image/svg+xml",
      width = "100%", height = "100%", ## actual size: 1040x800
      alt = "pathbank pathway SVG"
    )
  }

  return(destfile)
}

#' Download and color nodes in WikiPathways SVGs
#'
#' Downloads a WikiPathways SVG file, modifies its XML structure to color pathway nodes
#' according to regulation values, and returns the path to the modified file.
#'
#' @param wp A character string specifying the WikiPathways pathway ID.
#' @param val A named numeric vector where names correspond to pathway element IDs,
#'            and values represent the level of regulation (e.g., fold-change or intensity).
#'            Positive values indicate up-regulation, while negative values indicate down-regulation.
#'
#' @return A character string representing the path to the modified SVG file,
#'         or `NULL` if the download fails.
#'
#' @details
#' The function downloads the SVG representation of a WikiPathways pathway, modifies its
#' XML content to color pathway nodes based on the provided `val` parameter:
#' - **Positive values**: Nodes are colored red, with the intensity increasing as the value grows.
#' - **Negative values**: Nodes are colored blue, with the intensity increasing as the absolute value grows.
#'
#' Opacity is determined by the square root of the absolute value, scaled and capped at a maximum of 1.
#'
#' Nodes with IDs matching the names in `val` are colored accordingly. Nodes without matching IDs
#' or `val` values are left unchanged.
#' @examples
#' \dontrun{
#' wp <- "WP554"
#' val <- c("gene1" = 1.8, "gene2" = -2.2)
#' modified_svg <- wikipathwaysview(wp, val)
#' if (!is.null(modified_svg)) {
#'   print(paste("Modified SVG saved at:", modified_svg))
#' }
#' }
#'
#' @export
wikipathview <- function(wp, val, as.img = FALSE) {
  require(xml2)

  isClassic <- FALSE
  url <- paste0("https://www.wikipathways.org/wikipathways-assets/pathways/", wp, "/", wp, ".svg")
  destfile <- tempfile(fileext = ".svg")
  down <- tryCatch(
    {
      download.file(url, destfile)
    },
    error = function(w) {
      tryCatch(
        {
          isClassic <<- TRUE
          url <- paste0("https://classic.wikipathways.org//wpi/wpi.php?action=downloadFile&type=svg&pwTitle=Pathway:", wp)
          download.file(url, destfile)
        },
        error = function(w) {
          return(NULL)
        }
      )
    }
  ) |> is.null()

  if (down) {
    return(NULL)
  }

  # Read the file line by line
  lines <- readLines(destfile)

  # Use gsub to replace the line
  lines <- gsub(
    'xmlns="http://www.w3.org/2000/svg"',
    'xmlns:svg="http://www.w3.org/2000/svg"',
    lines
  )

  # Write the lines back to the file
  writeLines(lines, destfile)

  # Load the SVG
  doc <- tryCatch(
    {
      read_xml(destfile)
    },
    error = function(w) {
      NULL
    }
  )
  if (is.null(doc)) {
    return(NULL)
  }

  # Find all 'text' elements
  label_nodes <- xml_find_all(doc, ".//text")
  labels <- xml2::xml_text(label_nodes)

  # Find the 'a' parent nodes of the label nodes
  parent_nodes <- lapply(label_nodes, xml_parent)
  parent_paths <- lapply(parent_nodes, xml_path)
  duplicated_parents <- duplicated(parent_paths) | duplicated(parent_paths, fromLast = TRUE)
  if (any(duplicated_parents)) {
    duplicated_label_nodes_indices <- which(duplicated_parents)
  } else {
    duplicated_label_nodes_indices <- NULL
  }
  # Remove duplicated parents
  if (!is.null(duplicated_label_nodes_indices)) {
    label_nodes <- label_nodes[-duplicated_label_nodes_indices]
    labels <- labels[-duplicated_label_nodes_indices]
  }
  a_nodes <- xml_parent(label_nodes)

  # Find the 'rect' children of the 'a' nodes
  if (isClassic) {
    rect_nodes <- xml_find_first(a_nodes, ".//path")
  } else {
    rect_nodes <- xml_find_first(a_nodes, ".//rect")
  }

  if (all(is.na(rect_nodes))) {
    val <- NULL
  }

  if (!is.null(val)) {
    if (sum(names(val) %in% toupper(labels)) > 0) {
      found_indexes <- which(toupper(labels) %in% names(val))
      labels <- labels[found_indexes]
      rect_nodes <- rect_nodes[found_indexes]
      val <- val[toupper(labels)]
      rr <- as.character(round(66 * pmin(1, abs(val / 2.0))**0.5))
      rr <- stringr::str_pad(rr, width = 2, pad = "0")
      colors <- ifelse(val > 0, paste0("#ff0000", rr), paste0("#0055ff", rr))
      if (isClassic) {
        current_style <- xml_attr(rect_nodes, "style")
        new_style <- paste0(current_style, " fill:", colors, ";")
        lapply(seq_along(new_style), function(x) {
          xml_set_attr(rect_nodes[x], "style", new_style[x])
        })
      } else {
        xml_attr(rect_nodes, "fill") <- colors
      }
    }
  }

  write_xml(doc, destfile)

  # Read the file line by line
  lines <- readLines(destfile)

  # Use gsub to replace the line
  lines <- gsub(
    'xmlns:svg="http://www.w3.org/2000/svg"',
    'xmlns="http://www.w3.org/2000/svg" xmlns:svg="http://www.w3.org/2000/svg"',
    lines
  )

  # Write the lines back to the file
  writeLines(lines, destfile)

  if (as.img) {
    destfile <- list(
      src = normalizePath(destfile),
      contentType = "image/svg+xml",
      width = "100%", height = "100%", ## actual size: 1040x800
      alt = "wikipathway SVG"
    )
  }

  return(destfile)
}


#'
#'
#' @export
getReactomeSVG.SBGN <- function(pathway.id, val, sbgn.dir, as.img = FALSE) {
  suppressMessages(require(SBGNview)) ## slow!! but needed!!!

  dbg("[getReactomeSVG] pathway.id = ", pathway.id)

  ## this is a trick. the original object in SBGNview.data was 700MB!!
  #  sbgn.dir <- pgx.system.file("sbgn/", package = "pathway")
  #  sbgn.dir <- normalizePath(sbgn.dir) ## absolute path
  sbgn.xmls <- dir(sbgn.dir, ".sbgn")
  names(sbgn.xmls) <- sbgn.xmls

  ## We temporarily switch the working directory to always readable
  ## TMP folder
  curwd <- getwd()
  tmpdir <- tempdir()
  setwd(tmpdir)

  obj <- tryCatch(
    {
      SBGNview::SBGNview(
        gene.data = val,
        gene.id.type = "SYMBOL",
        sbgn.dir = sbgn.dir,
        input.sbgn = pathway.id,
        output.file = "reactome",
        output.formats = c("svg")
      )
    },
    error = function(w) {
      SBGNview::SBGNview(
        gene.data = NULL,
        gene.id.type = "SYMBOL",
        sbgn.dir = sbgn.dir,
        input.sbgn = pathway.id,
        output.file = "reactome",
        output.formats = c("svg")
      )
    }
  )
  if (class(obj) == "SBGNview") {
    try(print(obj))
  }
  Sys.sleep(0.2) ## wait for graph

  ## back to previous working folder
  setwd(curwd)

  #  imgfile <- "/tmp/hsa00010.png"
  #  imgfile <- file.path(tmpdir, paste0("reactome_", pathway.id, ".png"))
  svgfile <- file.path(tmpdir, paste0("reactome_", pathway.id, ".svg"))

  file.exists(svgfile)
  if (!file.exists(svgfile)) {
    return(NULL.IMG)
  }

  ## ## parse image dimensions from file
  ## img.dim <- NULL
  ## if (grepl("png|PNG", imgfile)) img.dim <- dim(png::readPNG(imgfile))[1:2]
  ## if (grepl("jpg|JPG", imgfile)) img.dim <- dim(jpeg::readJPEG(imgfile))[1:2]
  ## img.dim

  if (as.img) {
    imgfile <- list(
      ## src = imgfile,
      src = svgfile,
      contentType = "image/svg+xml",
      ## width = img.dim[2], height = img.dim[1], ## actual size
      alt = "reactome pathway (SVG)"
    )
  }

  imgfile
}
