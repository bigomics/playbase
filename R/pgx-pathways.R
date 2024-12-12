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

pathbankview <- function(pb, val) {
  require(xml2)

  url <- paste0("https://www.pathbank.org/view/", pb, "/download?type=simple_vector_image")
  destfile <- tempfile(fileext = ".svg")
  down <- download.file(url, destfile)

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

  # Find the 'rect' childs nodes of the g_nodes
  rect_nodes <- xml_find_first(g_nodes, ".//circle")

  if (!is.null(val)) {
    if (sum(names(val) %in% g_nodes_labels) > 0) {
      found_indexes <- which(g_nodes_labels %in% names(val))
      g_nodes_labels <- g_nodes_labels[found_indexes]
      rect_nodes <- rect_nodes[found_indexes]
      val <- val[g_nodes_labels]
      rr <- as.character(round(66 * pmin(1, abs(val / 2.0))**0.5))
      rr <- stringr::str_pad(rr, width = 2, pad = "0")
      colors <- ifelse(val > 0, paste0("#ff0000", rr), paste0("#0055ff", rr))
      xml_attr(rect_nodes, "fill") <- colors
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

wikipathview <- function(wp, val) {
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

  return(destfile)
}
