##
## This file is part of the Omics Playground project.
## Copyright (c) 2018-2026 BigOmics Analytics SA. All rights reserved.
##


IMAGE_MODELS <- c(
  "gemini-3-pro-image-preview",
  "gemini-2.5-flash-image",
  "xai:grok-imagine-image"
)

#' Return list of available remote models. 
#'
#' @export
ai.get_image_models <- function(models=NULL) {
  keys <- NULL
  if(!is.null(models) && "IMAGE_MODELS" %in% models) {
    models <- unique(c(models, IMAGE_MODELS))
  }
  if(is.null(models)) {
    models <- IMAGE_MODELS
  }
  if (Sys.getenv("OPENAI_API_KEY")!="") keys <- c(keys,"gpt-","openai:")
  if (Sys.getenv("XAI_API_KEY")!="") keys <- c(keys,"grok-","xai:")
  if (Sys.getenv("GROQ_API_KEY")!="") keys <- c(keys,"groq:")
  if (Sys.getenv("GEMINI_API_KEY")!="") keys <- c(keys,"gemini-","google:")
  if (Sys.getenv("OLLAMA_REMOTE")!="") keys <- c(keys,"remote:.*")
  
  if(is.null(models) || length(models)==0 || models[1]=="" ) {
    models <- keys
  } else if(!is.null(keys)) {
    regex <- paste0("^",keys,collapse="|")
    models <- grep(regex,models,value=TRUE)
  } else {
    models <- NULL
  }

  models
}



#' Generate image with Gemini (aka Nano Banana). Note this model
#' handles very large prompts correctly.
#'
#' @export
ai.create_image_gemini <- function(prompt,
                                   model = "gemini-2.5-flash-image",
                                   api_key = Sys.getenv("GEMINI_API_KEY"),
                                   format = c("file","base64","raw")[1], 
                                   filename = "image.png",
                                   aspectRatio = "16:9", imageSize = "1K",
                                   base_url = "https://generativelanguage.googleapis.com/v1beta"
                                   ) {

  assertthat::assert_that(assertthat::is.string(prompt), assertthat::noNA(prompt))
  assertthat::assert_that(assertthat::is.string(model), assertthat::noNA(model))
  assertthat::assert_that(assertthat::is.string(api_key), assertthat::noNA(api_key))
  require(dplyr)
  
  if (nchar(api_key) == 0) {
    stop("GEMINI_API_KEY environment variable is not set", call. = FALSE)
  }

  message("calling gemini image (warning: $0.134 per image)")
  model <- sub("^google:","",model)
  #url <- glue::glue("{base_url}/models/{model}:generateContent?key={api_key}")
  url <- glue::glue("{base_url}/models/{model}:generateContent")

  headers <- c(
    `x-goog-api-key` = api_key,
    `Content-Type` = "application/json"
  )

  body <- list(
    contents = list(
      list(
        parts = list(
          list(text = prompt)
        )
      )
    ),
    generationConfig = list(
      responseModalities = list("TEXT", "IMAGE"),
      imageConfig =  list( aspectRatio = aspectRatio, imageSize = imageSize )      
    )
  )

  if(grepl("gemini-2.5",model)) {
    ##body$generationConfig$imageConfig <- list( aspectRatio = "16:9", imageSize = "1K" )
    body$generationConfig$imageConfig <- list(aspectRatio = aspectRatio)
  }
  
  response <- httr::POST(
    url = url,
    httr::add_headers(.headers = headers),
    body = jsonlite::toJSON(body, auto_unbox = TRUE),
    encode = "raw"
  )

  httr::http_type(response)
  if (httr::http_type(response) != "application/json") {
    stop("Gemini API returned unexpected content type", call. = FALSE)
  }

  parsed <- response %>%
    httr::content(as = "text", encoding = "UTF-8") %>%
    jsonlite::fromJSON(flatten = TRUE)

  httr::http_error(response)
  if (httr::http_error(response)) {
    error_msg <- if (!is.null(parsed$error$message)) parsed$error$message else "Unknown error"
    stop(paste0("Gemini API request failed [", httr::status_code(response), "]: ", error_msg),
         call. = FALSE)
  }
  
  parts <- parsed$candidates$content.parts  
  b64 <- NULL
  mimetype <- NULL
  for (part in parts) {
    if (!is.null(part$inlineData.data)) {
      b64 <- part$inlineData.data
      b64 <- head(b64[!is.na(b64)],1)
      mimetype <- part$inlineData.mimeType
      mimetype <- head(mimetype[!is.na(mimetype)],1)
      break()
    }
  }
  
  if(is.null(b64) || length(b64)==0) stop("No image data found in response")

  if(format=="file") {
    raw_image <- base64enc::base64decode(b64)    
    filetype <- sub("jpeg","jpg",sub("image/","",mimetype))
    filename2 <- paste0(sub("[.](jpg|jpeg|png)$","",filename,ignore.case=TRUE),".",filetype)
    filename2
    writeBin(raw_image, filename2)
    message("Saved image to: ", filename2)
    return(invisible(filename2))
  }
  if(format=="raw") {
    raw_image <- base64enc::base64decode(b64)    
    return(invisible(raw_image))
  }
  if(format=="base64") {
    return(invisible(b64))
  }
  stop("return error")
}


#' Generate image with OpenAI's dallE. Note the limitation of the
#' prompt of about 1000 characters. 
#'
#' @export
ai.create_image_openai <- function (prompt, model=NULL, 
                                    size = c("auto","1024x1024")[1], 
                                    aspect_ratio = NULL,
                                    format = c("file","base64","raw"),
                                    filename = "image.png",
                                    api_key = Sys.getenv("OPENAI_API_KEY"),
                                    base_url = "https://api.openai.com/v1",
                                    response_format = NULL,
                                    user = NULL, organization = NULL) 
{
  ##    size <- match.arg(size)
  format <- match.arg(format)
  assertthat::assert_that(assertthat::is.string(prompt), assertthat::noNA(prompt))
  assertthat::assert_that(assertthat::is.string(size), assertthat::noNA(size))
  assertthat::assert_that(assertthat::is.string(format), 
    assertthat::noNA(format))
  if (!is.null(user)) {
    assertthat::assert_that(assertthat::is.string(user), 
      assertthat::noNA(user))
  }
  assertthat::assert_that(assertthat::is.string(api_key), 
    assertthat::noNA(api_key))
  if (!is.null(organization)) {
    assertthat::assert_that(assertthat::is.string(organization), 
      assertthat::noNA(organization))
  }

  if(is.null(model)) stop("must provide model")
  model <- sub("^openai:","",model)

  
  if(grepl("api.x.ai",base_url,fixed=TRUE)) {
    message("calling grok (warning: $0.07 per image)")    
    if(is.null(model)) model <- "grok-2-image-1212"
    size <- NULL
  } else if(grepl("api.openai.com",base_url,fixed=TRUE)) {
    ##if(is.null(model)) model <- NULL
    message("calling openai (warning: $0.05 per image)")    
  } else {
    stop("invalid base_url =",base_url)
  }
  
  url <- glue::glue("{base_url}/images/generations")
  headers <- c(Authorization = paste("Bearer", api_key), 
    `Content-Type` = "application/json")
  #    if (!is.null(organization)) {
  #        headers["Organization"] <- organization
  #    }

  body <- list()
  body[["model"]] <- model
  body[["prompt"]] <- prompt
  body[["n"]] <- 1
  if(!is.null(response_format)) body[["response_format"]] <- "b64_json"
  body[["size"]] <- size
  if(!is.null(aspect_ratio)) body[["aspect_ratio"]] <- aspect_ratio
  body[["user"]] <- user
  response <- httr::POST(url = url, httr::add_headers(.headers = headers), 
    body = body, encode = "json")

  httr::http_type(response)
  if (httr::http_type(response) != "application/json") {
    paste("OpenAI API probably has been changed. Please check online documentation.") %>% 
      stop()
  }

  parsed <- response %>% httr::content(as = "text", encoding = "UTF-8") %>% 
    jsonlite::fromJSON(flatten = TRUE)
  if (httr::http_error(response)) {
    error_msg <- parsed$error
    if(is.list(error_msg)) error_msg <- parsed$error$message
    paste0("OpenAI API request failed [", httr::status_code(response), 
      "]:\n\n", error_msg) %>% stop(call. = FALSE)
  }

  b64 <- parsed$data[['b64_json']]
  if(is.null(b64) || length(b64)==0) stop("No image data found in parsed response")

  if(format == "file") {
    raw_image <- base64enc::base64decode(b64)    
    writeBin(raw_image, filename)
    message("Saved image to: ", filename)
    return(invisible(filename))
  }
  if(format == "raw") {
    raw_image <- base64enc::base64decode(b64)    
    return(invisible(raw_image))
  }
  if(format == "base64") {
    return(invisible(b64))
  }
  stop("return error")
  
}


#' Generate image with Grok (which uses Flux). Note the limitation of
#' the prompt of about 1000 characters.
#'
#' @export
ai.create_image_grok <- function(prompt,
                                 model = "grok-imagine-image", 
                                 model2 = "xai:grok-4-1-fast-non-reasoning",
                                 format = c("file","base64","raw")[1], 
                                 api_key = Sys.getenv("XAI_API_KEY"),
                                 base_url = "https://api.x.ai/v1",
                                 filename = "image.png",
                                 aspect_ratio = "1:1",
                                 user = NULL, organization = NULL) {


  prompt2 <- prompt
  if(nchar(prompt) > 7800) {
    prompt2 <- ai.ask(
      paste("Summarize the following prompt to less than 8000 characters. Retain all information and instructions to faithfully create the image. This is the prompt: ",prompt),
      model = model2)
  }
  nchar(prompt2)

  model <- sub("^grok:|^xai:","",model)
  ai.create_image_openai (
    prompt2, 
    size = "default",
    aspect_ratio = aspect_ratio,
    format = format,
    filename = filename,
    model = model, ##size=NULL,
    base_url = base_url,
    response_format <- "b64_json",
    api_key = api_key,
    user = user,
    organization = organization) 

}


