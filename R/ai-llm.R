#'
#' @export
ai.get_ollama_models <- function(models=NULL, size=NULL) {
  available.models <- system("ollama list | tail -n +2 | cut -d' ' -f 1", intern=TRUE)

  models.sizes <- system("ollama list | tail -n +2 | tr -s ' ' | cut -d' ' -f 3", intern=TRUE)
  models.sizes <- as.numeric(models.sizes)
  models.sizes <- ifelse(models.sizes < 100, models.sizes, models.sizes/1000)
  names(models.sizes) <- available.models
  table(models.sizes < 5)

  ##hist(models.sizes, breaks=50)
  if(!is.null(models) && !any(models=="OLLAMA_MODELS")) {
    available.models <- intersect(models,available.models)
  }

  msize <- models.sizes[available.models] 
  if(!is.null(size) && size=="S") {
    sel <- which( msize <= 3 )
    available.models <- available.models[sel]
  }
  if(!is.null(size) && size=="M") {
    sel <- which( msize > 3 & msize <= 6)
    available.models <- available.models[sel]
  }
  if(!is.null(size) && size=="L") {
    msize <- models.sizes[available.models] 
    sel <- which( msize > 6 )
    available.models <- available.models[sel]
  }
  
  return(available.models)
}

OLLAMA_MODELS = ai.get_ollama_models()
DEFAULT_LLM = "gpt-5-nano"
DEFAULT_LLM = NULL

REMOTE_MODELS <- c(
  "openai:gpt-5-nano",
  "xai:grok-4-fast-non-reasoning", 
  "groq:llama-3.1-8b-instant",
  "groq:meta-llama/llama-4-scout-17b-16e-instruct",  
  "groq:openai/gpt-oss-20b",
  "groq:openai/gpt-oss-120b",
  "google:gemini-2.5-flash-lite"
)

if(0) {
  model="gpt-5-nano";prompt=NULL
  model="gemma3:1b";prompt=NULL
  model="grok-4-fast-non-reasoning";prompt=NULL
}

#' @export
ai.get_remote_models <- function(models=NULL) {
  keys <- NULL
  dbg("[ai.get_remote_models] input models = ",models)
  dbg("[ai.get_remote_models] len.models = ",length(models))
  if(!is.null(models) && "REMOTE_MODELS" %in% models) {
    models <- unique(c(models, REMOTE_MODELS))
  }
  if(is.null(models)) {
    models <- REMOTE_MODELS
  }
  if (Sys.getenv("OPENAI_API_KEY")!="") keys <- c(keys,"gpt-","openai:")
  if (Sys.getenv("XAI_API_KEY")!="") keys <- c(keys,"grok-","xai:")
  if (Sys.getenv("GROQ_API_KEY")!="") keys <- c(keys,"groq:")
  if (Sys.getenv("GEMINI_API_KEY")!="") keys <- c(keys,"gemini-","google:")
  if (Sys.getenv("OLLAMA_REMOTE")!="") keys <- c(keys,"remote:.*")

  dbg("[ai.get_remote_models] keys = ", keys)
  
  if(is.null(models) || length(models)==0 || models[1]=="" ) {
    models <- keys
  } else if(!is.null(keys)) {
    regex <- paste0("^",keys,collapse="|")
    models <- grep(regex,models,value=TRUE)
  } else {
    models <- NULL
  }

  dbg("[ai.get_remote_models] available models = ",models)
  
  models
}

#' @export
ai.get_models <- function(models=NULL) {
  local.models  <- sort(ai.get_ollama_models(models))
  remote.models <- sort(ai.get_remote_models(models))   
  models <- list( local = local.models, remote = remote.models)
  return(models)
}

#' @export
ai.model_is_available <- function(model) {
  model %in% ai.get_models(models=model) 
}

#' @export
ai.ask <- function(question, model, engine=c("ellmer","tidyprompt")[2]) {
  if(model == "ellmer" && grepl("grok",model)) model <- "tidyprompt"
  if(engine=="ellmer") {
    resp <- ai.ask_ellmer(question=question, model=model, prompt=NULL) 
  }
  if(engine=="tidyprompt") {
    resp <- ai.ask_tidyprompt(question=question, model=model) 
  }
  return(resp)
}

#' @export
ai.ask_ellmer <- function(question, model=DEFAULT_LLM, prompt=NULL) {
  chat <- NULL
  if(inherits(model, "Chat")) {
    chat <- model
  } else if(is.character(model)) {
    if (model %in% OLLAMA_MODELS || grepl("^ollama:",model) ) {
      model1 <- sub("^ollama:","",model)
      chat <- ellmer::chat_ollama(model = model1, system_prompt = prompt)
    } else if (grepl("^gpt|^openai:",model) && Sys.getenv("OPENAI_API_KEY")!="") {
      message("warning: using remote GPT model:", model)
      model1 <- sub("^openai:","",model)
      chat <- ellmer::chat_openai(
        model = model1, system_prompt = prompt,
        api_key = Sys.getenv("OPENAI_API_KEY") )
    } else if (grepl("^grok|^xai:",model) && Sys.getenv("XAI_API_KEY")!="") {
      model1 <- sub("^xai:","",model)
      chat <- ellmer::chat_openai(
        model = model1, system_prompt = prompt,
        api_key = Sys.getenv("XAI_API_KEY"),
        base_url="https://api.x.ai/v1/")
    } else if (grepl("^groq:",model) && Sys.getenv("GROQ_API_KEY")!="") {
      model1 <- sub("groq:","",model)
      chat <- ellmer::chat_groq(
        model = model1, system_prompt = prompt,
        api_key = Sys.getenv("GROQ_API_KEY")
      )
    } else if (grepl("^gemini|^google:",model) && Sys.getenv("GEMINI_API_KEY")!="") {
      model1 <- sub("^google:","",model)
      chat <- ellmer::chat_google_gemini(
        model = model1, system_prompt = prompt,
        api_key = Sys.getenv("GEMINI_API_KEY")
      )
    }    
  }
  
  if(is.null(chat)) {
    message("ERROR. could not create model ", model)
    return(NULL)
  }
  . <- chat$chat(question, echo=FALSE)  
  chat$last_turn()@text
}

ai.ask_tidyprompt <- function(question, model, verbose=0) {
  llm <- NULL
  if(model %in% OLLAMA_MODELS || grepl("^ollama:",model) ) {
    model1 <- sub("^ollama:","",model)
    llm <- tidyprompt::llm_provider_ollama(
      parameters = list(
        model = model1
      )
    )
  } else if(grepl("^remote:",model) ) {
    remotesrv <- Sys.getenv("OLLAMA_REMOTE")
    if(remotesrv=="") {
      message("error: please set OLLAMA_REMOTE")
    }
    if(remotesrv!="") {
      model1 <- sub("^remote:","",model)    
      if(verbose>0) {
        message("connecting to remote ollama server = ",remotesrv)
        message("remote model = ",model1)        
      }
      llm <- tidyprompt::llm_provider_ollama(
        parameters = list(
          model = model1
        ),
        url = paste0("http://",remotesrv,"/api/chat")
      )
    }
  } else if(grepl("^groq:",model)) {
    model2 <- sub("groq:","",model)
    llm <- tidyprompt::llm_provider_groq(
      parameters = list(model = model2)
    )
  } else if(grepl("^grok|^xai:",model)) {
    model2 <- sub("^xai:","",model)
    llm <- tidyprompt::llm_provider_xai(
      parameters = list(model = model2)
    )
  } else if(grepl("^gpt-|^openai:",model)) {
    model2 <- sub("^openai:","",model)
    llm <- tidyprompt::llm_provider_openai(
      parameters = list(model = model2)
    )
  } else if(grepl("^gemini-|^google:",model)) {
    model2 <- sub("^google:","",model)
    llm <- tidyprompt::llm_provider_google_gemini(
      parameters = list(model = model2),
      api_key = Sys.getenv("GEMINI_API_KEY")
    )
  }
  if(is.null(llm)) {
    message("warning. unsupported model: ", model)
    return(NULL)
  }
  
  ##question = "what is P53?"
  if(verbose>0) {
    message("model = ",model)
    message("question = ",question)    
  }

  resp <- NULL
  resp <- question |>
    tidyprompt::send_prompt(
      llm_provider = llm,
      clean_chat_history = TRUE,
      verbose = FALSE,
      return_mode = "only_response"
    )

  ## clean response
  resp <- sub("<think>.*</think>","",resp)

  return(resp)
}


##----------------------------------------------------------------------
##----------------------------------------------------------------------
##----------------------------------------------------------------------

#'
#'
#' @export
ai.genesets_summary <- function(gsets, pheno=NULL, model=DEFAULT_LLM,
                                detail=1, html=FALSE, verbose=1) {
  q <- "Extract the main biological function of this list of gene sets that were found by doing geneset enrichment. Just give the answer. Do not acknowledge."
  if(!is.null(pheno)) q <- paste0(q, "Discuss in relation with the phenotype: '",pheno,"'.")
  if(detail==0) q <- paste(q, "Be very very short.")
  if(detail==1) q <- paste(q, "Describe in one short paragraph.")
  if(detail>=2) q <- paste(q, "Describe in detail.")
  if(html) q <- paste(q, "Use HTML formatting.")  
  if(verbose>0) cat("Question:",q,"... \n")
  ss <- paste(gsets, collapse='; ')
  q <- paste(q, "These are the genesets: <list>",ss,"</list>. ")
  r <- ai.ask(q, model=model)
  #r <- ai.ask(q, model="gemma3:270m")
  #r <- ai.ask(q, model="gemma3:1b")    
  return(r)
}

num=3
#' @export
ai.genesets_keywords <- function(gsets, num=3, pheno=NULL, model=DEFAULT_LLM) {
  ss <- paste(gsets, collapse='; ')
  q <- paste0("Extract ",num," keywords describing the following collection of gene sets. ")
  q <- paste0(q, "These are the genesets: <list>",ss,"</list>. ")
  r <- ai.ask(q, model=model)
  return(r)
}


