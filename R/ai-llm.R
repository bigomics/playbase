#'
#' @export
ai.get_ollama_models <- function(models=NULL) {
  available.models <- system("ollama list | tail -n +2 | cut -d' ' -f 1", intern=TRUE)
  if(!is.null(models)) available.models <- intersect(models,available.models)
  return(available.models)
}

OLLAMA_MODELS = ai.get_ollama_models()
DEFAULT_LLM = "gpt-5-nano"
DEFAULT_LLM = NULL

if(0) {
  model="gpt-5-nano";prompt=NULL
  model="gemma3:1b";prompt=NULL
  model="grok-4-fast-non-reasoning";prompt=NULL
}

#' @export
ai.get_remote_models <- function(models=NULL) {
  keys <- NULL

  dbg("[ai.get_remote_models] models = ",models)
  dbg("[ai.get_remote_models] len.models = ",length(models))
  dbg("[ai.get_remote_models] OPENAI_API_KEY = ",Sys.getenv("OPENAI_API_KEY"))
  dbg("[ai.get_remote_models] XAI_API_KEY = ",Sys.getenv("XAI_API_KEY"))
  dbg("[ai.get_remote_models] GROQ_API_KEY = ",Sys.getenv("GROQ_API_KEY"))
  dbg("[ai.get_remote_models] GEMINI_API_KEY = ",Sys.getenv("GEMINI_API_KEY"))  

  if (Sys.getenv("OPENAI_API_KEY")!="") keys <- c(keys,"gpt-.*")
  if (Sys.getenv("XAI_API_KEY")!="") keys <- c(keys,"grok-.*")
  if (Sys.getenv("GROQ_API_KEY")!="") keys <- c(keys,"groq:.*")
  if (Sys.getenv("GEMINI_API_KEY")!="") keys <- c(keys,"gemini-.*")

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

#' @export
ai.get_models <- function(models=NULL) {
  local.models <- ai.get_ollama_models(models)
  remote.models <- ai.get_remote_models(models)   
  if(!is.null(models)) {
    models <- models[ models %in% c(local.models, remote.models)]
  } else {
    models <- c(local.models, remote.models)
  }
  return(models)
}

#' @export
ai.model_is_available <- function(model) {
  model %in% ai.get_models(models=model) 
}

#' @export
ai.ask <- function(question, model=DEFAULT_LLM, prompt=NULL) {
  chat <- NULL
  if(inherits(model, "Chat")) {
    chat <- model
  } else if(is.character(model)) {
    if (model %in% OLLAMA_MODELS) {
      chat <- ellmer::chat_ollama(model = model, system_prompt = prompt)
    } else if (grepl("^gpt",model) && Sys.getenv("OPENAI_API_KEY")!="") {
      message("warning: using remote GPT model:", model)
      chat <- ellmer::chat_openai(
        model = model, system_prompt = prompt,
        api_key = Sys.getenv("OPENAI_API_KEY") )
    } else if (grepl("^grok",model) && Sys.getenv("XAI_API_KEY")!="") {
      chat <- ellmer::chat_openai(
        model = model, system_prompt = prompt,
        api_key = Sys.getenv("XAI_API_KEY"),
        base_url="https://api.x.ai/v1/")
    } else if (grepl("^groq",model) && Sys.getenv("GROQ_API_KEY")!="") {
      model <- sub("groq:","",model)
      chat <- ellmer::chat_groq(
        model = model, system_prompt = prompt,
        api_key = Sys.getenv("GROQ_API_KEY")
      )
    } else if (grepl("^gemini",model) && Sys.getenv("GEMINI_API_KEY")!="") {
      chat <- ellmer::chat_google_gemini(
        model = model, system_prompt = prompt,
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


