
ai.get_ollama_models <- function() {
  models <- system("ollama list | tail -n +2 | cut -d' ' -f 1", intern=TRUE)
  return(models)
}

OLLAMA_MODELS = ai.get_ollama_models()
OLLAMA_MODELS
TINY_MODELS = c("llama3.2:1b","granite4:micro","gemma3:1b")
DEFAULT_LLM = "gpt-5-nano"

if(0) {
  model="gpt-5-nano";prompt=NULL
  model="gemma3:1b";prompt=NULL
  model="grok-4-fast-non-reasoning";prompt=NULL
}

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
ai.genesets_keywords <- function(gsets, num=3, pheno=NULL, model=DEFAULT_LLM) {
  ss <- paste(gsets, collapse='; ')
  q <- paste0("Extract ",num," keywords describing the following collection of gene sets. ")
  q <- paste0(q, "These are the genesets: <list>",ss,"</list>. ")
  r <- ai.ask(q, model=model)
  return(r)
}


