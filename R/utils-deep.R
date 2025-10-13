
ai.get_ollama_models <- function() {
  models <- system("ollama list | tail -n +2 | cut -d' ' -f 1", intern=TRUE)
  return(models)
}

OLLAMA_MODELS = ai.get_ollama_models()
OLLAMA_MODELS

GPT_MINI_MODELS = c("gpt-5-nano")
DEFAULT_LLM = "gpt-5-nano"

model="gpt-4o-mini";prompt=NULL
model="gpt-oss:20b";prompt=NULL
model="gemma3:1b";prompt=NULL
model="gemma3:270m";prompt=NULL
model="qwen3:0.6b";prompt=NULL
model="qwen3:1.7b";prompt=NULL

ai.ask <- function(question, model=DEFAULT_LLM, prompt=NULL) {
  if (grepl("^gpt",model)) {
    if(!model %in% GPT_MINI_MODELS) {
      message("warning: using alarge GPT model:", model)
    }
    chat <- ellmer::chat_openai(model = model, system_prompt = prompt)
  } else if (model %in% OLLAMA_MODELS) {
    chat <- ellmer::chat_ollama(model = model, system_prompt = prompt)
  } else {
    message("ERROR. unknown model", model)
    message("valid models: ", paste(c(OLLAMA_MODELS, GPT_MODELS), collapse = ", "))
    return(NULL)
  }
  . <- chat$chat(question, echo=FALSE)  
  chat$last_turn()@text
}

ai.genesets_summary <- function(gsets, pheno=NULL, model=DEFAULT_LLM) {
  ss <- paste(gsets, collapse='; ')
  q <- paste0("Be very short. Extract the main biological function of the following gene sets. ")
  q <- paste0(q, "These are the genesets: <list>",ss,"</list>. ")
  if(!is.null(pheno)) q <- paste0(q, "Discuss in relation with the phenotype: '",pheno,"'.")
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


