OLLAMA_MODELS <- c("gemma3:4b", "gemma3:27b", "deepseek-r1:1.5b", "deepseek-r1:7b")
GPT_MODELS <- c("gpt-4o-mini")
prompt <- NULL
deep.ask <- function(question, model = "gpt-4o-mini", prompt = NULL) {
  if (model == "gpt-4o-mini") {
    chat <- ellmer::chat_openai(model = "gpt-4o-mini", system_prompt = prompt)
  } else if (model %in% OLLAMA_MODELS) {
    chat <- ellmer::chat_ollama(model = model, system_prompt = prompt)
  } else {
    message("ERROR. unknown model", model)
    message("valid models: ", paste(c(OLLAMA_MODELS, GPT_MODELS), collapse = ", "))
    return(NULL)
  }
  . <- chat$chat(question, echo = FALSE)
  chat$last_turn()@text
}
