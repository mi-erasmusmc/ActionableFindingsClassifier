# Copyright 2021 Jan Kors, Erasmus University Medical Center, Rotterdam.

# ----- Script to build classifiers that detect actionable findings in Dutch radiology notes

library(tidyverse)
library(stringr)
library(text2vec)
library(tokenizers)
library(stopwords)
library(caret)
library(pROC)

read_data <- function(input_file) {
  # ----- Read radiology notes
  #       Assumes that ID, indication, clinical history, description, conclusion, and class (indicating absence (NEG)
  #       or presence (POS) of an actionable finding), are stored in a csv file
  data <- read_csv(input_file)
  data <- data %>%
    select(id, indication, clinical, description, conclusion, class)
  
  # ----- Set text field
  data <- data %>%
    mutate(text = if_else(is.na(conclusion), description, conclusion))  # Insert description if conclusion is NA
  
  data$text <- gsub("-\\s+", " ", data$text)  # Remove hyphens as they may disturb sentence splitting
  
  # ----- Set text field without negation
  data <- data %>%
    mutate(text_no_neg = text)
  for (index in 1:length(data$text_no_neg)) {
    if (!is.na(data$text_no_neg[index])) {
      lines <- unlist(tokenize_sentences(data$text_no_neg[index]))
      for (index2 in 1:length(lines)) {
        line <- lines[index2]
        tokens <-
          unlist(tokenize_words(line, lowercase = FALSE))
        if (sum(i <- grep("^geen$|^niet$|^zonder$|^ongewijzigde?$|^onveranderde?$|^bekende?$", tolower(tokens))) > 0) {
          if (min(i) > 1) {
            tokens <- tokens[1:(min(i) - 1)]
            line <- paste(tokens, collapse = " ")
            line <- paste0(line, ".")
          } else {
            line <- ""
          }
          lines[index2] <- line
        }
      }
      note <- paste(lines, collapse = " ")
      if (any(grepl("^\\s+$", note))) {
        note <- ""
      }
      data$text_no_neg[index] <- note
    }
  }
  return(data)
}

# ----- Read training and test data sets
data <- read_data("training_rad_notes.csv")
test <- read_data("test_rad_notes.csv")

# ----- Randomize order of training data
set.seed(123)
rows <- sample(nrow(data))
train <- data[rows, ]

# ----- Results file
datetime <- strftime(Sys.time(), "%Y%m%d_%H%M%S")
ofile <- paste0("models_", datetime, "_Summary.csv")

hdr <- "rm_neg,rm_stop,stem,ngram,AUC_train,AUC_test\n"
capture.output(cat(hdr), file = ofile, append = TRUE)

# ----- Training settings
trctrl <- trainControl(
  summaryFunction = twoClassSummary,
  classProbs = TRUE,
  verboseIter = TRUE,
  savePredictions = TRUE,
  method = "cv",
  number = 10
)

# ----- Batch options
#opt_rm_negation <- c(FALSE, TRUE)
#opt_rm_stopwords <- c(FALSE, TRUE)
#opt_stemming <- c(FALSE, TRUE)
#opt_ngram <- c("unigram", "bigram", "trigram")
# Settings for best training model
opt_rm_negation <- c(TRUE)
opt_rm_stopwords <- c(FALSE)
opt_stemming <- c(TRUE)
opt_ngram <- c("unigram")

# ----- Parameter: remove negation
for (rm_negation in opt_rm_negation) {
  if (rm_negation) {
    text <- train$text_no_comm_no_neg  # TO BE CHANGED!!
    text_test <- test$text_no_neg
  } else {
    text <- train$text
    text_test <- test$text
  }
  
  # ----- Parameter: remove stopwords
  for (rm_stopwords in opt_rm_stopwords) {
    if (rm_stopwords) {
      stopword_list = stopwords(language = "nl")
    } else {
      stopword_list = ""
    }
    
    # ----- Parameter: stemming
    for (stemming in opt_stemming) {
      if (stemming) {
        tok_fun = function(x) {
          tokenize_words(x, stopwords = stopword_list) %>%
            lapply(function(x)
              SnowballC::wordStem(x, language = "nl"))
        }
      } else {
        tok_fun = function(x) {
          tokenize_words(x, stopwords = stopword_list)
        }
      }
      
      # ----- Parameter: n-grams
      for (ngram in opt_ngram) {
        if (ngram == "unigram") {
          ngr <- c(1L, 1L)
        } else if (ngram == "bigram") {
          ngr <- c(1L, 2L)
        } else {
          ngr <- c(1L, 3L)
        }
        
        # ----- Create vocabulary
        it_train <-
          itoken(
            text,
            preprocessor = tolower,
            tokenizer = tok_fun,
            ids = train$id
          )
        vocab <-
          prune_vocabulary(create_vocabulary(it_train, ngram = ngr),
                           term_count_min = 5)
        
        # ----- Construct document-term matrix and rename column names for training set
        vectorizer <- vocab_vectorizer(vocab)
        dtm_train <- create_dtm(it_train, vectorizer)
        dtm_train <- as.matrix(dtm_train)
        colmap <-
          tibble(colname = colnames(dtm_train),
                 colindex = 1:length(colnames(dtm_train)))
        colmap[[2]] <- paste0("V", colmap[[2]])
        colnames(dtm_train) <- colmap[[2]]
        
        # ----- Construct document-term matrix and rename column names for test set
        it_test <-
          itoken(
            text_test,
            preprocessor = tolower,
            tokenizer = tok_fun,
            ids = test$id
          )
        dtm_test <- create_dtm(it_test, vectorizer)
        dtm_test <- as.matrix(dtm_test)
        colmap <-
          tibble(colname = colnames(dtm_test),
                 colindex = 1:length(colnames(dtm_test)))
        colmap[[2]] <- paste0("V", colmap[[2]])
        colnames(dtm_test) <- colmap[[2]]
        
        # Fit random forest model
        set.seed(123)
        model_trained <- train(
          x = dtm_train,
          y = train$class,
          method = "ranger",
          tuneGrid = expand.grid(
            .mtry = seq(5, 30, 5),
            .splitrule = "gini",
            .min.node.size = 1
          ),
          trControl = trctrl,
          metric = "ROC"
        )
        
        
        # ----- Test on validation set
        predictions <-
          predict(model_trained, dtm_test, type = "prob")
        auc_test <- auc(test$class, predictions[, "POS"])
        
        # ----- Save results
        line <-
          paste(
            rm_negation,
            rm_stopwords,
            stemming,
            ngram,
            round(mean(model_trained$resample$ROC), digits = 4),
            round(auc_test, digits = 4),
            sep = ","
          )
        capture.output(cat(paste0(line, "\n")),
                       file = ofile,
                       append = TRUE)
        
        settings <- paste0(
          if_else(rm_negation, "noneg_", "neg_"),
          if_else(rm_stopwords, "nostop_", "stop_"),
          if_else(stemming, "stem_", "nostem_"),
          ngram
        )
        saveRDS(model_trained, file = paste0("model_", settings, ".RDS"))
        saveRDS(vocab, file = paste0("vocab_", settings, ".RDS"))
      }
    }
  }
}