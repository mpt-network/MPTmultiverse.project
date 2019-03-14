

summarize_study <- function(
  model
  , dataset
  , data
  , id = NULL
  , condition = NULL
) {
  
  # catch the function call that was used,
  # and other stuff that should be save along the results
  matched_call <- match.call()
  used_model <- utils::read.table(model, skip = 1, stringsAsFactors = FALSE)
  

  # prepare data ----
  if (missing(data)) {
    data <- as.data.frame(readr::read_csv(dataset))
  }
  
  if(is.null(condition)) {
    data$ExpCond <- "no_condition"
    condition <- "ExpCond"
  }
  
  if(is.null(id)) {
    data$Subject <- 1:nrow(data)
    id <- "Subject"
  }
  
  # Ensure that all variables are character
  data$ExpCond <- as.character(data[[condition]])
  data$Subject <- as.character(data[[id]])
  
  
  # check MPT file
  mpt_model <- TreeBUGS::readEQN(model)
  
  if(!is.data.frame(mpt_model)) {
    "I can't comprehend your .eqn file."
  }
  
  
  # remove extraneous colums and check if all specified columns are present
  # in data
  freq_cols <- MPTmultiverse:::get_eqn_categories(model)
  valid_cols <- c(id, condition, freq_cols)
  check_cols <- valid_cols %in% colnames(data)
  
  if(!all(check_cols)) {
    stop("Variable \"", paste(valid_cols[!check_cols], collapse = ", "), "\" not found in data.frame.")
  }
  
  data <- data[, valid_cols]
  
  
  # Check NAs ----
  nas_found <- unlist(lapply(X = data, FUN = anyNA))
  
  if(any(nas_found)) {
    stop("Variable \"", paste(valid_cols[nas_found], collapse = ", "), "\" contains missing values.")
  }
  
  # Check whether freqencies are integer ----
  not_integer <- unlist(lapply(X = data[, freq_cols], FUN = function(x) {
    any(as.integer(x)!=x)
  }
  ))
  
  if(any(not_integer)) {
    stop("Variable \"", paste(freq_cols[not_integer], collapse = ", "), "\" contains non-integer values.")
  }
  
  # Ensure that id and condition are character, also drops unused levels
  data[[id]] <- as.character(data[[id]])
  data[[condition]] <- as.character(data[[condition]])
  
  # summarize the design of your study
  participants <- table(data[[condition]])
  trees <- MPTmultiverse:::get_eqn_trees(model)
  names(trees) <- freq_cols
  n_per_tree <- vector(mode = "list", length = length(unique(trees)))
  names(n_per_tree) <- unname(unique(trees))
  
  for(i in trees) {
    n_per_tree[[i]] <- rowSums(data[, names(trees[trees == i])])
  }
  
  list(
    model = model
    , dataset = dataset
    , data = data
    , participants = participants
    , n_per_tree = n_per_tree
  )
  
}