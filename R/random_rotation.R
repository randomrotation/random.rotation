#' Random rotation (including potential flip)
#' 
#' Generate random member of orthogonal group O(n)
#' @param n desired dimension of rotation
#' @keywords random rotation
#' @export
#' @examples
#' random_rotation_matrix_incl_flip(3)
random_rotation_matrix_incl_flip <- function(n)
{
  QR <- qr(matrix(rnorm(n^2), ncol=n))          # A = QR
  M <- qr.Q(QR) %*% diag(sign(diag(qr.R(QR))))  # ensure diag(R) > 0
  return(M)
}

#' Random rotation (without flip)
#' 
#' Generate random member of special orthogonal group SO(n)
#' @param n desired dimension of rotation
#' @keywords random rotation
#' @export
#' @examples
#' random_rotation_matrix(3)
random_rotation_matrix <- function(n)
{
  M <- random_rotation_matrix_incl_flip(n)
  if(det(M)<0) M[,1] <- -M[,1]                  # ensure det(M) = 1
  return(M)
}

#' Equal-weight top-h rotations
#' 
#' Weight of rotation based on cutoff. Name in paper: cut
#' @param R total number of available rotations
#' @param h half-life tuning parameter (cutoff)
#' @keywords random rotation
#' @export
#' @examples
#' weights_cut(20,5)
weights_cut <- function(R, h)
{
  (1:R<=h)/h  # equal-weight for top-h
}

#' Exponential decay in weight for rotations
#' 
#' Weight of rotation based on an exponential decay. Name in paper: exp
#' @param R total number of available rotations
#' @param h half-life tuning parameter
#' @keywords random rotation
#' @export
#' @examples
#' weights_exp(20,5)
weights_exp <- function(R, h)
{
  2^(-(1:R)/h) * (2^(1/h)-1) / ( 1 - 2^(-R/h)) # exponential decay with parameter h	
}

#' Equal-weight all rotations
#' 
#' This is just a special case of weights_cut with h=R. Name in paper: rre. This matches the weighting of the earlier paper.
#' @param R total number of available rotations
#' @param h half-life tuning parameter (ignored for equal-weight)
#' @keywords random rotation
#' @export
#' @examples
#' weights_rre(20,0)
weights_rre <- function(R, h)
{
  weights_cut(R,R)
}

#' Determine height of RF decision tree
#' 
#' Recursively compute tree height from the randomForest R package tree representation
#' @param mytree tree representation from randomForest output rf$forest$treemap
#' @param r internal recursion parameter (not used in user call)
#' @param d internal recursion parameter (not used in user call), represents depth
#' @keywords random rotation
#' @export
#' @examples
#' rf <- randomForest(Species ~ ., iris)
#' tree_height(rf$forest$treemap[,,1])
tree_height <- function(mytree, r=1,d=1)
{
  if(mytree[r,1]==0) return(d);
  return(max(tree_height(mytree,mytree[r,1],d+1), tree_height(mytree,mytree[r,2],d+1)))
}

#' Compute number of trees per rotation
#' 
#' Compute the correct number of trees to grow per rotation, given h and the distribution function
#' @param N number of trees in total
#' @param R number of rotations in total
#' @param h value of tuning parameter
#' @param fn distribution function, currently in weights_{cut, exp, rre}
#' @keywords random rotation
#' @export
#' @examples
#' num_trees_per_rot(5000, 100, 10, fn=weights_cut)
num_trees_per_rot <- function(N, R, h, fn=weights_exp)
{
  a <- round(N*fn(R,h))
  b <- N-sum(a)
  if(b > 0) a[1] <- a[1] + b
  if(b < 0)
  {
    for(i in R:1)
    {
      if(a[i] > 0)
      {
      	a[i] <- a[i] + max(b,-a[i])
      	b <- N-sum(a)
      }
    }
  }
  return(a)
}

#' Check if sorted
#'
#' Determine if an input vector v is sorted (boolean)
#' @param v input vector
#' @keywords sorting
#' @export
#' @examples
#' is_sorted(c(1,2,4,3))
is_sorted <- function(v)
{
  l <- length(v)
  if(l > 1)
    for(i in 2:l)
    {
      if(v[i] < v[i-1]) 
        return(FALSE)
    }
  TRUE
}

#' Create map
#'
#' Create map from input vector v in \[-Inf..Inf\] to a rank in \[0..1\]
#' @param v numeric input vector
#' @keywords mapping
#' @export
#' @examples
#' zero_one_rank_map(c(-9, 0, 4, 1000))
zero_one_rank_map <- function(v)
{
  if(is.unsorted(v, na.rm=TRUE))
    v <- sort(v)

  # enforce -Inf => 0, Inf => 1
  s <- v
  if(v[1] > -Inf)
    s <- c(-Inf, s)
  if(v[length(v)] < Inf)
    s <- c(s, Inf)

  p <- rank(s) - 1
  r <- p/max(p)
  list(s,r)
}

#' Apply map
#'
#' Use an in-sample rank map (v_map) to obtain the ranks of an out-of-sample input vector v
#' @param v numeric input vector
#' @param v_map list of inputs and corresponding mapped values
#' @keywords ranking
#' @export
#' @examples
#' uni_rank(c(1.7, 2.1, 3.5), zero_one_rank_map(0, 1, 4, 9))
uni_rank <- function(v, v_map)
{
  i_a <- sapply(v,FUN=function(x){which(x<=v_map[[1]])[1]})
  i_b <- sapply(v,FUN=function(x){a<-which(x>=v_map[[1]]);l<-length(a);if(l) a[length(a)] else NA})
  .5*(v_map[[2]][i_a]+v_map[[2]][i_b])
}

#' Impute NA values
#'
#' Replace NA values in a vector v with a given value x (default 0.5)
#' @param v numeric input vector
#' @param x numeric replacement value (default 0.5)
#' @keywords imputation
#' @export
#' @examples
#' impute(c(0.2, 0.3, NA, 0.7, 0.9))
impute <- function(v, x=0.5)
{
  if(sum(is.na(v))) 
    v[which(is.na(v))] <- x
  v
}

#' Obtain integer indices
#'
#' Matches list/vector mix of names and/or indices vc to columns in a data frame
#' @param df data frame
#' @param vc vector or list of indices (integer) and/or column names
#' @keywords indexing
#' @export
#' @examples
#' match_columns(df, list("col1", 4, "col9"))
match_columns <- function(df, vc)
{
      c_name <- sapply(vc,FUN=function(x)1:ncol(df)==x)
      c_index <- sapply(vc,FUN=function(x)colnames(df)==x)
      c_comb <- cbind(c_name, c_index)

      if(length(c_comb)) which(rowSums(c_comb)>0) else numeric(0)
}

#' Get X and Y columns
#'
#' Returns columns for X matrix and Y vector (and ignored cols) from the input, given set of in-/exclusions
#' @param df data frame
#' @param response index (integer) or column name of class (response)
#' @param exclude vector or list of indices (integer) and/or column names that should be excluded (default: NULL)
#' @param include  vector or list of indices (integer) and/or column names that should ONLY be included (default: NULL)
#' @keywords helper 
#' @export
#' @examples
#' determine_relevant_columns(df, response="myY")
#' determine_relevant_columns(df, response="myY", exclude=c("col1", "col9"))
#' determine_relevant_columns(df, response=5, include=list(2,3,4, "col10"))
determine_relevant_columns <- function(df, response, exclude=NULL, include=NULL)
{
  c_complete <- 1:ncol(df)
  c_response <- match_columns(df, response)
  
  if(length(c_response) > 1)
    warning("Only single column response is supported; using first matching column and ignoring others.")

  if(length(c_response) < 1)
  {
      warning("At least one valid response column is required; defaulting to last column in data frame.")
      c_response <- ncol(df)
  }

  c_exclude <- integer(0)
  if(length(exclude)>0)
  {
    c_exclude <- match_columns(df, exclude)  
    if(length(c_exclude) < length(exclude))
      warning("Some of the supplied exclusions are duplicates or could not be found in the data frame, will ignore these.")
  }
  c_exclude <- union(c_exclude, setdiff(c_response, c_response[1]))

  # only include predictors that are explicitly included (default: all)
  c_include <- c_complete
  if(length(include)>0)
    c_include <- match_columns(df, include)

  c_include <- setdiff(setdiff(c_include, c_exclude), c_response[1]) 

  if(length(c_include) < 1)
  {
      c_include <- setdiff(setdiff(c_complete, c_exclude), c_response[1]) 
      
      if(length(c_include) < 1)
      {
        warning("At least one (not excluded) predictor is required; defaulting to all columns, excl. response.")
        c_include <- setdiff(c_complete, c_response)
        if(length(c_include) < 1)
          c_include <- setdiff(c_complete, c_response[1])
      }
      else
        warning("At least one (not excluded) predictor is required; defaulting to all columns, excl. response and exclusions")
  }

  # here we know the columns for Y, X and the exclusions (if needed)
  # Y <- as.factor(df[,c_response[1]])  # classification problem

  list(c_response[1], c_include, setdiff(setdiff(c_complete,c_include),c_response[1]))
}

#' Mode
#'
#' Compute the most frequently occurring element (mode) of a column or vector
#' @param v input vector
#' @keywords statistics
#' @export
#' @examples
#' compute_mode(c(1,2,3,3,3,3,4))
#' compute_mode(c("this", "and", "that", "or", "this"))
compute_mode <- function(v) 
{
   uniqv <- unique(v)
   uniqv[which.max(tabulate(match(v, uniqv)))]
}

#' Auto-detect suspicious columns
#'
#' Detect quasi-unique string (character) and factor columns with a threshold tl,
#' as well as columns that include too much repetition of one value with a threshold th
#' @param df data frame
#' @param th threshold for percentage of unique values
#' @param tl threshold for percentage of most frequent identical value
#' @keywords pre-processing
#' @export
#' @examples
#' df_uniqueness_properties(df)
df_uniqueness_properties <- function(df, th=0.9, tl=0.85) 
{
  c_class <- which(sapply(df,class) %in% c("character","factor"))

  # strings and factor columns with many unique values
  c_high <- as.numeric(which(sapply(df, function(x) length(unique(x))/length(x)) > th))
  c_high <- intersect(c_class, c_high)  # only for strings and factors

  # columns of any type with many repeated values (few unique values)
  c_low <- as.numeric(which(sapply(df, function(x) sum(x==compute_mode(x))/length(x)) > tl))

  list(c_high, c_low)
}

#' Auto-detect fake string columns that should be numeric
#'
#' Returns string columns that include mostly numbers with few exceptions (e.g. "N.A.")
#' @param df data frame
#' @param tn threshold for percentage of values that can successfully converted to numeric
#' @keywords pre-processing
#' @export
#' @examples
#' df_numeric_strings(df)
df_numeric_strings <- function(df, tn=0.6)
{
  c_char <-as.numeric(which(sapply(df,class)=="character"))   
  c_num <- as.numeric(which(sapply(df, function(x) 1-sum(is.na(suppressWarnings(as.numeric(x))))/length(x)) > tn))
  intersect(c_char, c_num)
}

#' Auto-detect non-factor columns that should be factors
#'
#' Returns columns for which there are very few unique values (should be factor)
#' @param df data frame
#' @param tf threshold for percentage of unique values under which it should be a factor
#' @keywords pre-processing
#' @export
#' @examples
#' df_missed_factors(df)
df_missed_factors <- function(df, tf=0.075)
{
  as.numeric(which(sapply(df, function(x) length(unique(x))/length(x)) < tf))
}


#' Convert cols to numeric
#'
#' Convert specific columns to numeric, coercing non-numerics to NA
#' @param df data frame
#' @param c_conv column indices of cols to be converted to numeric
#' @keywords pre-processing
#' @export
#' @examples
#' df_convert_cols_to_numeric(df, c(1,4,9))
df_convert_cols_to_numeric <- function(df, c_conv)
{
  for(i in c_conv)
    df[,i] <- suppressWarnings(as.numeric(df[,i]))
  df
}


#' Convert cols to factor
#'
#' Convert specific columns to factor
#' @param df data frame
#' @param c_conv column indices of cols to be converted to factor
#' @keywords pre-processing
#' @export
#' @examples
#' df_convert_cols_to_factor(df, c(1,4,9))
df_convert_cols_to_factor <- function(df, c_conv)
{
  for(i in c_conv)
    df[,i] <- suppressWarnings(as.factor(df[,i]))
  df
}

#' Remove cols 
#'
#' Exclude specific columns from data frame
#' @param df data frame
#' @param c_excl column indices of cols to be removed from data frame
#' @keywords pre-processing
#' @export
#' @examples
#' df_exclude_cols(df, c(1,4))
df_exclude_cols <- function(df, c_excl)
{
  for(i in c_excl)
    df[,i] <- NULL
  df
}

#' Get numeric cols
#'
#' Returns list of columns that are of numeric type
#' @param df data frame
#' @keywords pre-processing
#' @export
#' @examples
#' numeric_cols(df)
numeric_cols <- function(df)
{
  as.numeric(which(sapply(df, class)=="numeric"))
}

#' Split into train and test 1/2: training (in-sample) indices
#'
#' Create a random split of indices with proportion ptrain in the training set. Also see: generate_testing_row_indices(n, training_rows)
#' @param n numer of data points in total (training + testing)
#' @param ptrain percentage of rows in the training set (fraction)
#' @keywords splitting
#' @export
#' @examples
#' train <- generate_training_row_indices(nrow(df), 0.6)
#' test <- generate_testing_row_indices(nrow(df), train)
generate_training_row_indices <- function(n, ptrain)
{
  sample(n, ptrain*n)
}

#' Split into train and test 2/2: testing (out-of-sample) indices
#'
#' Create a random split of indices with proportion ptrain in the training set. Also see: generate_training_row_indices(n, training_rows) 
#' @param n numer of data points in total (training + testing)
#' @param ptrain percentage of rows in the training set (fraction)
#' @keywords splitting
#' @export
#' @examples
#' train <- generate_training_row_indices(nrow(df), 0.6)
#' test <- generate_testing_row_indices(nrow(df), train)
generate_testing_row_indices <- function(n, training_rows)
{
  sample((1:n)[!(1:n) %in% training_rows],n-length(training_rows))
}


#' Pre-processing of raw input file (BEFORE in-/out-of-sample split)
#'
#' This function takes a data frame and response column (target) and performs the following tasks:
#'   (1) For string (character) columns it checks if many of them are numbers and converts the column to numeric if necessary
#'   (2) For columns with very few unique values, it converts these columns to factors
#'   (3) For string (character) or factor columns it checks if most rows are unique and drops columns for which this is the case
#'   (4) For string (character) or factor columns it checks if one specific string or factor is very frequently present and drops those columns
#'
#' Additional pre-processing functions could be added here.
#'
#' It should be noted that these steps occur BEFORE the file is split into training and testing data and as such care must be taken to 
#' operations that lead to a bias when used to train a classifier (e.g. ranking of all values, imputation based on all values, etc). 
#'
#' Also see: post_sample_preprocessing() for a function that only gets called AFTER the file is split into training and testing data.
#'
#' The trade-off is speed vs bias. Putting everything here leads to a significant speed-up but with a potential bias and vice versa.
#'
#' @param df data frame
#' @param target column name or index of the response
#' @param tn_in threshold to determine if a column should be numeric
#' @param th_in threshold to determine if a string or factor column contains too many unique values and should be dropped
#' @param tl_in threshold to determine if a specific string or factor is repeated too many times and the column should be dropped 
#' @param tf_in threshold to determine if a column should be converted to a factor
#' @keywords pre-processing
#' @export
#' @examples
#' pre_sample_preprocessing <- function(df, 5)
pre_sample_preprocessing <- function(df, target, tn_in=0.6, th_in=0.9, tl_in=0.85, tf_in=0.075)
{
  # convert string columns that are mostly numeric to numeric (R imports cols with NAs as strings)
  d <- df_convert_cols_to_numeric(df, df_numeric_strings(df, tn=tn_in))
  d <- df_convert_cols_to_factor(d, df_missed_factors(d, tf=tf_in))

  # find text / factor columns with extreme uniqueness or extreme repetition
  u <- df_uniqueness_properties(d, th=th_in, tl=tl_in)

  # obtain indices for Y and X while excluding extreme columns as above
  r <- determine_relevant_columns(d, response=target, exclude=unique(unlist(u)), include=NULL)

  # X, Y 
  list(d[,r[[2]]], as.factor(d[,r[[1]]]))
}

#' Pre-processing of input file (AFTER in-/out-of-sample split)
#'
#' This function takes a data frame and performs the following tasks:
#'   (1) For each numeric column, it creates a ranking function based only on in-sample data
#'   (2) It applies this function to all numeric columns and to both in- and out-of-sample data (could also be applied to online data)
#'   (3) For each numeric column, it computes the median of only in-sample data 
#'   (4) It imputes missing values in numeric columns with these in-sample medians (could also be applied to online data)
#'
#' Additional pre-processing functions could be added here. For now, Ypre and r_test are not used.
#'
#' It should be noted that these steps occur AFTER the file is split into training and testing data. As long as only in-sample data 
#' is used to create the transformations, there is not bias when training a classifier. 
#'
#' Also see: pre_sample_preprocessing() for a function that already gets called BEFORE the file is split into training and testing data.
#'
#' The trade-off is speed vs bias. Putting everything here leads to slower run-times but without any potential for bias and vice versa.
#'
#' @param Xpre data frame
#' @param Ypre data frame, currently unused (but will make it easier later to add processing that depends on it)
#' @param r_train vector of in-sample indices into the data frames 
#' @param r_test vector of out-of-sample indices into the data frames, currently unused 
#' @keywords pre-processing
#' @export
#' @examples
#' r_train <- generate_training_row_indices(nrow(df), 0.6)
#' post_sample_preprocessing <- function(dfX, dfY, r_train, 0)
post_sample_preprocessing <- function(Xpre, Ypre, r_train, r_test)
{
  c_num <- numeric_cols(Xpre)

  # create ranking function from in-sample data ONLY and then apply to full column 
  X <- Xpre
  X[,c_num] <- sapply(Xpre[,c_num], function(x) uni_rank(x, zero_one_rank_map(x[r_train])))

  # impute all numeric missing values to median of in-sample data (not unknown overall median)
  X[,c_num] <- sapply(X[,c_num], function(x) impute(x, median(x[r_train])))
  X
}

#' Add noise dimensions
#'
#' Add n columns of white noise to data frame. Useful to compute the SNR of a classifier.
#' @param df data frame
#' @param n number of noise dimensions to add
#' @keywords noise
#' @export
#' @examples
#' add_noise_columns(df, 8)
add_noise_columns <- function(df, n)
{
  if(n>0)
    for(i in 1:n)
      	df[,paste("N",i,sep="")] <- rnorm(nrow(df),0,1)
  df
}

#' Create an array of random rotations (plus identity)
#'
#' Create an array containing the identity rotation plus numRR random rotations of (square) dimension dimRR
#' @param numRR number of random rotations to create (in addition to the identity rotation)
#' @param dimRR dimension of the generated random rotations
#' @keywords random rotation
#' @export
#' @examples
#' create_random_rotation_matrix_array(10,3)
create_random_rotation_matrix_array <- function(numRR, dimRR)
{
  mArray <- sapply(rep(dimRR,numRR),random_rotation_matrix,simplify=FALSE)
  c(list(diag(1,dimRR,dimRR)),mArray) # add identity rotation	
}	


#' Apply an array of random rotations to an input data frame
#'
#' Create an array containing a rotated version of the input data frame for each rotation in the matrix.
#' The dimension of the rotation matrices in mArray must match the number of numeric columns in df.
#' @param df data frame
#' @param mArray array of rotations
#' @keywords random rotation
#' @export
#' @examples
#' mArray <- create_random_rotation_matrix_array(5, dimRR=3)
#' df_apply_random_rotations(df, mArray)
df_apply_random_rotations <- function(df, mArray)
{
  c_num <- numeric_cols(df)
  if(sum(length(c_num)!=sapply(mArray,ncol)) && sum(length(c_num)!=sapply(mArray,nrow)))
    warning("Dimension of the rotation matrices in mArray must match the number of numeric columns in df")
    
  # multiply the matrix extracted from the numeric columns of df with each rotation matrix in the array
  m <- sapply(1:length(mArray),FUN=function(x){as.matrix(df[,c_num]) %*% mArray[[x]]},simplify=FALSE)
  
  # create the fully rotated data frames (incl. unrotated / non-numeric columns)
  tArray <- replicate(length(mArray), df, simplify=FALSE)
  for(i in 1:length(mArray))
    tArray[[i]][,c_num] <- m[[i]]
  
  tArray
}

#' Optimize tuning parameter
#'
#' First sort rotations by complexity (primary criterium) and then perform a grid search over h
#' to find the value that produces the rotation weighting with lowest OOB error.
#' @param v_complex vector of rotation complexities
#' @param v_oob vector of rotation OOB errors
#' @param fn weight function to use, currently in weights_{cut, exp, rre}
#' @param min_step minimal value of h to test (maximum is always R, the number of rotations)
#' @param step_size determines the granularity of the grid
#' @export
#' @examples
#' compute_tuning_parameter(v_complex=c(10,4,5,4,6), v_oob=c(0.05, 0.03, 0.1, 0.06, 0.02), fn=weights_exp, 0.1, 0.1)
compute_tuning_parameter <- function(v_complex, v_oob, fn=weights_cut, min_step=1, step_size=1)
{
  n <- length(v_oob)
  oob_error <- v_oob[order(v_complex)] # OOB error vector, sorted by tree complexity

  h_grid <- seq(min_step,n,step_size)
  w_curve <- sapply(h_grid,FUN=function(x) fn(n,x))
  w_error <- colSums(oob_error * w_curve)
  
  v_match <- which(w_error==min(w_error))
  h_grid[v_match[length(v_match)]]    # return first optimal h within grid
}
