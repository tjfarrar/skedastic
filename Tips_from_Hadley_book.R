# Don't use for (i in 1:length(somevector)) ; rather use for (i in seq_along(somevector))

# Avoid is.numeric
# Lists are generally more time-efficient than data frames for element-wise computation, especially for large (>800) columns/elements


#Individual attributes can be retrieved and modified with attr(), or retrieved en masse with attributes(), and set en masse with structure().

#NROW, NCOL are more general than nrow, ncol

# unclass

# c() will combine several lists into one. If given a combination of atomic vectors and lists, c() will coerce the vectors to lists before combining them.

# tibble instead of data frame
# as.matrix vs. data.matrix

# View()

# Subsetting: Negative integers exclude elements at the specified positions (use RATHER THAN setdiff)

# drop=F preserves dimensionality when subsetting. Use it!

# While you must use [[ when working with lists, I’d also recommend using it with atomic vectors whenever you want to extract a single value... Doing so reinforces the expectation that you are getting and setting individual values.

#  I highly recommend setting the global option warnPartialMatchDollar to TRUE:

#options(warnPartialMatchDollar = TRUE)

# Tip: The broom-package provides a very useful approach to work with models in a tidy way

# With lists, you can use x[[i]] <- NULL to remove a component. To add a literal NULL, use x[i] <- list(NULL)

# Subsetting with nothing  e.g. x[]<- can be useful with assignment because it preserves the structure of the original object.
# subset(), merge(), dplyr::arrange()
# merge is for table joins
# arrange is for sorting a data frame or tibble
# match does same thing as %in%
# switch

# Sys.setenv("_R_CHECK_LENGTH_1_CONDITION_" = "true")
# Gives error if if statement has condition of length > 1
#dplyr::case_when() as more general version of ifelse

#srcref gives source code for a function

# %>% "and then" piping

# codetools::findGlobals() lists external dependencies in a function
# formals() gives arguments of a function; body() gives code; environment() gives environment

# include expensive computations in function arguments as it will only be evaluated if argument is used

# range

# use ... in functions

# Functionals: instead of for loops use purrr:map and map_dbl, map_int etc. These are preferable to lapply, vapply

