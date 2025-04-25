Taylors_Law_Chapter_04_Timecourse
================
Soren Maret
2025-04-07

``` r
Mouse_counts <- read_tsv("/scratch/Shares/rinnclass/MASTER_CLASS/STUDENTS/genehomies/DATA/MOUSE/salmon.merged.gene_counts.tsv")
```

    ## Rows: 55401 Columns: 17
    ## ── Column specification ────────────────────────────────────────────────────────
    ## Delimiter: "\t"
    ## chr  (2): gene_id, gene_name
    ## dbl (15): WT_0_1, WT_0_2, WT_0_3, WT_12_1, WT_12_2, WT_12_3, WT_24_1, WT_24_...
    ## 
    ## ℹ Use `spec()` to retrieve the full column specification for this data.
    ## ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.

``` r
counts_0 <- Mouse_counts[, grep("^WT_0_", colnames(Mouse_counts))]
counts_12 <- Mouse_counts[, grep("^WT_12_", colnames(Mouse_counts))]
counts_24 <- Mouse_counts[, grep("^WT_24_", colnames(Mouse_counts))]
counts_48 <- Mouse_counts[, grep("^WT_48_", colnames(Mouse_counts))]
counts_96 <- Mouse_counts[, grep("^WT_96_", colnames(Mouse_counts))]
```

``` r
g2s <- Mouse_counts[, c(1,2)]

Mouse_counts <- Mouse_counts[, -c(1,2)]
```

``` r
taylor_func <- function(df) {
  
  df_numeric <- df[sapply(df, is.numeric)]
  
  row_mean <- rowMeans(df_numeric, na.rm = TRUE)
  row_sd <- apply(df_numeric, 1, sd, na.rm = TRUE)
  
  ln_mean <- log(row_mean)
  ln_sd <- log(row_sd)
  
  # Combine results into dataframe
  result_df <- data.frame(
    row_mean = row_mean,
    row_sd = row_sd,
    log_m = ln_mean,
    log_sd = ln_sd
  )
  
  # Remove rows where row_mean is 0
  result_df <- result_df[result_df$row_mean !=0, ]
  
  return(result_df)
}
```

``` r
df_names <- ls()[sapply(ls(), function(x) is.data.frame(get(x)))]

regression_results <- tibble(
  Dataset = character(),            
  Slope = numeric(),                # regression slope (a in y = ax + b)
  Intercept = numeric(),            # regression intercept (b in y = ax + b)
  R_squared = numeric(),            # R² value 
  P_value = numeric(),              # p-value from the regression slope
  T_value_vs_first = numeric(),     # t-stat comparing slope to that of first model
  P_value_vs_first = numeric(),     # p-value of the t-test comparing slopes
  time = numeric()
)


# Initialize variables to store the first model's parameters
first_model <- NULL  # used to store the first model object
b1 <- NULL           # slope from the first model
se1 <- NULL          # standard error of slope from first model
df1 <- NULL          # degrees of freedom of first model

# Loop through each dataframe name in df_names
for (i in seq_along(df_names)) {
  df_name <- df_names[i]          # current dataframe name as string
  df <- get(df_name)              # retrieve dataframe object from the environment
  df <- taylor_func(df)           # apply your preprocessing/statistics function

  # Check if the expected columns exist before running regression
  if (all(c("log_m", "log_sd") %in% colnames(df))) {
    # Filter out incomplete or infinite values
    valid <- complete.cases(df$log_m, df$log_sd) &
             !is.infinite(df$log_m) &
             !is.infinite(df$log_sd)

    # Continue only if enough valid data points exist
    if (sum(valid) >= 2) {
      # Perform linear regression: log_sd ~ log_m
      model <- lm(log_sd ~ log_m, data = df[valid, ])

      # Extract regression slope (a), intercept (b), R², and p-value
      a <- coef(model)[["log_m"]]
      b <- coef(model)[["(Intercept)"]]
      r2 <- summary(model)$r.squared
      p_val <- coef(summary(model))["log_m", "Pr(>|t|)"]

      # If this is the first model, store its values for later comparison
      if (is.null(first_model)) {
        first_model <- model
        b1 <- a
        se1 <- coef(summary(model))["log_m", "Std. Error"]
        df1 <- df.residual(model)

        # No comparison yet for first model
        t_val <- NA
        p_val_vs_first <- NA
      } else {
        # For subsequent models, extract slope and SE
        b2 <- a
        se2 <- coef(summary(model))["log_m", "Std. Error"]
        df2 <- df.residual(model)

        # Calculate t-statistic to compare slope with first model's slope
        t_val <- (b1 - b2) / sqrt(se1^2 + se2^2)
        df_comp <- min(df1, df2)  # conservative degrees of freedom
        p_val_vs_first <- 2 * pt(-abs(t_val), df = df_comp)  # two-tailed p-value
      }

      # Add results to the output tibble
      regression_results <- add_row(
        regression_results,
        Dataset = df_name,
        Slope = a,
        Intercept = b,
        R_squared = r2,
        P_value = p_val,
        T_value_vs_first = t_val,
        P_value_vs_first = p_val_vs_first,
        time = i
      )
    } else {
      # Warn if there aren't enough points to perform regression
      warning(paste("Skipping", df_name, "- not enough valid points."))
    }

  } else {
    # Warn if expected columns are missing
    warning(paste("Skipping", df_name, "- missing log_m or log_sd columns."))
  }
}
```

    ## Warning: Skipping g2s - not enough valid points.

``` r
# Print out the full regression summary table
regression_results <- regression_results %>% filter(row_number() <= n()-1)

print(regression_results)
```

    ## # A tibble: 5 × 8
    ##   Dataset   Slope Intercept R_squared P_value T_value_vs_first P_value_vs_first
    ##   <chr>     <dbl>     <dbl>     <dbl>   <dbl>            <dbl>            <dbl>
    ## 1 counts_0  0.717    -0.240     0.938       0            NA          NA        
    ## 2 counts_12 0.812    -0.195     0.958       0           -63.3         0        
    ## 3 counts_24 0.668    -0.172     0.921       0            29.9         4.19e-193
    ## 4 counts_48 0.710    -0.220     0.932       0             4.05        5.11e-  5
    ## 5 counts_96 0.785    -0.205     0.953       0           -44.6         0        
    ## # ℹ 1 more variable: time <dbl>

``` r
ggplot(regression_results, mapping = aes(x = time, y = Slope)) + 
  geom_point() +
  geom_smooth(formula = 'y ~ x', method = "loess")
```

    ## Warning in simpleLoess(y, x, w, span, degree = degree, parametric = parametric,
    ## : span too small.  fewer data values than degrees of freedom.

    ## Warning in simpleLoess(y, x, w, span, degree = degree, parametric = parametric,
    ## : pseudoinverse used at 0.98

    ## Warning in simpleLoess(y, x, w, span, degree = degree, parametric = parametric,
    ## : neighborhood radius 2.02

    ## Warning in simpleLoess(y, x, w, span, degree = degree, parametric = parametric,
    ## : reciprocal condition number 0

    ## Warning in simpleLoess(y, x, w, span, degree = degree, parametric = parametric,
    ## : There are other near singularities as well. 4.0804

    ## Warning in predLoess(object$y, object$x, newx = if (is.null(newdata)) object$x
    ## else if (is.data.frame(newdata))
    ## as.matrix(model.frame(delete.response(terms(object)), : span too small.  fewer
    ## data values than degrees of freedom.

    ## Warning in predLoess(object$y, object$x, newx = if (is.null(newdata)) object$x
    ## else if (is.data.frame(newdata))
    ## as.matrix(model.frame(delete.response(terms(object)), : pseudoinverse used at
    ## 0.98

    ## Warning in predLoess(object$y, object$x, newx = if (is.null(newdata)) object$x
    ## else if (is.data.frame(newdata))
    ## as.matrix(model.frame(delete.response(terms(object)), : neighborhood radius
    ## 2.02

    ## Warning in predLoess(object$y, object$x, newx = if (is.null(newdata)) object$x
    ## else if (is.data.frame(newdata))
    ## as.matrix(model.frame(delete.response(terms(object)), : reciprocal condition
    ## number 0

    ## Warning in predLoess(object$y, object$x, newx = if (is.null(newdata)) object$x
    ## else if (is.data.frame(newdata))
    ## as.matrix(model.frame(delete.response(terms(object)), : There are other near
    ## singularities as well. 4.0804

    ## Warning in max(ids, na.rm = TRUE): no non-missing arguments to max; returning
    ## -Inf

![](Taylors_Law_Chapter_04_files/figure-gfm/unnamed-chunk-6-1.png)<!-- -->

``` r
ggplot(regression_results, mapping = aes(x = time, y = Intercept)) + 
  geom_point() +
  geom_smooth(formula = 'y ~ x', method = "loess")
```

    ## Warning in simpleLoess(y, x, w, span, degree = degree, parametric = parametric,
    ## : span too small.  fewer data values than degrees of freedom.

    ## Warning in simpleLoess(y, x, w, span, degree = degree, parametric = parametric,
    ## : pseudoinverse used at 0.98

    ## Warning in simpleLoess(y, x, w, span, degree = degree, parametric = parametric,
    ## : neighborhood radius 2.02

    ## Warning in simpleLoess(y, x, w, span, degree = degree, parametric = parametric,
    ## : reciprocal condition number 0

    ## Warning in simpleLoess(y, x, w, span, degree = degree, parametric = parametric,
    ## : There are other near singularities as well. 4.0804

    ## Warning in predLoess(object$y, object$x, newx = if (is.null(newdata)) object$x
    ## else if (is.data.frame(newdata))
    ## as.matrix(model.frame(delete.response(terms(object)), : span too small.  fewer
    ## data values than degrees of freedom.

    ## Warning in predLoess(object$y, object$x, newx = if (is.null(newdata)) object$x
    ## else if (is.data.frame(newdata))
    ## as.matrix(model.frame(delete.response(terms(object)), : pseudoinverse used at
    ## 0.98

    ## Warning in predLoess(object$y, object$x, newx = if (is.null(newdata)) object$x
    ## else if (is.data.frame(newdata))
    ## as.matrix(model.frame(delete.response(terms(object)), : neighborhood radius
    ## 2.02

    ## Warning in predLoess(object$y, object$x, newx = if (is.null(newdata)) object$x
    ## else if (is.data.frame(newdata))
    ## as.matrix(model.frame(delete.response(terms(object)), : reciprocal condition
    ## number 0

    ## Warning in predLoess(object$y, object$x, newx = if (is.null(newdata)) object$x
    ## else if (is.data.frame(newdata))
    ## as.matrix(model.frame(delete.response(terms(object)), : There are other near
    ## singularities as well. 4.0804

    ## Warning in max(ids, na.rm = TRUE): no non-missing arguments to max; returning
    ## -Inf

![](Taylors_Law_Chapter_04_files/figure-gfm/unnamed-chunk-6-2.png)<!-- -->

``` r
ggplot(regression_results, mapping = aes(x = time, y = R_squared)) + 
  geom_point() +
  geom_smooth(formula = 'y ~ x', method = "loess")
```

    ## Warning in simpleLoess(y, x, w, span, degree = degree, parametric = parametric,
    ## : span too small.  fewer data values than degrees of freedom.

    ## Warning in simpleLoess(y, x, w, span, degree = degree, parametric = parametric,
    ## : pseudoinverse used at 0.98

    ## Warning in simpleLoess(y, x, w, span, degree = degree, parametric = parametric,
    ## : neighborhood radius 2.02

    ## Warning in simpleLoess(y, x, w, span, degree = degree, parametric = parametric,
    ## : reciprocal condition number 0

    ## Warning in simpleLoess(y, x, w, span, degree = degree, parametric = parametric,
    ## : There are other near singularities as well. 4.0804

    ## Warning in predLoess(object$y, object$x, newx = if (is.null(newdata)) object$x
    ## else if (is.data.frame(newdata))
    ## as.matrix(model.frame(delete.response(terms(object)), : span too small.  fewer
    ## data values than degrees of freedom.

    ## Warning in predLoess(object$y, object$x, newx = if (is.null(newdata)) object$x
    ## else if (is.data.frame(newdata))
    ## as.matrix(model.frame(delete.response(terms(object)), : pseudoinverse used at
    ## 0.98

    ## Warning in predLoess(object$y, object$x, newx = if (is.null(newdata)) object$x
    ## else if (is.data.frame(newdata))
    ## as.matrix(model.frame(delete.response(terms(object)), : neighborhood radius
    ## 2.02

    ## Warning in predLoess(object$y, object$x, newx = if (is.null(newdata)) object$x
    ## else if (is.data.frame(newdata))
    ## as.matrix(model.frame(delete.response(terms(object)), : reciprocal condition
    ## number 0

    ## Warning in predLoess(object$y, object$x, newx = if (is.null(newdata)) object$x
    ## else if (is.data.frame(newdata))
    ## as.matrix(model.frame(delete.response(terms(object)), : There are other near
    ## singularities as well. 4.0804

    ## Warning in max(ids, na.rm = TRUE): no non-missing arguments to max; returning
    ## -Inf

![](Taylors_Law_Chapter_04_files/figure-gfm/unnamed-chunk-6-3.png)<!-- -->

``` r
ggplot(regression_results, mapping = aes(x = time, y = P_value_vs_first)) + 
  geom_point() +
  geom_smooth(formula = 'y ~ x', method = "loess")
```

    ## Warning: Removed 1 row containing non-finite outside the scale range
    ## (`stat_smooth()`).

    ## Warning in simpleLoess(y, x, w, span, degree = degree, parametric = parametric,
    ## : span too small.  fewer data values than degrees of freedom.

    ## Warning in simpleLoess(y, x, w, span, degree = degree, parametric = parametric,
    ## : pseudoinverse used at 1.985

    ## Warning in simpleLoess(y, x, w, span, degree = degree, parametric = parametric,
    ## : neighborhood radius 2.015

    ## Warning in simpleLoess(y, x, w, span, degree = degree, parametric = parametric,
    ## : reciprocal condition number 0

    ## Warning in simpleLoess(y, x, w, span, degree = degree, parametric = parametric,
    ## : There are other near singularities as well. 4.0602

    ## Warning in predLoess(object$y, object$x, newx = if (is.null(newdata)) object$x
    ## else if (is.data.frame(newdata))
    ## as.matrix(model.frame(delete.response(terms(object)), : span too small.  fewer
    ## data values than degrees of freedom.

    ## Warning in predLoess(object$y, object$x, newx = if (is.null(newdata)) object$x
    ## else if (is.data.frame(newdata))
    ## as.matrix(model.frame(delete.response(terms(object)), : pseudoinverse used at
    ## 1.985

    ## Warning in predLoess(object$y, object$x, newx = if (is.null(newdata)) object$x
    ## else if (is.data.frame(newdata))
    ## as.matrix(model.frame(delete.response(terms(object)), : neighborhood radius
    ## 2.015

    ## Warning in predLoess(object$y, object$x, newx = if (is.null(newdata)) object$x
    ## else if (is.data.frame(newdata))
    ## as.matrix(model.frame(delete.response(terms(object)), : reciprocal condition
    ## number 0

    ## Warning in predLoess(object$y, object$x, newx = if (is.null(newdata)) object$x
    ## else if (is.data.frame(newdata))
    ## as.matrix(model.frame(delete.response(terms(object)), : There are other near
    ## singularities as well. 4.0602

    ## Warning: Removed 1 row containing missing values or values outside the scale range
    ## (`geom_point()`).

    ## Warning in max(ids, na.rm = TRUE): no non-missing arguments to max; returning
    ## -Inf

![](Taylors_Law_Chapter_04_files/figure-gfm/unnamed-chunk-6-4.png)<!-- -->

``` r
ggplot(regression_results, mapping = aes(x = time, y = T_value_vs_first)) + 
  geom_point() +
  geom_smooth(formula = 'y ~ x', method = "loess")
```

    ## Warning: Removed 1 row containing non-finite outside the scale range
    ## (`stat_smooth()`).

    ## Warning in simpleLoess(y, x, w, span, degree = degree, parametric = parametric,
    ## : span too small.  fewer data values than degrees of freedom.

    ## Warning in simpleLoess(y, x, w, span, degree = degree, parametric = parametric,
    ## : pseudoinverse used at 1.985

    ## Warning in simpleLoess(y, x, w, span, degree = degree, parametric = parametric,
    ## : neighborhood radius 2.015

    ## Warning in simpleLoess(y, x, w, span, degree = degree, parametric = parametric,
    ## : reciprocal condition number 0

    ## Warning in simpleLoess(y, x, w, span, degree = degree, parametric = parametric,
    ## : There are other near singularities as well. 4.0602

    ## Warning in predLoess(object$y, object$x, newx = if (is.null(newdata)) object$x
    ## else if (is.data.frame(newdata))
    ## as.matrix(model.frame(delete.response(terms(object)), : span too small.  fewer
    ## data values than degrees of freedom.

    ## Warning in predLoess(object$y, object$x, newx = if (is.null(newdata)) object$x
    ## else if (is.data.frame(newdata))
    ## as.matrix(model.frame(delete.response(terms(object)), : pseudoinverse used at
    ## 1.985

    ## Warning in predLoess(object$y, object$x, newx = if (is.null(newdata)) object$x
    ## else if (is.data.frame(newdata))
    ## as.matrix(model.frame(delete.response(terms(object)), : neighborhood radius
    ## 2.015

    ## Warning in predLoess(object$y, object$x, newx = if (is.null(newdata)) object$x
    ## else if (is.data.frame(newdata))
    ## as.matrix(model.frame(delete.response(terms(object)), : reciprocal condition
    ## number 0

    ## Warning in predLoess(object$y, object$x, newx = if (is.null(newdata)) object$x
    ## else if (is.data.frame(newdata))
    ## as.matrix(model.frame(delete.response(terms(object)), : There are other near
    ## singularities as well. 4.0602

    ## Warning: Removed 1 row containing missing values or values outside the scale range
    ## (`geom_point()`).

    ## Warning in max(ids, na.rm = TRUE): no non-missing arguments to max; returning
    ## -Inf

![](Taylors_Law_Chapter_04_files/figure-gfm/unnamed-chunk-6-5.png)<!-- -->

``` r
ggplot(regression_results, mapping = aes(x = R_squared, y = Slope)) +
  geom_point() +
  geom_smooth(formula = 'y ~ x', method = "lm")
```

![](Taylors_Law_Chapter_04_files/figure-gfm/unnamed-chunk-6-6.png)<!-- -->

``` r
  # Loop over all objects with names that contain "counts_"
for (obj_name in ls(pattern = "^counts_")) {
  
  # Get the object
  df <- get(obj_name)
  
  # Check it's a data frame
  if (is.data.frame(df)) {
    
    # Identify only the numeric columns (e.g., expression counts)
    numeric_cols <- sapply(df, is.numeric)
    
    # Compute stats only on numeric columns
    row_mean <- rowMeans(df[, numeric_cols], na.rm = TRUE)
    row_sd   <- apply(df[, numeric_cols], 1, sd, na.rm = TRUE)
    row_se   <- row_sd / sqrt(rowSums(!is.na(df[, numeric_cols])))
    
    # Append new columns to the original dataframe
    df$row_mean <- row_mean
    df$row_sd   <- row_sd
    df$row_se   <- row_se

    # Assign the modified dataframe back to the original name
    assign(obj_name, df)
  }
}
```

``` r
save(g2s, regression_results, taylor_func, counts_0, counts_12, counts_24, counts_48, counts_96, file = "/scratch/Shares/rinnclass/MASTER_CLASS/STUDENTS/genehomies/RESULTS/TPL_ANALYSIS/RData/chapter_04.RData")
```
