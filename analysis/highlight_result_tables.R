# Highlight results at the displayed precision, preserving all numerical entries.
highlight_result_tables <- function(path, highlight_coverage = FALSE,
                                    highlight_variance_mse = FALSE) {
  lines <- readLines(path, warn = FALSE)
  starts <- grep("\\begin{tabular}", lines, fixed = TRUE)
  ends <- grep("\\end{tabular}", lines, fixed = TRUE)
  unbold <- function(x) gsub("\\\\textbf\\{([^{}]*)\\}", "\\1", x, perl = TRUE)
  numeric_cell <- function(x) {
    x <- unbold(x)
    x <- gsub("\\\\rev\\{([^{}]*)\\}", "\\1", x, perl = TRUE)
    x <- sub("\\\\\\\\.*$", "", x)
    suppressWarnings(as.numeric(trimws(x)))
  }
  mark_cell <- function(x, winner, blue_table) {
    clean <- unbold(x)
    if (!winner) return(clean)
    match <- regexpr("[+-]?[0-9]+(?:\\.[0-9]+)?(?:[eE][+-]?[0-9]+)?", clean, perl = TRUE)
    stopifnot(match[1] > 0L)
    value <- regmatches(clean, match)
    bold <- paste0("\\textbf{", value, "}")
    if (!blue_table && !grepl("\\rev{", clean, fixed = TRUE) &&
        !grepl("\\textbf{", x, fixed = TRUE)) bold <- paste0("\\rev{", bold, "}")
    regmatches(clean, match) <- bold
    clean
  }
  tables <- cells_marked <- 0L
  previous_end <- 0L
  for (start in starts) {
    end <- min(ends[ends > start])
    block <- seq.int(start + 1L, end - 1L)
    headers <- block[grepl(" & ", lines[block], fixed = TRUE)]
    if (!length(headers)) { previous_end <- end; next }
    header <- headers[1]
    names <- trimws(strsplit(lines[header], " & ", fixed = TRUE)[[1]])
    bias <- grepl("Bias(", names, fixed = TRUE)
    rmse <- grepl("RMSE", names, fixed = TRUE)
    sd <- startsWith(names, "SD(")
    variance_mse <- highlight_variance_mse & names == "Var\\_MSE"
    coverage <- highlight_coverage & startsWith(names, "CP(")
    selected <- which(bias | rmse | sd | variance_mse | coverage)
    if (!length(selected)) { previous_end <- end; next }
    rows <- list()
    groups <- integer()
    group <- 0L
    for (index in seq.int(header + 1L, end - 1L)) {
      if (grepl("\\\\(c?midrule|multicolumn)", lines[index], perl = TRUE)) group <- group + 1L
      row <- strsplit(lines[index], " & ", fixed = TRUE)[[1]]
      if (length(row) != length(names) ||
          !any(is.finite(vapply(row[selected], numeric_cell, numeric(1))))) next
      rows[[as.character(index)]] <- row
      groups <- c(groups, group)
    }
    stopifnot(length(rows) > 0L)
    prefix <- seq.int(previous_end + 1L, start - 1L)
    blue_table <- any(grepl("\\ifshowrevisions\\color{blue}\\fi", lines[prefix], fixed = TRUE))
    for (g in unique(groups)) {
      members <- which(groups == g)
      for (column in selected) {
        values <- vapply(rows[members], function(row) numeric_cell(row[column]), numeric(1))
        metric <- if (bias[column]) abs(values) else values
        if (coverage[column]) metric <- round(abs(values - .95), 12L)
        stopifnot(any(is.finite(metric)))
        best <- min(metric[is.finite(metric)])
        for (i in seq_along(members)) {
          member <- members[i]
          winner <- is.finite(metric[i]) && metric[i] == best
          rows[[member]][column] <- mark_cell(rows[[member]][column], winner, blue_table)
          cells_marked <- cells_marked + as.integer(winner)
        }
      }
    }
    for (index in names(rows)) lines[as.integer(index)] <- paste(rows[[index]], collapse = " & ")
    captions <- prefix[grepl("\\caption{", lines[prefix], fixed = TRUE)]
    stopifnot(length(captions) > 0L)
    caption <- tail(captions, 1L)
    minima <- c(if (any(bias)) "absolute bias", if (any(rmse)) "RMSE",
                if (any(sd)) "SD", if (any(variance_mse)) "variance MSE")
    description <- function(items) {
      if (length(items) == 1L) return(items)
      paste(paste(head(items, -1L), collapse = ", "), tail(items, 1L), sep = " and ")
    }
    note <- if (any(coverage) && any(variance_mse)) {
      paste0("Within each condition, bold marks the smallest ", description(minima),
             ", as well as the CP closest to 0.95. Ties at the displayed precision are included.")
    } else if (any(coverage)) {
      paste0("Within each condition, bold marks the smallest ", description(minima),
             " and the CP closest to 0.95, including ties at the displayed precision.")
    } else {
      paste0("Bold entries indicate the smallest ", description(minima),
             " within each condition, including ties at the displayed precision.")
    }
    if (!grepl(note, lines[caption], fixed = TRUE)) {
      legacy <- " Bold entries indicate the smallest absolute bias and the smallest RMSE among the main-text methods for each scenario and parameter."
      lines[caption] <- gsub(legacy, "", lines[caption], fixed = TRUE)
      revised_caption <- grepl("\\caption{\\rev{", lines[caption], fixed = TRUE)
      braces <- if (revised_caption) 2L else 1L
      text <- if (revised_caption) note else paste0("\\rev{", note, "}")
      position <- nchar(lines[caption]) - braces
      lines[caption] <- paste0(substr(lines[caption], 1L, position), " ", text,
                               substr(lines[caption], position + 1L, nchar(lines[caption])))
    }
    tables <- tables + 1L
    previous_end <- end
  }
  writeLines(lines, path)
  invisible(c(tables = tables, highlighted_cells = cells_marked))
}
