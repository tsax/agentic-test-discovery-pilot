#!/usr/bin/env Rscript
# Run the real statcheck R package (Nuijten et al.) over paragraph texts.
#
# Usage:
#   Rscript run_statcheck.R <input> <output.csv>
#
# <input> is one of:
#   - a .tsv with header paper_id<TAB>para_id<TAB>text (tabs/newlines in text replaced
#     by spaces; quoting disabled) -- this is what statcheck_check.py writes;
#   - a .jsonl of {"paper_id","para_id","text"} (requires jsonlite);
#   - a directory of parsed-paper JSON files (corpus/parsed format; requires jsonlite).
#
# Output CSV columns: paper_id, para_id, raw, test_type, reported_p, computed_p,
#   error, decision_error, plus every other column statcheck returns (prefixed sc_).
# statcheck reports no character offsets; statcheck_check.py locates `raw` in the text.
# A CSV with only the header is written when statcheck finds nothing.

suppressPackageStartupMessages(library(statcheck))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) stop("usage: run_statcheck.R <input> <output.csv>")
inp <- args[1]; outp <- args[2]

read_input <- function(inp) {
  if (dir.exists(inp)) {
    suppressPackageStartupMessages(library(jsonlite))
    files <- list.files(inp, pattern = "\\.json$", full.names = TRUE)
    rows <- lapply(files, function(f) {
      p <- fromJSON(f, simplifyVector = FALSE)
      if (length(p$paragraphs) == 0) return(NULL)
      data.frame(paper_id = p$paper_id,
                 para_id = vapply(p$paragraphs, function(x) x$para_id, ""),
                 text = vapply(p$paragraphs, function(x) x$text, ""),
                 stringsAsFactors = FALSE)
    })
    return(do.call(rbind, rows))
  }
  if (grepl("\\.jsonl$", inp)) {
    suppressPackageStartupMessages(library(jsonlite))
    return(stream_in(file(inp), verbose = FALSE))
  }
  read.delim(inp, sep = "\t", quote = "", comment.char = "", header = TRUE,
             stringsAsFactors = FALSE, encoding = "UTF-8", na.strings = character(0),
             colClasses = "character")
}

df <- read_input(inp)
empty <- data.frame(paper_id = character(0), para_id = character(0), raw = character(0),
                    test_type = character(0), reported_p = character(0),
                    computed_p = character(0), error = character(0),
                    decision_error = character(0))
if (is.null(df) || nrow(df) == 0) { write.csv(empty, outp, row.names = FALSE); quit(status = 0) }

key <- paste(df$paper_id, df$para_id, sep = "\x1f")
texts <- enc2utf8(as.character(df$text))
names(texts) <- key

res <- suppressMessages(suppressWarnings(statcheck(texts, messages = FALSE)))
if (is.null(res) || nrow(res) == 0) { write.csv(empty, outp, row.names = FALSE); quit(status = 0) }

res <- as.data.frame(res, stringsAsFactors = FALSE)
nm <- tolower(names(res))
pick <- function(cands) { i <- match(cands, nm); i <- i[!is.na(i)]; if (length(i)) res[[i[1]]] else NA }
src <- as.character(pick(c("source")))
parts <- strsplit(src, "\x1f", fixed = TRUE)
out <- data.frame(
  paper_id = vapply(parts, `[`, "", 1),
  para_id = vapply(parts, `[`, "", 2),
  raw = as.character(pick(c("raw"))),
  test_type = as.character(pick(c("test_type", "statistic"))),
  reported_p = as.character(pick(c("reported_p", "reported.p.value"))),
  computed_p = as.character(pick(c("computed_p", "computed"))),
  error = as.character(pick(c("error"))),
  decision_error = as.character(pick(c("decision_error", "decisionerror"))),
  stringsAsFactors = FALSE)
extra <- res[, setdiff(names(res), c("source", "Source")), drop = FALSE]
names(extra) <- paste0("sc_", names(extra))
write.csv(cbind(out, extra), outp, row.names = FALSE, fileEncoding = "UTF-8")
