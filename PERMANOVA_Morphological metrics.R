# Load packages
need <- c(
  "readr",
  "dplyr",
  "tibble",
  "tidyr",
  "vegan",
  "ggplot2",
  "ggrepel",
  "ragg"
)

to_install <- setdiff(need, rownames(installed.packages()))

if (length(to_install)) {
  install.packages(to_install, dependencies = TRUE)
}

invisible(lapply(need, library, character.only = TRUE))

# Read data
df <- readr::read_csv(
  "Morphological metrics.csv",
  show_col_types = FALSE
) |>
  
  dplyr::mutate(
    species = factor(
      species,
      levels = c("S. acapulcoensis", "S. flavilatus")
    )
  )

vars <- c("2", "3", "4", "5", "6", "7")

stopifnot(all(vars %in% names(df)))

# Remove rows with NA values in morphological variables
df <- df |>
  tidyr::drop_na(dplyr::all_of(vars))

X <- as.data.frame(
  dplyr::select(df, dplyr::all_of(vars))
)

# PERMANOVA using Euclidean distance
set.seed(123)

perm <- adonis2(
  X ~ species,
  data = df,
  method = "euclidean",
  permutations = 9999
)

print(perm)

