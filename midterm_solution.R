library(here)
library(tidyverse)

###############################
# Reading contaminants data
###############################

## Option 1
readData <- function(path){
    files = paste0(path, "/", dir(path))
    df <- read.csv(files[1])
    for(f in 2:length(files)){
        df <- rbind(df, read.csv(files[f]))
    }
    df
}

## Option 2
readData <- function(path){
    files = dir(path, full.names = TRUE)
    list_of_csv <- vector("list", length(files))
    for(f in seq_along(files))
        list_of_csv[[f]] <- read.csv(files[f])
    df <- do.call(rbind, list_of_csv)
}

## Option 3 - my choice
readData <- function(path){
    files = dir(path, full.names = TRUE)
    df <- map(files, read.csv) |> bind_rows()
    df
}

## Option 4 - overly smart maybe?
readData <- function(path){
    read_csv(dir(path, full.names = TRUE))
}

path = here("contaminants")
contaminants = readData(path)

###############################
# Reading species data
###############################

species_data = read_csv("species_counts.csv")

###############################
# Creating Month column and summarizing
###############################

## Option 1
contaminants |>
    separate(Date, into = c("Year", "Month", "Date")) |>
    mutate(YearMonth = paste(Year, Month, sep = "-"))  |>
    group_by(YearMonth, Location) |>
    na.omit() |> # probably a bad idea, could remove useful rows
    summarize(avgPb = mean(Lead), 
              avgHg = mean(Mercury))


## Option 2
contaminants_summary <- contaminants |>
    separate_wider_regex(Date, c(Year = "\\d{4}", "-", 
                                 Month = "\\d{2}", "-\\d{2}") ) |>       
    mutate(YearMonth = paste(Year, Month, sep = "-")) |>
    group_by(YearMonth, Location) |>
    summarize(avgPb = mean(Lead, na.rm = T), 
              avgHg = mean(Mercury, na.rm = T))

###############################
# Creating diversity fuctions
###############################

getProb = function(x) x/sum(x)
zero.omit = function(x) x[x!=0]
Simpson <- function(x){
  x <- as.numeric(x)
  x <- na.omit(x)
  p = getProb(x)
  D = sum(p^2)
  return(c(Simpson = D))
}
Shannon <- function(x){
  x <- as.numeric(x)
  x <- zero.omit(na.omit(x))
  p = getProb(x)
  H = -sum(p*log(p))
  return(c(Shannon = H))
}

###############################
# Calculating diversity indexes
###############################

# Option 1 - long data, a bit cleaner

diversity_indexes = species_data |> 
    pivot_longer(Species_1:Species_10, 
                 names_to = "species", 
                 values_to = "counts") |>
    group_by(Month, Location) |> 
    summarize(Simpson = Simpson(counts), 
              Shannon = Shannon(counts))

# Option 2 - wide data

species_data |> 
    mutate(Shannon = apply(species_data[,-c(1, 2)], 1, Shannon),
           Simpson = apply(species_data[,-c(1, 2)], 1, Simpson)) |>
    select(Month, Location, Simpson, Shannon)

###############################
# Extra point: creating and using the Diversity function
###############################

## Easy version for wide data

Diversity = function(x, index){
    apply(x[,-c(1:2)], 1, index)
}

species_data |> 
    mutate(Shannon = Diversity(species_data, Shannon),
           Simpson = Diversity(species_data, Simpson)) |>
    select(Month, Location, Simpson, Shannon)


## The brutally complicated version that works well with the long data

Diversity = function(x, col, ...){
    index = rlang::list2(...)                     # Getting a list of index functions
    index_name <- match.call(expand.dots = FALSE) # Getting a list of index function names
    index_name <- as.character(index_name$...)    # Setting the function names to character
    map2(index, index_name,                       # Mapping over functions and function names
      function(f, name){ x |> 
        summarize(current_index = f({{col}}))  |> # Calculating the indexes on the col argument
        rename_at("current_index", ~ name)}) |>   # Setting the summary column name
        reduce(inner_join)                        # Joining all indexes
}

species_data |> 
    pivot_longer(Species_1:Species_10, 
                 names_to = "species", 
                 values_to = "counts") |>
    group_by(Month, Location) |>
    Diversity(counts, Simpson, Shannon)

# We can even add more diversity functions without changing anything!

Gini <- function(x) {
  x <- as.numeric(x)
  x <- zero.omit(na.omit(x))
  n <- length(x)
  gini <- 1 + (1 / n) - 2 * sum(cumsum(sort(x)) / sum(x)) / n
  return(gini)
}

species_data |> 
    pivot_longer(Species_1:Species_10, 
                 names_to = "species", 
                 values_to = "counts") |>
    group_by(Month, Location) |>
    Diversity(counts, Simpson, Shannon, Gini)

#########################################
# Merging contaminants and diversity data
#########################################

diversity_contaminants = inner_join(diversity_indexes,
                                    contaminants_summary, 
                                    by = c(Month = "YearMonth", "Location")) 

########################################
# Plots
########################################

# Temporal trends

## By day
contaminants  |> 
    pivot_longer(Lead:Mercury,
                 names_to = "contaminants", 
                 values_to  = "concentration", 
                 values_drop_na = TRUE) |>
    ggplot(aes(Date, concentration, color = contaminants, 
               group = contaminants)) +
        facet_wrap(~Location) + 
        geom_point() + 
        theme_bw()

## By month
contaminants_summary |> 
    pivot_longer(avgPb:avgHg,
                 names_to = "contaminants", 
                 values_to  = "concentration", 
                 values_drop_na = TRUE) |>
    ggplot(aes(YearMonth, concentration, color = contaminants, 
               group = contaminants)) +
        facet_wrap(~Location) + 
        geom_point() + 
        theme_bw()

# Relation between diversity and contaminants

library(patchwork)
diversity_contaminants |> 
    pivot_longer(avgPb:avgHg, 
                 names_to = "contaminants", 
                 values_to  = "concentration") |>
    pivot_longer(Shannon:Simpson, 
                 names_to = "index", 
                 values_to  = "diversity") |>
    ggplot(aes(concentration, diversity, color = index)) + 
        geom_vline(xintercept = 300) +           
        geom_point() + 
        facet_grid(contaminants~Location, scales = "free") + 
        theme_bw() + 
    plot_annotation(caption = 'Lead concentration above 300 leads to a decrease in diversity')
