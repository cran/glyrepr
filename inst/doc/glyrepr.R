## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(glyrepr)

## -----------------------------------------------------------------------------
# Just tell R what you have
glycan_composition(c(Hex = 5, HexNAc = 2), c(Gal = 1, GalNAc = 1))

## -----------------------------------------------------------------------------
# Perfect when you're processing data from files or databases
comp_list <- list(c(Hex = 5, HexNAc = 2), c(Gal = 1, GalNAc = 1))
as_glycan_composition(comp_list)

## -----------------------------------------------------------------------------
# Copy-paste from your mass spec software? No problem!
as_glycan_composition(c("Hex(5)HexNAc(2)", "H1N1"))

## -----------------------------------------------------------------------------
comp <- glycan_composition(
  c(Hex = 5, HexNAc = 2),          # generic sugars
  c(Gal = 1, Man = 1, GalNAc = 1)  # concrete sugars
)

# How many galactose residues?
count_mono(comp, "Gal")

# How many hexose residues? (This includes Gal and Man!)
count_mono(comp, "Hex")

## -----------------------------------------------------------------------------
iupacs <- c(
  "Man(a1-3)[Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-",  # The famous N-glycan core
  "Gal(b1-3)GalNAc(a1-",                                  # O-glycan core 1
  "Gal(b1-3)[GlcNAc(b1-6)]GalNAc(a1-",                    # O-glycan core 2
  "Man(a1-3)[Man(a1-6)]Man(a1-3)[Man(a1-6)]Man(a1-",      # A branched mannose tree
  "GlcNAc6Ac(b1-4)Glc3Me(a1-"                             # With some decorations
)

struc <- as_glycan_structure(iupacs)
struc

## -----------------------------------------------------------------------------
# Create a big dataset with lots of repetition
large_struc <- rep(struc, 1000)  # 5,000 structures total
large_struc

## -----------------------------------------------------------------------------
library(tictoc)

tic("Converting 5 structures")
result_small <- convert_to_generic(struc)
toc()

tic("Converting 5,000 structures")
result_large <- convert_to_generic(large_struc)
toc()

## -----------------------------------------------------------------------------
glycans <- as_glycan_structure(c(
  "Gal(b1-3)GalNAc(a1-",
  "Gal(b1-?)GalNAc(a1-",
  "Gal(??-?)GalNAc(??-",
  "Hex(??-?)HexNAc(??-",
  "Hex(b1-3)HexNAc(a1-"
))
get_structure_level(glycans)

## -----------------------------------------------------------------------------
remove_linkages(struc)

## -----------------------------------------------------------------------------
# Let's look at our decorated structure first
struc[5]

# Now remove the decorations (6Ac and 3Me)
remove_substituents(struc[5])

## -----------------------------------------------------------------------------
comp <- as_glycan_composition(struc)
comp

## -----------------------------------------------------------------------------
# Get the original string representations
as.character(struc)
as.character(comp)

## -----------------------------------------------------------------------------
suppressPackageStartupMessages(library(tibble))
suppressPackageStartupMessages(library(dplyr))

df <- tibble(
  id = seq_along(struc),
  structures = struc,
  names = c("N-glycan core", "Core 1", "Core 2", "Branched Man", "Decorated")
)

df %>% 
  mutate(n_man = count_mono(structures, "Man")) %>%
  filter(n_man > 1)

## -----------------------------------------------------------------------------
sessionInfo()

