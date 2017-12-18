pkgs <- c("psych","cvTools","robustbase","svd","tidyverse","rlang","svd","sigmoid","foreach","parallel","data.table","fanc","")

pkgs.required <- pkgs[!pkgs %in% utils::installed.packages()]
if (length(pkgs.required) != 0){
  utils::install.packages(pkgs.required)
}

source("simulation/functions.R")
source("simulation/simulate data LE.R")

library(fanc)
filter <- dplyr::filter