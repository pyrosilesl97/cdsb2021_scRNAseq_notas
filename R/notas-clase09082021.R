library("gitcreds", "gert")
library('gh')

#Revisar la configuración de Git
gh::gh_whoami() # para checar cómo quedó la configuración

usethis::use_r("02-Ejercicios_10082021.R")
