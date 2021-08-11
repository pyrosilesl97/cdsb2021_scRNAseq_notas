library("gitcreds", "gert")
library('gh')

#Revisar la configuración de Git
gh::gh_whoami() # para checar cómo quedó la configuración

gert::git_commit_all("Nuevos")
gert::git_push()

usethis::use_r('05-control-de-calidad.R')
