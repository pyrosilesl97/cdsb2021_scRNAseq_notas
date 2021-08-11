library("gitcreds", "gert")
library('gh')

#Revisar la configuración de Git
gh::gh_whoami() # para checar cómo quedó la configuración

gert::git_commit_all("Nuevos")
gert::git_push()

usethis::use_r('07-seleccion-de-genes.R')

##Change the error
#options(error= function() starwarssay::say("Do or do not... There is no try", by = "baby_yoda"))

options(prompt = "🦖 ~ 🧬 >")

starwarssay::say("Join the dark side.... #CODE \n Welcome back Pablo", by = "darth_vader")


if (interactive()) {
  .First <- function() {
    # aquí van los mensajes o las funciones para la bienvenida
  }

  .Last <- function() {
    # aquí van los mensajes o las funciones para la despedida
  }
}
