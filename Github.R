# https://happygitwithr.com/hello-git


#### 1.  Introduce yourself to Git #### 

## in shell 
> git config --global user.name "So-Yeon Kim"
> git config --global user.name "kim.soyeon.bio@gmail.com"
  #  email associated with your GitHub account.


## in R 
# install.packages("usethis")
library(usethis)
use_git_config(user.name = "So-Yeon Kim", user.email = "kim.soyeon.bio@gmail.com")

## in shell 
> git config --global core.editor "emacs"
> git config --global init.defaultBranch main

## in R 
usethis::git_default_branch_configure()


#### Personal access token for HTTPS
usethis::create_github_token()
usethis::github_token()
gh::gh_token()
