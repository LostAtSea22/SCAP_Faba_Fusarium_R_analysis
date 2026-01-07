# configuring git in RStudio ----------------------------------------------
library(usethis)
use_git_config(user.name = "Laura Los",
               user.email = "losl@myumanitoba.ca")

# GitHub Set Up -----------------------------------------------------------
library(gitcreds)
library(usethis)

# Connect RStudio to your GitHub
create_github_token()
gitcreds::gitcreds_set()

# Put a project onto GitHub
use_github()

# Create a readme
use_readme_md()
