# One-time package restoration. Run deliberately, not before every lesson.
# Open the RStudio project first. Restart R after restoration completes.
if (!requireNamespace("renv", quietly = TRUE)) {
  install.packages("renv", repos = "https://cloud.r-project.org")
}
renv::load(project = ".")         # Restore this project's library, not a global one.
renv::restore(prompt = TRUE)
renv::status()
