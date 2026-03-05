pkg_root <- normalizePath(
  file.path(testthat::test_path("..", "..")),
  winslash = "/",
  mustWork = TRUE
)

sys.source(file.path(pkg_root, "R", "XgeneR_class.R"), envir = globalenv())
sys.source(file.path(pkg_root, "R", "XgeneR_plotting.R"), envir = globalenv())
