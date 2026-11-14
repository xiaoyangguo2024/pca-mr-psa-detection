# Installs the MR packages into a project-local library, Rlib/.
lib <- file.path(getwd(), "Rlib")
.libPaths(c(lib, .libPaths()))
repos <- c(mrcieu = "https://mrcieu.r-universe.dev", CRAN = "https://cloud.r-project.org")
need <- c("TwoSampleMR", "MendelianRandomization", "MVMR", "mr.raps", "MRPRESSO", "coloc", "metafor")
miss <- need[!sapply(need, requireNamespace, quietly = TRUE)]
if (length(miss)) install.packages(miss, lib = lib, repos = repos, type = "binary", dependencies = c("Depends", "Imports", "LinkingTo"))
miss2 <- need[!sapply(need, requireNamespace, quietly = TRUE)]
if (length(miss2)) install.packages(miss2, lib = lib, repos = repos, dependencies = c("Depends", "Imports", "LinkingTo"))
print(sapply(need, requireNamespace, quietly = TRUE))
for (p in need) if (requireNamespace(p, quietly = TRUE)) cat(p, as.character(packageVersion(p)), "\n")
