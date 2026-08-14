# R options for docker building
options(Ncpus = 8L)
options(timeout = 99999)  ## download time.out

# Configure BioCManager to use Posit Package Manager:
options(BioC_mirror = "https://packagemanager.posit.co/bioconductor")
options(BIOCONDUCTOR_CONFIG_FILE = "https://packagemanager.posit.co/bioconductor/config.yaml")

# Configure a CRAN snapshot compatible with Bioconductor 3.22.
# NOTE: the __linux__/<distro> segment must match the base image, otherwise
# PPM silently serves source packages instead of binaries.
#
# Pinned to a date, NOT `latest`. On `latest` the binary channel served lattice
# 0.23-1, which dropped the `parallel` export that methylumi 2.56.0 still
# imports; that cascades methylumi -> lumi -> wateRmelon, and wateRmelon is a
# playbase Import, so playbase itself then fails to install. This snapshot
# still serves binaries, so the PPM speedup is retained.
# Keep this date in step with omicsplayground's dev/Rprofile -- if the two
# disagree, the opg build stages upgrade packages this base was compiled against.
options(repos = c(CRAN = "https://packagemanager.posit.co/cran/__linux__/trixie/2026-07-15"))

# PPM picks the binary build from the User-Agent. R does not send its own
# identity unless HTTPUserAgent is set explicitly -- it sends "libcurl/x.y.z",
# which PPM answers with source packages. Without this line the __linux__/
# trixie URL above has no effect.
options(HTTPUserAgent = sprintf(
  "R/%s R (%s)", getRversion(),
  paste(getRversion(), R.version["platform"], R.version["arch"], R.version["os"])
))
