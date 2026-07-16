# R options for docker building
options(Ncpus = 8L)
options(timeout = 99999)  ## download time.out

# Configure BioCManager to use Posit Package Manager:
options(BioC_mirror = "https://packagemanager.posit.co/bioconductor")
options(BIOCONDUCTOR_CONFIG_FILE = "https://packagemanager.posit.co/bioconductor/config.yaml")

# Configure a CRAN snapshot compatible with Bioconductor 3.22.
# NOTE: the __linux__/<distro> segment must match the base image, otherwise
# PPM silently serves source packages instead of binaries.
options(repos = c(CRAN = "https://packagemanager.posit.co/cran/__linux__/trixie/latest"))

# PPM picks the binary build from the User-Agent. R does not send its own
# identity unless HTTPUserAgent is set explicitly -- it sends "libcurl/x.y.z",
# which PPM answers with source packages. Without this line the __linux__/
# trixie URL above has no effect.
options(HTTPUserAgent = sprintf(
  "R/%s R (%s)", getRversion(),
  paste(getRversion(), R.version["platform"], R.version["arch"], R.version["os"])
))
