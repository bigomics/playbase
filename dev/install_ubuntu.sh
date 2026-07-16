export DEBIAN_FRONTEND=noninteractive
export LANG=en_US.UTF-8
export LANGUAGE=en_US:en
export LC_ALL=en_US.UTF-8 

apt update && apt install -y \
    locales apt-utils \
    libcurl4-gnutls-dev libnode-dev libv8-dev \
    libssl-dev libxml2-dev libjpeg-dev \
    libgl-dev libglu-dev tk-dev libhdf5-dev \
    libgit2-dev libssh2-1-dev libnetcdf-dev \
    libudunits2-dev libgdal-dev libproj-dev libbz2-dev \
    jags cmake git procps htop \
    python3 python3-pip python-is-python3 \
    pdftk vim-tiny less wget gdebi-core \
    pandoc imagemagick libmagick++-dev libfftw3-dev libglpk-dev \
    libgsl-dev librsvg2-dev libgsl-dev curl \
    libsodium-dev libnlopt-dev \
    libharfbuzz-dev libfribidi-dev \
    fonts-lato

# Install Chrome (apt-key no longer exists on trixie, so use a keyring)
apt-get update && apt-get install -y wget gnupg2
install -d -m 0755 /etc/apt/keyrings
wget -q -O - https://dl-ssl.google.com/linux/linux_signing_key.pub \
    | gpg --dearmor -o /etc/apt/keyrings/google-chrome.gpg
echo "deb [arch=amd64 signed-by=/etc/apt/keyrings/google-chrome.gpg] http://dl.google.com/linux/chrome/deb/ stable main" > /etc/apt/sources.list.d/google-chrome.list
apt-get update && apt-get install -y google-chrome-stable

# disable c++ version warning (msa package). Glob the gcc version rather
# than hardcoding it: noble had gcc-13, trixie has gcc-14.
for h in /usr/include/c++/*/bits/c++0x_warning.h; do
    [ -f "$h" ] || continue
    mv "$h" "$h.DISABLED" && : > "$h"
done

## basic R
apt install -y r-base r-base-dev 
