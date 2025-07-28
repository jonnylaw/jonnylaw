FROM rocker/r-apt:disco

# 
RUN apt-get update \
	&& apt-get install -y --no-install-recommends r-cran-rstan \
	gdebi-core \
  libxml2-dev \
  libssl-dev \
  libcurl4-openssl-dev \
  && apt-get clean \
  && rm -rf /var/lib/apt/lists/ \
  && wget https://github.com/jgm/pandoc/releases/download/2.7.3/pandoc-2.7.3-1-amd64.deb \
  && gdebi --non-interactive pandoc-2.7.3-1-amd64.deb

# Copy contents of repo to container
COPY . .

# Install Dependencies
RUN Rscript -e "install.packages('remotes')"
RUN Rscript -e "remotes::install_deps(dependencies = TRUE)"

# Test package
RUN Rscript  -e "testthat::test_dir(path = 'tests/testthat')"

# Install package
RUN Rscript -e "pkgbuild::build(path = '.', dest_path = '.'); install.packages(pkgs = 'jonnylaw_0.1.0.tar.gz', repos = NULL, type = 'source')"

# Build pkgdown website
RUN Rscript -e "pkgdown::build()"
