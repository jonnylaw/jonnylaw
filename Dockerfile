FROM rocker/r-apt:disco

RUN apt-get update \
	&& apt-get install -y --no-install-recommends r-cran-rstan \
  libxml2-dev \
  libssl-dev \
  libcurl4-openssl-dev \
  && apt-get clean \
  && rm -rf /var/lib/apt/lists/

# Copy contents of repo to container
COPY . .

# Install Dependencies
RUN Rscript -e "install.packages('remotes', repos = 'https://demo.rstudiopm.com/all/__linux__/bionic/latest')"
RUN Rscript -e "remotes::install_deps(dependencies = TRUE, repos = 'https://demo.rstudiopm.com/all/__linux__/bionic/latest')"

# Compile Blog
RUN Rscript -e "blogdown::build_site()"
