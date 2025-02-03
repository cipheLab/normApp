FROM rocker/rstudio:4.2.2



RUN R -e 'install.packages("BiocManager")'
RUN R -e 'BiocManager::install("flowCore")'
RUN R -e 'BiocManager::install("flowDensity")'

RUN R -e "install.packages(c('shiny','mclust','colorspace','DT','RColorBrewer','shinybusy','shinydashboard','ggplot2','car','parallel','future.apply'), repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('devtools')"


RUN mkdir /root/normApp
COPY normApp /root/normApp


RUN ls -la /root/normApp/


EXPOSE 3838


CMD ["R", "-e", "shiny::runApp('/root/normApp', host = '0.0.0.0', port = 3838)"]
