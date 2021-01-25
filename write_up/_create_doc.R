
detach("package:sa4ss", unload = TRUE)
remotes::install_github("nwfsc-assess/sa4ss")
library(sa4ss)
library(r4ss)

devtools::load_all("C:/Users/Chantel.Wetzel/Documents/GitHub/r4ss")


# Specify the directory for the document
setwd("C:/Assessments/2021/squarespot_rockfish_2021/write_up")

# Create the needed items to generate the "right" template that would be based on the inputs here:
sa4ss::draft(authors = c("Chantel R. Wetzel", "Brian J. Langseth", "Jason M. Cope", "John E. Budrick"),
  			 species = "Squarespot Rockfish",
  			 latin = "Sebastes hopkinsi",
  			 coast = "US West - California",
  			 type = c("sa"),
  			 create_dir = FALSE,
  			 edit = FALSE)

# Create a model Rdata object
# Create a model Rdata object
sa4ss::read_model(
          		  mod_loc = "C:/Assessments/2021/squarespot_rockfish_2021/models/2.3_hessian",
				  create_plots = TRUE, 
				  save_loc = file.path(getwd(), "tex_tables"),
				  verbose = TRUE)
#source("C:/Users/Chantel.Wetzel/Documents/GitHub/sa4ss/R/es_table_tex.R")

sa4ss::es_table_tex(dir = mod_loc, 
            save_loc = file.path(getwd(), "tex_tables"), 
            csv_name = "table_labels.csv")

# Read and create tex files for tables listed in "table" folder in the doc
sa4ss::es_table_tex(dir = file.path(getwd(), 'tables'), 
            save_loc = file.path(getwd(), "tex_tables"), 
            csv_name = "all_tables.csv")

if(file.exists("_main.Rmd")){
	file.remove("_main.Rmd")
}
# Render the pdf
bookdown::render_book("00a.Rmd", clean=FALSE, output_dir = getwd())




# Use to only render a specific section which can be quicker
bookdown::preview_chapter("01executive.Rmd", preview = TRUE, clean = FALSE)