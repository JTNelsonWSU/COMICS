\name{ui}
\alias{ui}
\docType{data}
\title{
%%   ~~ data name/kind ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of the dataset. ~~
}
\usage{data("ui")}
\format{
  The format is:
List of 3
 $ : NULL
 $ :List of 3
  ..$ name    : chr "head"
  ..$ attribs : Named list()
  ..$ children:List of 1
  .. ..$ :List of 3
  .. .. ..$ name    : chr "link"
  .. .. ..$ attribs :List of 3
  .. .. .. ..$ rel : chr "stylesheet"
  .. .. .. ..$ type: chr "text/css"
  .. .. .. ..$ href: chr "shinythemes/css/superhero.min.css"
  .. .. ..$ children: list()
  .. .. ..- attr(*, "class")= chr "shiny.tag"
  ..- attr(*, "class")= chr "shiny.tag"
 $ :List of 1
  ..$ :List of 3
  .. ..$ name    : chr "div"
  .. ..$ attribs :List of 1
  .. .. ..$ class: chr "container-fluid"
  .. ..$ children:List of 2
  .. .. ..$ :List of 2
  .. .. .. ..$ :List of 3
  .. .. .. .. ..$ name    : chr "head"
  .. .. .. .. ..$ attribs : Named list()
  .. .. .. .. ..$ children:List of 1
  .. .. .. .. .. ..$ :List of 3
  .. .. .. .. .. .. ..$ name    : chr "title"
  .. .. .. .. .. .. ..$ attribs : Named list()
  .. .. .. .. .. .. ..$ children:List of 1
  .. .. .. .. .. .. .. ..$ : chr "Thanks for using C.O.M.I.C.S"
  .. .. .. .. .. .. ..- attr(*, "class")= chr "shiny.tag"
  .. .. .. .. ..- attr(*, "class")= chr "shiny.tag"
  .. .. .. ..$ :List of 3
  .. .. .. .. ..$ name    : chr "h2"
  .. .. .. .. ..$ attribs : Named list()
  .. .. .. .. ..$ children:List of 1
  .. .. .. .. .. ..$ : chr "Thanks for using C.O.M.I.C.S"
  .. .. .. .. ..- attr(*, "class")= chr "shiny.tag"
  .. .. .. ..- attr(*, "class")= chr [1:2] "shiny.tag.list" "list"
  .. .. ..$ :List of 3
  .. .. .. ..$ name    : chr "div"
  .. .. .. ..$ attribs :List of 1
  .. .. .. .. ..$ class: chr "row"
  .. .. .. ..$ children:List of 2
  .. .. .. .. ..$ :List of 3
  .. .. .. .. .. ..$ name    : chr "div"
  .. .. .. .. .. ..$ attribs :List of 1
  .. .. .. .. .. .. ..$ class: chr "col-sm-4"
  .. .. .. .. .. ..$ children:List of 1
  .. .. .. .. .. .. ..$ :List of 3
  .. .. .. .. .. .. .. ..$ name    : chr "form"
  .. .. .. .. .. .. .. ..$ attribs :List of 1
  .. .. .. .. .. .. .. .. ..$ class: chr "well"
  .. .. .. .. .. .. .. ..$ children:List of 8
  .. .. .. .. .. .. .. .. ..$ :List of 3
  .. .. .. .. .. .. .. .. .. ..$ name    : chr "div"
  .. .. .. .. .. .. .. .. .. ..$ attribs :List of 1
  .. .. .. .. .. .. .. .. .. .. ..$ class: chr "form-group shiny-input-container"
  .. .. .. .. .. .. .. .. .. ..$ children:List of 3
  .. .. .. .. .. .. .. .. .. .. ..$ :List of 3
  .. .. .. .. .. .. .. .. .. .. .. ..$ name    : chr "label"
  .. .. .. .. .. .. .. .. .. .. .. ..$ attribs :List of 2
  .. .. .. .. .. .. .. .. .. .. .. .. ..$ class: chr "control-label"
  .. .. .. .. .. .. .. .. .. .. .. .. ..$ for  : chr "file1"
  .. .. .. .. .. .. .. .. .. .. .. ..$ children:List of 1
  .. .. .. .. .. .. .. .. .. .. .. .. ..$ : chr "Data file"
  .. .. .. .. .. .. .. .. .. .. .. ..- attr(*, "class")= chr "shiny.tag"
  .. .. .. .. .. .. .. .. .. .. ..$ :List of 3
  .. .. .. .. .. .. .. .. .. .. .. ..$ name    : chr "div"
  .. .. .. .. .. .. .. .. .. .. .. ..$ attribs :List of 1
  .. .. .. .. .. .. .. .. .. .. .. .. ..$ class: chr "input-group"
  .. .. .. .. .. .. .. .. .. .. .. ..$ children:List of 2
  .. .. .. .. .. .. .. .. .. .. .. .. ..$ :List of 3
  .. .. .. .. .. .. .. .. .. .. .. .. .. ..$ name    : chr "label"
  .. .. .. .. .. .. .. .. .. .. .. .. .. ..$ attribs :List of 1
  .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..$ class: chr "input-group-btn"
  .. .. .. .. .. .. .. .. .. .. .. .. .. ..$ children:List of 1
  .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..$ :List of 3
  .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..$ name    : chr "span"
  .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..$ attribs :List of 1
  .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..$ class: chr "btn btn-default btn-file"
  .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..$ children:List of 2
  .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..$ : chr "Browse"
  .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..$ :List of 3
  .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..$ name    : chr "input"
  .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..$ attribs :List of 5
  .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..$ id    : chr "file1"
  .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..$ name  : chr "file1"
  .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..$ type  : chr "file"
  .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..$ style : chr "display: none;"
  .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..$ accept: chr "text/csv,text/tab-seperated-vlaues,text/plain,.csv"
  .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..$ children: list()
  .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..- attr(*, "class")= chr "shiny.tag"
  .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..- attr(*, "class")= chr "shiny.tag"
  .. .. .. .. .. .. .. .. .. .. .. .. .. ..- attr(*, "class")= chr "shiny.tag"
  .. .. .. .. .. .. .. .. .. .. .. .. ..$ :List of 3
  .. .. .. .. .. .. .. .. .. .. .. .. .. ..$ name    : chr "input"
  .. .. .. .. .. .. .. .. .. .. .. .. .. ..$ attribs :List of 4
  .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..$ type       : chr "text"
  .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..$ class      : chr "form-control"
  .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..$ placeholder: chr "No file selected"
  .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..$ readonly   : chr "readonly"
  .. .. .. .. .. .. .. .. .. .. .. .. .. ..$ children: list()
  .. .. .. .. .. .. .. .. .. .. .. .. .. ..- attr(*, "class")= chr "shiny.tag"
  .. .. .. .. .. .. .. .. .. .. .. ..- attr(*, "class")= chr "shiny.tag"
  .. .. .. .. .. .. .. .. .. .. ..$ :List of 3
  .. .. .. .. .. .. .. .. .. .. .. ..$ name    : chr "div"
  .. .. .. .. .. .. .. .. .. .. .. ..$ attribs :List of 2
  .. .. .. .. .. .. .. .. .. .. .. .. ..$ id   : chr "file1_progress"
  .. .. .. .. .. .. .. .. .. .. .. .. ..$ class: chr "progress progress-striped active shiny-file-input-progress"
  .. .. .. .. .. .. .. .. .. .. .. ..$ children:List of 1
  .. .. .. .. .. .. .. .. .. .. .. .. ..$ :List of 3
  .. .. .. .. .. .. .. .. .. .. .. .. .. ..$ name    : chr "div"
  .. .. .. .. .. .. .. .. .. .. .. .. .. ..$ attribs :List of 1
  .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..$ class: chr "progress-bar"
  .. .. .. .. .. .. .. .. .. .. .. .. .. ..$ children: list()
  .. .. .. .. .. .. .. .. .. .. .. .. .. ..- attr(*, "class")= chr "shiny.tag"
  .. .. .. .. .. .. .. .. .. .. .. ..- attr(*, "class")= chr "shiny.tag"
  .. .. .. .. .. .. .. .. .. ..- attr(*, "class")= chr "shiny.tag"
  .. .. .. .. .. .. .. .. ..$ :List of 3
  .. .. .. .. .. .. .. .. .. ..$ name    : chr "div"
  .. .. .. .. .. .. .. .. .. ..$ attribs :List of 1
  .. .. .. .. .. .. .. .. .. .. ..$ class: chr "form-group shiny-input-container"
  .. .. .. .. .. .. .. .. .. ..$ children:List of 3
  .. .. .. .. .. .. .. .. .. .. ..$ :List of 3
  .. .. .. .. .. .. .. .. .. .. .. ..$ name    : chr "label"
  .. .. .. .. .. .. .. .. .. .. .. ..$ attribs :List of 2
  .. .. .. .. .. .. .. .. .. .. .. .. ..$ class: chr "control-label"
  .. .. .. .. .. .. .. .. .. .. .. .. ..$ for  : chr "configuration"
  .. .. .. .. .. .. .. .. .. .. .. ..$ children:List of 1
  .. .. .. .. .. .. .. .. .. .. .. .. ..$ : chr "Genome Configuration File"
  .. .. .. .. .. .. .. .. .. .. .. ..- attr(*, "class")= chr "shiny.tag"
  .. .. .. .. .. .. .. .. .. .. ..$ :List of 3
  .. .. .. .. .. .. .. .. .. .. .. ..$ name    : chr "div"
  .. .. .. .. .. .. .. .. .. .. .. ..$ attribs :List of 1
  .. .. .. .. .. .. .. .. .. .. .. .. ..$ class: chr "input-group"
  .. .. .. .. .. .. .. .. .. .. .. ..$ children:List of 2
  .. .. .. .. .. .. .. .. .. .. .. .. ..$ :List of 3
  .. .. .. .. .. .. .. .. .. .. .. .. .. ..$ name    : chr "label"
  .. .. .. .. .. .. .. .. .. .. .. .. .. ..$ attribs :List of 1
  .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..$ class: chr "input-group-btn"
  .. .. .. .. .. .. .. .. .. .. .. .. .. ..$ children:List of 1
  .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..$ :List of 3
  .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..$ name    : chr "span"
  .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..$ attribs :List of 1
  .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..$ class: chr "btn btn-default btn-file"
  .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..$ children:List of 2
  .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..$ : chr "Browse"
  .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..$ :List of 3
  .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..$ name    : chr "input"
  .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..$ attribs :List of 5
  .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..$ id    : chr "configuration"
  .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..$ name  : chr "configuration"
  .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..$ type  : chr "file"
  .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..$ style : chr "display: none;"
  .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..$ accept: chr "text/csv,text/tab-seperated-values,text/plain,.csv"
  .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..$ children: list()
  .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..- attr(*, "class")= chr "shiny.tag"
  .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..- attr(*, "class")= chr "shiny.tag"
  .. .. .. .. .. .. .. .. .. .. .. .. .. ..- attr(*, "class")= chr "shiny.tag"
  .. .. .. .. .. .. .. .. .. .. .. .. ..$ :List of 3
  .. .. .. .. .. .. .. .. .. .. .. .. .. ..$ name    : chr "input"
  .. .. .. .. .. .. .. .. .. .. .. .. .. ..$ attribs :List of 4
  .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..$ type       : chr "text"
  .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..$ class      : chr "form-control"
  .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..$ placeholder: chr "No file selected"
  .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..$ readonly   : chr "readonly"
  .. .. .. .. .. .. .. .. .. .. .. .. .. ..$ children: list()
  .. .. .. .. .. .. .. .. .. .. .. .. .. ..- attr(*, "class")= chr "shiny.tag"
  .. .. .. .. .. .. .. .. .. .. .. ..- attr(*, "class")= chr "shiny.tag"
  .. .. .. .. .. .. .. .. .. .. ..$ :List of 3
  .. .. .. .. .. .. .. .. .. .. .. ..$ name    : chr "div"
  .. .. .. .. .. .. .. .. .. .. .. ..$ attribs :List of 2
  .. .. .. .. .. .. .. .. .. .. .. .. ..$ id   : chr "configuration_progress"
  .. .. .. .. .. .. .. .. .. .. .. .. ..$ class: chr "progress progress-striped active shiny-file-input-progress"
  .. .. .. .. .. .. .. .. .. .. .. ..$ children:List of 1
  .. .. .. .. .. .. .. .. .. .. .. .. ..$ :List of 3
  .. .. .. .. .. .. .. .. .. .. .. .. .. ..$ name    : chr "div"
  .. .. .. .. .. .. .. .. .. .. .. .. .. ..$ attribs :List of 1
  .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..$ class: chr "progress-bar"
  .. .. .. .. .. .. .. .. .. .. .. .. .. ..$ children: list()
  .. .. .. .. .. .. .. .. .. .. .. .. .. ..- attr(*, "class")= chr "shiny.tag"
  .. .. .. .. .. .. .. .. .. .. .. ..- attr(*, "class")= chr "shiny.tag"
  .. .. .. .. .. .. .. .. .. ..- attr(*, "class")= chr "shiny.tag"
  .. .. .. .. .. .. .. .. ..$ :List of 3
  .. .. .. .. .. .. .. .. .. ..$ name    : chr "div"
  .. .. .. .. .. .. .. .. .. ..$ attribs :List of 1
  .. .. .. .. .. .. .. .. .. .. ..$ class: chr "form-group shiny-input-container"
  .. .. .. .. .. .. .. .. .. ..$ children:List of 2
  .. .. .. .. .. .. .. .. .. .. ..$ :List of 3
  .. .. .. .. .. .. .. .. .. .. .. ..$ name    : chr "label"
  .. .. .. .. .. .. .. .. .. .. .. ..$ attribs :List of 2
  .. .. .. .. .. .. .. .. .. .. .. .. ..$ class: chr "control-label"
  .. .. .. .. .. .. .. .. .. .. .. .. ..$ for  : chr "Cutoff"
  .. .. .. .. .. .. .. .. .. .. .. ..$ children:List of 1
  .. .. .. .. .. .. .. .. .. .. .. .. ..$ : chr "Statistical Cutoff"
  .. .. .. .. .. .. .. .. .. .. .. ..- attr(*, "class")= chr "shiny.tag"
  .. .. .. .. .. .. .. .. .. .. ..$ :List of 3
  .. .. .. .. .. .. .. .. .. .. .. ..$ name    : chr "input"
  .. .. .. .. .. .. .. .. .. .. .. ..$ attribs :List of 6
  .. .. .. .. .. .. .. .. .. .. .. .. ..$ id   : chr "Cutoff"
  .. .. .. .. .. .. .. .. .. .. .. .. ..$ type : chr "number"
  .. .. .. .. .. .. .. .. .. .. .. .. ..$ class: chr "form-control"
  .. .. .. .. .. .. .. .. .. .. .. .. ..$ value: chr "5"
  .. .. .. .. .. .. .. .. .. .. .. .. ..$ min  : num 0
  .. .. .. .. .. .. .. .. .. .. .. .. ..$ max  : num 10
  .. .. .. .. .. .. .. .. .. .. .. ..$ children: list()
  .. .. .. .. .. .. .. .. .. .. .. ..- attr(*, "class")= chr "shiny.tag"
  .. .. .. .. .. .. .. .. .. ..- attr(*, "class")= chr "shiny.tag"
  .. .. .. .. .. .. .. .. ..$ :List of 3
  .. .. .. .. .. .. .. .. .. ..$ name    : chr "div"
  .. .. .. .. .. .. .. .. .. ..$ attribs :List of 1
  .. .. .. .. .. .. .. .. .. .. ..$ class: chr "form-group shiny-input-container"
  .. .. .. .. .. .. .. .. .. ..$ children:List of 2
  .. .. .. .. .. .. .. .. .. .. ..$ :List of 3
  .. .. .. .. .. .. .. .. .. .. .. ..$ name    : chr "label"
  .. .. .. .. .. .. .. .. .. .. .. ..$ attribs :List of 2
  .. .. .. .. .. .. .. .. .. .. .. .. ..$ class: chr "control-label"
  .. .. .. .. .. .. .. .. .. .. .. .. ..$ for  : chr "dataset"
  .. .. .. .. .. .. .. .. .. .. .. ..$ children:List of 1
  .. .. .. .. .. .. .. .. .. .. .. .. ..$ : chr "Choose a dataset"
  .. .. .. .. .. .. .. .. .. .. .. ..- attr(*, "class")= chr "shiny.tag"
  .. .. .. .. .. .. .. .. .. .. ..$ :List of 3
  .. .. .. .. .. .. .. .. .. .. .. ..$ name    : chr "div"
  .. .. .. .. .. .. .. .. .. .. .. ..$ attribs : Named list()
  .. .. .. .. .. .. .. .. .. .. .. ..$ children:List of 2
  .. .. .. .. .. .. .. .. .. .. .. .. ..$ :List of 3
  .. .. .. .. .. .. .. .. .. .. .. .. .. ..$ name    : chr "select"
  .. .. .. .. .. .. .. .. .. .. .. .. .. ..$ attribs :List of 1
  .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..$ id: chr "dataset"
  .. .. .. .. .. .. .. .. .. .. .. .. .. ..$ children:List of 1
  .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..$ : 'html' chr "<option value=\"ICS Distance\" selected>ICS Distance</option>"
  .. .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..- attr(*, "html")= logi TRUE
  .. .. .. .. .. .. .. .. .. .. .. .. .. ..- attr(*, "class")= chr "shiny.tag"
  .. .. .. .. .. .. .. .. .. .. .. .. ..$ :List of 3
  .. .. .. .. .. .. .. .. .. .. .. .. .. ..$ name    : chr "script"
  .. .. .. .. .. .. .. .. .. .. .. .. .. ..$ attribs :List of 3
  .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..$ type         : chr "application/json"
  .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..$ data-for     : chr "dataset"
  .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..$ data-nonempty: chr ""
  .. .. .. .. .. .. .. .. .. .. .. .. .. ..$ children:List of 1
  .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..$ : chr "{}"
  .. .. .. .. .. .. .. .. .. .. .. .. .. ..- attr(*, "class")= chr "shiny.tag"
  .. .. .. .. .. .. .. .. .. .. .. ..- attr(*, "class")= chr "shiny.tag"
  .. .. .. .. .. .. .. .. .. ..- attr(*, "class")= chr "shiny.tag"
  .. .. .. .. .. .. .. .. .. ..- attr(*, "html_dependencies")=List of 1
  .. .. .. .. .. .. .. .. .. .. ..$ :List of 10
  .. .. .. .. .. .. .. .. .. .. .. ..$ name      : chr "selectize"
  .. .. .. .. .. .. .. .. .. .. .. ..$ version   : chr "0.11.2"
  .. .. .. .. .. .. .. .. .. .. .. ..$ src       :List of 1
  .. .. .. .. .. .. .. .. .. .. .. .. ..$ href: chr "shared/selectize"
  .. .. .. .. .. .. .. .. .. .. .. ..$ meta      : NULL
  .. .. .. .. .. .. .. .. .. .. .. ..$ script    : NULL
  .. .. .. .. .. .. .. .. .. .. .. ..$ stylesheet: chr "css/selectize.bootstrap3.css"
  .. .. .. .. .. .. .. .. .. .. .. ..$ head      : chr "<!--[if lt IE 9]>\n<script src=\"shared/selectize/js/es5-shim.min.js\"></script>\n<![endif]-->\n<script src=\"s"| __truncated__
  .. .. .. .. .. .. .. .. .. .. .. ..$ attachment: NULL
  .. .. .. .. .. .. .. .. .. .. .. ..$ package   : NULL
  .. .. .. .. .. .. .. .. .. .. .. ..$ all_files : logi TRUE
  .. .. .. .. .. .. .. .. .. .. .. ..- attr(*, "class")= chr "html_dependency"
  .. .. .. .. .. .. .. .. ..$ :List of 3
  .. .. .. .. .. .. .. .. .. ..$ name    : chr "div"
  .. .. .. .. .. .. .. .. .. ..$ attribs :List of 1
  .. .. .. .. .. .. .. .. .. .. ..$ class: chr "form-group shiny-input-container"
  .. .. .. .. .. .. .. .. .. ..$ children:List of 2
  .. .. .. .. .. .. .. .. .. .. ..$ :List of 3
  .. .. .. .. .. .. .. .. .. .. .. ..$ name    : chr "label"
  .. .. .. .. .. .. .. .. .. .. .. ..$ attribs :List of 2
  .. .. .. .. .. .. .. .. .. .. .. .. ..$ class: chr "control-label"
  .. .. .. .. .. .. .. .. .. .. .. .. ..$ for  : chr "TestOfInterest"
  .. .. .. .. .. .. .. .. .. .. .. ..$ children:List of 1
  .. .. .. .. .. .. .. .. .. .. .. .. ..$ : chr "Test of interest"
  .. .. .. .. .. .. .. .. .. .. .. ..- attr(*, "class")= chr "shiny.tag"
  .. .. .. .. .. .. .. .. .. .. ..$ :List of 3
  .. .. .. .. .. .. .. .. .. .. .. ..$ name    : chr "input"
  .. .. .. .. .. .. .. .. .. .. .. ..$ attribs :List of 6
  .. .. .. .. .. .. .. .. .. .. .. .. ..$ id   : chr "TestOfInterest"
  .. .. .. .. .. .. .. .. .. .. .. .. ..$ type : chr "number"
  .. .. .. .. .. .. .. .. .. .. .. .. ..$ class: chr "form-control"
  .. .. .. .. .. .. .. .. .. .. .. .. ..$ value: chr "1"
  .. .. .. .. .. .. .. .. .. .. .. .. ..$ min  : num 0
  .. .. .. .. .. .. .. .. .. .. .. .. ..$ max  : num 100
  .. .. .. .. .. .. .. .. .. .. .. ..$ children: list()
  .. .. .. .. .. .. .. .. .. .. .. ..- attr(*, "class")= chr "shiny.tag"
  .. .. .. .. .. .. .. .. .. ..- attr(*, "class")= chr "shiny.tag"
  .. .. .. .. .. .. .. .. ..$ :List of 3
  .. .. .. .. .. .. .. .. .. ..$ name    : chr "div"
  .. .. .. .. .. .. .. .. .. ..$ attribs :List of 1
  .. .. .. .. .. .. .. .. .. .. ..$ class: chr "form-group shiny-input-container"
  .. .. .. .. .. .. .. .. .. ..$ children:List of 2
  .. .. .. .. .. .. .. .. .. .. ..$ :List of 3
  .. .. .. .. .. .. .. .. .. .. .. ..$ name    : chr "label"
  .. .. .. .. .. .. .. .. .. .. .. ..$ attribs :List of 2
  .. .. .. .. .. .. .. .. .. .. .. .. ..$ class: chr "control-label"
  .. .. .. .. .. .. .. .. .. .. .. .. ..$ for  : chr "Chromosomes"
  .. .. .. .. .. .. .. .. .. .. .. ..$ children:List of 1
  .. .. .. .. .. .. .. .. .. .. .. .. ..$ : chr "Chromosome of interest"
  .. .. .. .. .. .. .. .. .. .. .. ..- attr(*, "class")= chr "shiny.tag"
  .. .. .. .. .. .. .. .. .. .. ..$ :List of 3
  .. .. .. .. .. .. .. .. .. .. .. ..$ name    : chr "input"
  .. .. .. .. .. .. .. .. .. .. .. ..$ attribs :List of 6
  .. .. .. .. .. .. .. .. .. .. .. .. ..$ id   : chr "Chromosomes"
  .. .. .. .. .. .. .. .. .. .. .. .. ..$ type : chr "number"
  .. .. .. .. .. .. .. .. .. .. .. .. ..$ class: chr "form-control"
  .. .. .. .. .. .. .. .. .. .. .. .. ..$ value: chr "1"
  .. .. .. .. .. .. .. .. .. .. .. .. ..$ min  : num 0
  .. .. .. .. .. .. .. .. .. .. .. .. ..$ max  : chr "configuration"
  .. .. .. .. .. .. .. .. .. .. .. ..$ children: list()
  .. .. .. .. .. .. .. .. .. .. .. ..- attr(*, "class")= chr "shiny.tag"
  .. .. .. .. .. .. .. .. .. ..- attr(*, "class")= chr "shiny.tag"
  .. .. .. .. .. .. .. .. ..$ :List of 3
  .. .. .. .. .. .. .. .. .. ..$ name    : chr "a"
  .. .. .. .. .. .. .. .. .. ..$ attribs :List of 5
  .. .. .. .. .. .. .. .. .. .. ..$ id      : chr "downloadData"
  .. .. .. .. .. .. .. .. .. .. ..$ class   : chr "btn btn-default shiny-download-link "
  .. .. .. .. .. .. .. .. .. .. ..$ href    : chr ""
  .. .. .. .. .. .. .. .. .. .. ..$ target  : chr "_blank"
  .. .. .. .. .. .. .. .. .. .. ..$ download: logi NA
  .. .. .. .. .. .. .. .. .. ..$ children:List of 2
  .. .. .. .. .. .. .. .. .. .. ..$ :List of 3
  .. .. .. .. .. .. .. .. .. .. .. ..$ name    : chr "i"
  .. .. .. .. .. .. .. .. .. .. .. ..$ attribs :List of 1
  .. .. .. .. .. .. .. .. .. .. .. .. ..$ class: chr "fa fa-download"
  .. .. .. .. .. .. .. .. .. .. .. ..$ children: list()
  .. .. .. .. .. .. .. .. .. .. .. ..- attr(*, "class")= chr "shiny.tag"
  .. .. .. .. .. .. .. .. .. .. .. ..- attr(*, "html_dependencies")=List of 1
  .. .. .. .. .. .. .. .. .. .. .. .. ..$ :List of 10
  .. .. .. .. .. .. .. .. .. .. .. .. .. ..$ name      : chr "font-awesome"
  .. .. .. .. .. .. .. .. .. .. .. .. .. ..$ version   : chr "5.3.1"
  .. .. .. .. .. .. .. .. .. .. .. .. .. ..$ src       :List of 1
  .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..$ file: chr "www/shared/fontawesome"
  .. .. .. .. .. .. .. .. .. .. .. .. .. ..$ meta      : NULL
  .. .. .. .. .. .. .. .. .. .. .. .. .. ..$ script    : NULL
  .. .. .. .. .. .. .. .. .. .. .. .. .. ..$ stylesheet: chr [1:2] "css/all.min.css" "css/v4-shims.min.css"
  .. .. .. .. .. .. .. .. .. .. .. .. .. ..$ head      : NULL
  .. .. .. .. .. .. .. .. .. .. .. .. .. ..$ attachment: NULL
  .. .. .. .. .. .. .. .. .. .. .. .. .. ..$ package   : chr "shiny"
  .. .. .. .. .. .. .. .. .. .. .. .. .. ..$ all_files : logi TRUE
  .. .. .. .. .. .. .. .. .. .. .. .. .. ..- attr(*, "class")= chr "html_dependency"
  .. .. .. .. .. .. .. .. .. .. .. ..- attr(*, "browsable_html")= logi TRUE
  .. .. .. .. .. .. .. .. .. .. ..$ : chr "Download ICS output"
  .. .. .. .. .. .. .. .. .. ..- attr(*, "class")= chr "shiny.tag"
  .. .. .. .. .. .. .. .. ..$ :List of 3
  .. .. .. .. .. .. .. .. .. ..$ name    : chr "a"
  .. .. .. .. .. .. .. .. .. ..$ attribs :List of 5
  .. .. .. .. .. .. .. .. .. .. ..$ id      : chr "downloadPlot"
  .. .. .. .. .. .. .. .. .. .. ..$ class   : chr "btn btn-default shiny-download-link "
  .. .. .. .. .. .. .. .. .. .. ..$ href    : chr ""
  .. .. .. .. .. .. .. .. .. .. ..$ target  : chr "_blank"
  .. .. .. .. .. .. .. .. .. .. ..$ download: logi NA
  .. .. .. .. .. .. .. .. .. ..$ children:List of 2
  .. .. .. .. .. .. .. .. .. .. ..$ :List of 3
  .. .. .. .. .. .. .. .. .. .. .. ..$ name    : chr "i"
  .. .. .. .. .. .. .. .. .. .. .. ..$ attribs :List of 1
  .. .. .. .. .. .. .. .. .. .. .. .. ..$ class: chr "fa fa-download"
  .. .. .. .. .. .. .. .. .. .. .. ..$ children: list()
  .. .. .. .. .. .. .. .. .. .. .. ..- attr(*, "class")= chr "shiny.tag"
  .. .. .. .. .. .. .. .. .. .. .. ..- attr(*, "html_dependencies")=List of 1
  .. .. .. .. .. .. .. .. .. .. .. .. ..$ :List of 10
  .. .. .. .. .. .. .. .. .. .. .. .. .. ..$ name      : chr "font-awesome"
  .. .. .. .. .. .. .. .. .. .. .. .. .. ..$ version   : chr "5.3.1"
  .. .. .. .. .. .. .. .. .. .. .. .. .. ..$ src       :List of 1
  .. .. .. .. .. .. .. .. .. .. .. .. .. .. ..$ file: chr "www/shared/fontawesome"
  .. .. .. .. .. .. .. .. .. .. .. .. .. ..$ meta      : NULL
  .. .. .. .. .. .. .. .. .. .. .. .. .. ..$ script    : NULL
  .. .. .. .. .. .. .. .. .. .. .. .. .. ..$ stylesheet: chr [1:2] "css/all.min.css" "css/v4-shims.min.css"
  .. .. .. .. .. .. .. .. .. .. .. .. .. ..$ head      : NULL
  .. .. .. .. .. .. .. .. .. .. .. .. .. ..$ attachment: NULL
  .. .. .. .. .. .. .. .. .. .. .. .. .. ..$ package   : chr "shiny"
  .. .. .. .. .. .. .. .. .. .. .. .. .. ..$ all_files : logi TRUE
  .. .. .. .. .. .. .. .. .. .. .. .. .. ..- attr(*, "class")= chr "html_dependency"
  .. .. .. .. .. .. .. .. .. .. .. ..- attr(*, "browsable_html")= logi TRUE
  .. .. .. .. .. .. .. .. .. .. ..$ : chr "Download ICS figure"
  .. .. .. .. .. .. .. .. .. ..- attr(*, "class")= chr "shiny.tag"
  .. .. .. .. .. .. .. ..- attr(*, "class")= chr "shiny.tag"
  .. .. .. .. .. ..- attr(*, "class")= chr "shiny.tag"
  .. .. .. .. ..$ :List of 3
  .. .. .. .. .. ..$ name    : chr "div"
  .. .. .. .. .. ..$ attribs :List of 1
  .. .. .. .. .. .. ..$ class: chr "col-sm-8"
  .. .. .. .. .. ..$ children:List of 10
  .. .. .. .. .. .. ..$ :List of 3
  .. .. .. .. .. .. .. ..$ name    : chr "h4"
  .. .. .. .. .. .. .. ..$ attribs : Named list()
  .. .. .. .. .. .. .. ..$ children:List of 1
  .. .. .. .. .. .. .. .. ..$ : chr "Summary.ICS"
  .. .. .. .. .. .. .. ..- attr(*, "class")= chr "shiny.tag"
  .. .. .. .. .. .. ..$ :List of 3
  .. .. .. .. .. .. .. ..$ name    : chr "pre"
  .. .. .. .. .. .. .. ..$ attribs :List of 2
  .. .. .. .. .. .. .. .. ..$ id   : chr "summary.ICS"
  .. .. .. .. .. .. .. .. ..$ class: chr "shiny-text-output noplaceholder"
  .. .. .. .. .. .. .. ..$ children: list()
  .. .. .. .. .. .. .. ..- attr(*, "class")= chr "shiny.tag"
  .. .. .. .. .. .. ..$ :List of 3
  .. .. .. .. .. .. .. ..$ name    : chr "div"
  .. .. .. .. .. .. .. ..$ attribs :List of 3
  .. .. .. .. .. .. .. .. ..$ id   : chr "ICS.hist"
  .. .. .. .. .. .. .. .. ..$ class: chr "shiny-plot-output"
  .. .. .. .. .. .. .. .. ..$ style: chr "width: 100% ; height: 400px"
  .. .. .. .. .. .. .. ..$ children: list()
  .. .. .. .. .. .. .. ..- attr(*, "class")= chr "shiny.tag"
  .. .. .. .. .. .. ..$ :List of 3
  .. .. .. .. .. .. .. ..$ name    : chr "div"
  .. .. .. .. .. .. .. ..$ attribs :List of 3
  .. .. .. .. .. .. .. .. ..$ id   : chr "ICS.chromosome.hist"
  .. .. .. .. .. .. .. .. ..$ class: chr "shiny-plot-output"
  .. .. .. .. .. .. .. .. ..$ style: chr "width: 100% ; height: 400px"
  .. .. .. .. .. .. .. ..$ children: list()
  .. .. .. .. .. .. .. ..- attr(*, "class")= chr "shiny.tag"
  .. .. .. .. .. .. ..$ :List of 3
  .. .. .. .. .. .. .. ..$ name    : chr "div"
  .. .. .. .. .. .. .. ..$ attribs :List of 3
  .. .. .. .. .. .. .. .. ..$ id   : chr "GenomeScan"
  .. .. .. .. .. .. .. .. ..$ class: chr "shiny-plot-output"
  .. .. .. .. .. .. .. .. ..$ style: chr "width: 100% ; height: 400px"
  .. .. .. .. .. .. .. ..$ children: list()
  .. .. .. .. .. .. .. ..- attr(*, "class")= chr "shiny.tag"
  .. .. .. .. .. .. ..$ :List of 3
  .. .. .. .. .. .. .. ..$ name    : chr "h3"
  .. .. .. .. .. .. .. ..$ attribs : Named list()
  .. .. .. .. .. .. .. ..$ children:List of 1
  .. .. .. .. .. .. .. .. ..$ : chr "Summary.SingleTest"
  .. .. .. .. .. .. .. ..- attr(*, "class")= chr "shiny.tag"
  .. .. .. .. .. .. ..$ :List of 3
  .. .. .. .. .. .. .. ..$ name    : chr "pre"
  .. .. .. .. .. .. .. ..$ attribs :List of 2
  .. .. .. .. .. .. .. .. ..$ id   : chr "summary.SingleTest"
  .. .. .. .. .. .. .. .. ..$ class: chr "shiny-text-output noplaceholder"
  .. .. .. .. .. .. .. ..$ children: list()
  .. .. .. .. .. .. .. ..- attr(*, "class")= chr "shiny.tag"
  .. .. .. .. .. .. ..$ :List of 3
  .. .. .. .. .. .. .. ..$ name    : chr "div"
  .. .. .. .. .. .. .. ..$ attribs :List of 3
  .. .. .. .. .. .. .. .. ..$ id   : chr "GenomeTest"
  .. .. .. .. .. .. .. .. ..$ class: chr "shiny-plot-output"
  .. .. .. .. .. .. .. .. ..$ style: chr "width: 100% ; height: 400px"
  .. .. .. .. .. .. .. ..$ children: list()
  .. .. .. .. .. .. .. ..- attr(*, "class")= chr "shiny.tag"
  .. .. .. .. .. .. ..$ :List of 3
  .. .. .. .. .. .. .. ..$ name    : chr "div"
  .. .. .. .. .. .. .. ..$ attribs :List of 3
  .. .. .. .. .. .. .. .. ..$ id   : chr "GenomeMelt"
  .. .. .. .. .. .. .. .. ..$ class: chr "shiny-plot-output"
  .. .. .. .. .. .. .. .. ..$ style: chr "width: 100% ; height: 400px"
  .. .. .. .. .. .. .. ..$ children: list()
  .. .. .. .. .. .. .. ..- attr(*, "class")= chr "shiny.tag"
  .. .. .. .. .. .. ..$ :List of 3
  .. .. .. .. .. .. .. ..$ name    : chr "pre"
  .. .. .. .. .. .. .. ..$ attribs :List of 2
  .. .. .. .. .. .. .. .. ..$ id   : chr "nText"
  .. .. .. .. .. .. .. .. ..$ class: chr "shiny-text-output noplaceholder"
  .. .. .. .. .. .. .. ..$ children: list()
  .. .. .. .. .. .. .. ..- attr(*, "class")= chr "shiny.tag"
  .. .. .. .. .. ..- attr(*, "class")= chr "shiny.tag"
  .. .. .. ..- attr(*, "class")= chr "shiny.tag"
  .. ..- attr(*, "class")= chr "shiny.tag"
 - attr(*, "class")= chr [1:2] "shiny.tag.list" "list"
 - attr(*, "html_dependencies")=List of 1
  ..$ :List of 10
  .. ..$ name      : chr "bootstrap"
  .. ..$ version   : chr "3.4.1"
  .. ..$ src       :List of 2
  .. .. ..$ href: chr "shared/bootstrap"
  .. .. ..$ file: chr "C:/Users/joelt/Desktop/RStudio/R-3.5.3/library/shiny/www/shared/bootstrap"
  .. ..$ meta      :List of 1
  .. .. ..$ viewport: chr "width=device-width, initial-scale=1"
  .. ..$ script    : chr [1:3] "js/bootstrap.min.js" "shim/html5shiv.min.js" "shim/respond.min.js"
  .. ..$ stylesheet: chr "css/bootstrap.min.css"
  .. ..$ head      : NULL
  .. ..$ attachment: NULL
  .. ..$ package   : NULL
  .. ..$ all_files : logi TRUE
  .. ..- attr(*, "class")= chr "html_dependency"
}
\details{
%%  ~~ If necessary, more details than the __description__ above ~~
}
\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
}
\references{
%%  ~~ possibly secondary sources and usages ~~
}
\examples{
data(ui)
## maybe str(ui) ; plot(ui) ...
}
\keyword{datasets}
