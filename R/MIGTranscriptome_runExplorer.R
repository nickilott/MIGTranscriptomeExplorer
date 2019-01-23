#' Run the shiny app
#'
#' Run the shiny app
#' @import shiny
#' @examples
#' runExplorer()
#' @export

runExplorer <- function(){

    app <- system.file("R/app.R", package="MIGTranscriptomeExplorer")
    shiny(app)
}