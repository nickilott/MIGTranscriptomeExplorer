#' Run the shiny app
#'
#' Run the shiny app
#' @import shiny
#' @examples
#' runExplorer()
#' @export

runExplorer <- function(){

    app <- system.file("exec/app.R", package="MIGTranscriptomeExplorer")
    runApp(app)
}