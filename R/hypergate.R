{
    library(knitr)
    library(rmarkdown)
    if(Sys.info()['sysname']=="Linux"){
        pkg.path="/media/etienne/Data/FCtools/Packages/hypergate.dev/hypergate/vignettes/"
    } else {
        pkg.path="D:/FCtools/Packages/hypergate.dev/hypergate/vignettes"
    }
    setwd(pkg.path)
    ##knitr::opts_chunk$set(fig.pos = 'h',tidy.opts=list(width.cutoff=60),tidy=TRUE)
    knit("./rmd/hypergate_vignette.Rmd")
    rmarkdown::render("./hypergate_vignette.md")
    rmarkdown::render("./hypergate_vignette.md",output_format="pdf_document",output_file="hypergate_vignette.pdf")
}
