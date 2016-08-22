## loading library
print("Loading required R packages...")
tryCatch({
    library(data.table)
    library(ROCR)
    library(ggplot2)
    library(rPython)
    ## library(UniProt.ws)
    library(gProfileR)
    library(rJava)
    library(RWeka)
    WPM('install-package','alternatingDecisionTrees')
    WPM('load-package','alternatingDecisionTrees')
    LADTree = make_Weka_classifier('weka/classifiers/trees/LADTree')
}, error=function(e) {
    print("At least one required package failed to load. Make sure you have all the depencies.")
    stop(conditionMessage(e))
    }
)

########## MAIN ##########
## part 1: gather all data
## task 1.a: establish reference proteome
print("")

## decide which organism
taxId = readline(prompt = "Type the taxonomical ID of the organism:")

uniprotTab = paste("data/", taxId, ".tab", sep="")

if (file.exists(uniprotTab)){
    proteome = fread(uniprotTab)
} else {
    stop("Download Uniprot ref proteome in tab-delimited file first!")
}


## Part 2: fitting LADtree

## part 3: visualize results
