## loading library
print("Loading required R packages...")
tryCatch({
    library(data.table)
    library(reshape2)
    library(RMySQL); mysql = dbDriver(drvName = "MySQL")
    library(ROCR)
    library(ggplot2)
    library(rPython)
    library(UniProt.ws)
    library(gProfileR)
    library(rJava)
    library(RWeka)
    WPM('install-package','alternatingDecisionTrees')
    WPM('load-package','alternatingDecisionTrees')
    LADTree = make_Weka_classifier('weka/classifiers/trees/LADTree')
    library(h2o)
}, error=function(e) {
    print("At least one required package failed to load.")
    print("Make sure you have installed all the depencies.")
    stop(conditionMessage(e))
    }
)

########## MAIN ##########
## part 1: gather all data
## task 1: establish reference proteome

## decide which organism
taxId = readline(prompt = "Type the taxonomical ID of the organism: ")
orgMap = setNames(c("hsapiens"), c("9606"))

## download ref proteome from Uniprot
uniprotTab = paste("data/", taxId, ".tab", sep="")

if (file.exists(uniprotTab)){
    up = UniProt.ws(taxId = as.numeric(taxId))
    proteome = fread(uniprotTab)
    
    ## munging: unify column names
    if (taxId == "9606") {
        names(proteome) = gsub(pattern = "Entry",
                               replacement = "uniprotKb", x = names(proteome))
        names(proteome) = gsub(pattern = "Cross-reference",
                               replacement = "", x = names(proteome))
        names(proteome) = gsub(pattern = " (Ensembl)",
                               replacement = "ensembl",
                               x = names(proteome), fixed = T)
        names(proteome) = gsub(pattern = " (GeneID)",
                               replacement = "geneId",
                               x = names(proteome), fixed = T)
        names(proteome) = gsub(pattern = "Protein names",
                               replacement = "protein", x = names(proteome))
        names(proteome) = gsub(pattern = "Gene names  (primary )",
                               replacement = "symbol",
                               x = names(proteome), fixed = T)
        names(proteome) = gsub(pattern = "Gene names",
                               replacement = "gene", x = names(proteome))
        
        proteomeAll = proteome
        ## filter1: full annotation score
        proteome = proteome[Annotation=="5 out of 5"]
        ## convert numeric value
        proteome[, Mass := as.numeric(gsub(",","",Mass))]
        proteome[, Length := as.numeric(Length)]

        ## filter2: non-empty HGNC symbol
        proteinWoSymbol = proteome[symbol=="", protein]
        saveRDS(proteinWoSymbol,
                paste("data/", taxId, ".proteinWoSymbol.rds", sep=""))
        proteome = proteome[symbol != ""]
        
        ## filter3: deduplicate by hgnc, keep the ones with geneId
        symbolDup = proteome[which(duplicated(symbol)), unique(symbol)]
        proteomeSymbolDedup = 
            proteome[symbol %in% symbolDup &
                         geneId != "" &
                         !duplicated(geneId)][!duplicated(symbol)]
        proteome = rbind(proteome[!(hgnc %in% hgncDup)], proteomeHgncDedup)
        
        ## save the data
        setkey(proteome, "uniprotKb")
        write.table(proteome, 
                    gsub(".tab", ".mod.tab", uniprotTab),
                    quote=F, sep = "\t", row.names = F)
        saveRDS(proteome,
                paste("./data/",
                      paste(taxId, "proteome", "rds", sep="."), sep = ""))
    } else if (taxId == "559292") {
        print("analyzing S cerevisiae (S288c) data.")
    }
} else {
    stop("Download Uniprot ref proteome in tab-delimited file first!")
}

## task 2: assemble training labels of SUMO substrates
sumo = fread(paste("data/", taxId, ".sumo.txt", sep = ""))
if (taxId=="9606") {

}





## task 3: use gProfileR to find significant terms
gprofiler()

## task 4: retrieve GO-gene association for all genes and all selected terms
## set up connection
# goConn = dbConnect(mysql, user='go_select', password='amigo',
#                    host='mysql-amigo.ebi.ac.uk', database='go_latest',
#                    port=4085)
goConn2 = dbConnect(mysql, user='go_select', password='',
                   host='spitz.lbl.gov', dbname='go_latest')

## task 3: retrieve STRING database
## download protein links file from website
stringDb = fread("data/9606.protein.links.v10.txt")
## mapping from ENSP to UniprotKb
ensp = 

## task 4: pre-compute GPS-SUMO
seq = select(x = up, columns = c("UNIPROTKB", "SEQUENCE"),
             keytype = "UNIPROTKB", keys = proteome[, uniprotKb])
seq = data.table(seq)
writeLines(seq[, paste(">", UNIPROTKB, "\n", SEQUENCE, "\n", sep="")], con = paste("data/", taxId, ".seq.fa", sep=""), sep = "\n")

## Manual step: feed the sequences into GPS-SUMO 2.0 web server
## download the result as text and rename it ./data/[taxId].gps.txt
gps = fread(paste("data/",taxId,".gps.txt", sep=""))
gpsCount = data.table(dcast(gps, ID ~ Type))


## task 5: retrieve labels

## Part 2: fitting LADtree

## part 3: visualize results
