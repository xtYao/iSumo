## TODO: in the final release, remove these two lines
setwd("gitLcl/iSumo/")
load(".RData")
## loading library
print("Loading required R packages...")
tryCatch({
    library(data.table)
    library(reshape2)
    library(RMySQL)
    library(ROCR)
    library(ggplot2)
    library(rPython)
    library(UniProt.ws)
    library(gProfileR)
    library(h2o)
    library(knncat)
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
orgMap = setNames(c("hsapiens","scerevisiae"),
				  c("9606","559292"))

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
    setkey(sumo, uniprotKb)
    sumo[, hits := rowSums(sumo[, -("uniprotKb"), with=F])]
    sumo[, isSumo := hits>0]
} else if (taxId == "559292") {
    print("gathering yeast sumo labels.")
}

## task 3: use gProfileR to find significant terms
enrich = gprofiler(query = sumo[isSumo==T, uniprotKb],
                   organism = orgMap[taxId], ordered_query = F,
                   exclude_iea = T, custom_bg = proteome$uniprotKb,
                   significant = F)
enrich = data.table(enrich)
if(taxId == "9606"){
	sigGo = enrich[domain %in% c("MF", "CC", "BP") & 
				   significant==T & term.size<2000 &
				   	!grepl("sumo", term.name),
				   .(term.id, term.name, term.size, p.value)]
} else if (taxId == "599292") {
	print("yeast GO term enrichment.")
}
setkey(sigGo, term.name)

## task 4: retrieve GO-gene association for all genes and all selected terms
## set up connection
getUniprotKbByTermId = function(term.id, taxId="9606"){
    ## Note: this function is vectorized in respect to term.id
    ## TODO: check if taxId is of length 1
    
    ## set up GO MySQL connection
    mysql = dbDriver("MySQL")
    goConn2 = dbConnect(mysql, user='go_select', password='',
                        host='spitz.lbl.gov', dbname='go_latest')

    ## set up query
    stmtUniprotKbByTermId = sprintf(
        "SELECT DISTINCT dbxref.xref_key AS uniprotKb
        FROM term
        INNER JOIN graph_path ON (term.id = graph_path.term1_id)
        INNER JOIN association ON (graph_path.term2_id = association.term_id)
        INNER JOIN gene_product ON (association.gene_product_id = gene_product.id)
        INNER JOIN species ON (gene_product.species_id = species.id)
        INNER JOIN dbxref ON (gene_product.dbxref_id = dbxref.id)
        WHERE
        acc = '%s'
        AND
        ncbi_taxa_id = '%s'
        AND
        dbxref.xref_dbname = 'UniprotKB'",
        term.id, taxId)
    
    ## execute the query and return a list
    ## of the same length of term.id, each element is a char vec of UniprotKb
    res = lapply(stmtUniprotKbByTermId,
             function(x){
                 dbGetQuery(goConn2, x)$uniprotKb
             })
    res = setNames(res, term.id)
    dbDisconnect(goConn2)
    return(res)
}

## expand proteome with GO assocaition
goMat = as.data.table(lapply(getUniprotKbByTermId(sigGo$term.id, taxId), 
                             function(x) proteome$uniprotKb %in% x))
goMat = proteome[, .(uniprotKb)]
for (id in sigGo$term.id){
    res = getUniprotKbByTermId(id, taxId)
    goMat[[sigGo[term.id==id, term.name]]] = goMat$uniprotKb %in% res[[id]]
}
write.table(goMat, paste("data/", taxId, ".goMat.txt", sep = ""),
            quote = F, row.names = F, sep = "\t")

## task 5: retrieve STRING database
## download protein links file from website
protInteract = fread("data/9606.protein.actions.v10.txt")
protInteract = protInteract[mode=="binding",
							.(p1=item_id_a, p2=item_id_b)]
if (any(grepl("\\.", protInteract$p1))){
	protInteract$p1 = sapply(strsplit(protInteract$p1, split = "\\."),
							 function(x) x[2])
	protInteract$p2 = sapply(strsplit(protInteract$p2, split = "\\."),
							 function(x) x[2])
}
## mapping from ENSP to UniprotKb
ensp = fread(paste("data/", taxId, ".uniprotKb2ensp.txt", sep=""))
setkey(ensp, ensp)
protInteract$u1 = ensp[protInteract$p1, uniprotKb]
protInteract$u2 = ensp[protInteract$p2, uniprotKb]
protInteract = protInteract[!is.na(u1) & !is.na(u2), .(u1, u2)]
write.table(protInteract, paste("data/", taxId, ".stringInt.txt", sep=""),
			quote = F, sep = "\t", row.names = F)

## TODO: calculate degrees or find if they are recorded somewhere
ppiDegree = sapply(proteome$uniprotKb,
				   function(x){
				   	nrow(protInteract[u1==x | u2==x,])
				   })

## task 6: retrieve CORUM database, only if for human data
## summarize if a protein is within a complex
## AND how large is that complex
if (taxId == "9606"){
	corum = fread("data/coreCORUM.txt", sep = ";")
	corum = corum[organism == "Human", c(2, 5), with=F]
	corumSubunits = setNames(
		sapply(corum$`subunits (UniProt IDs)`, strsplit, ","),
		corum$`Complex name`)
	humanComplex = t(sapply(proteome$uniprotKb,
						  function(x){
						  	inCorum = which(sapply(corumSubunits,
						  						   function(y) x %in% y))
						  	nCorum = length(inCorum) # n comp has the prot
						  	avgCorumSz = 0
						  	if (nCorum != 0){
						  		avgCorumSz = # avg n subunit per comp
						  		mean(sapply(corumSubunits[inCorum], length))
						  	}
						  	return(c(nCorum, avgCorumSz))
						  }))
	humanComplex = as.data.table(humanComplex)
	colnames(humanComplex) = c("nCorum", "avgCorumSz")
}

## task 7: phosphosite
## include number of P/M/A/U modification sites
if (taxId == "9606"){
	## read in phosphosite data
	phos = fread("data/psp/Phosphorylation_site_dataset", skip = 3)
	meth = fread("data/psp/Methylation_site_dataset", skip=3)
	acet = fread("data/psp/Acetylation_site_dataset", skip=3)
	ubiq = fread("data/psp/Ubiquitination_site_dataset", skip = 3)
	
	## count how many each protein has
	phosCount = phos[, table(ACC_ID)]
	methCount = meth[, table(ACC_ID)]
	acetCount = acet[, table(ACC_ID)]
	ubiqCount = ubiq[, table(ACC_ID)]
	
	## construct feature vectors
	ptmMat = proteome[, .(uniprotKb)]
	ptmMat[, nPhos := 
		   	ifelse(uniprotKb %in% names(phosCount), phosCount[uniprotKb], 0)]
	ptmMat[, nMeth := 
		   	ifelse(uniprotKb %in% names(methCount), methCount[uniprotKb], 0)]
	ptmMat[, nAcet := 
		   	ifelse(uniprotKb %in% names(acetCount), acetCount[uniprotKb], 0)]
	ptmMat[, nUbiq := 
		   	ifelse(uniprotKb %in% names(ubiqCount), ubiqCount[uniprotKb], 0)]
}

## task 8: pre-compute GPS-SUMO
seq = select(x = up, columns = c("UNIPROTKB", "SEQUENCE"),
             keytype = "UNIPROTKB", keys = proteome[, uniprotKb])
seq = data.table(seq)
writeLines(seq[, paste(">", UNIPROTKB, "\n", SEQUENCE, "\n", sep="")],
           con = paste("data/", taxId, ".seq.fa", sep=""), sep = "\n")

## Manual step: feed the sequences into GPS-SUMO 2.0 web server
## download the result as text and rename it ./data/[taxId].gps.txt
gps = fread(paste("data/",taxId,".gps.txt", sep=""))
gpsCount = data.table(dcast(gps, ID ~ Type)); setkey(gpsCount, ID)
gpsSumo = do.call(rbind, 
				  lapply(proteome$uniprotKb,
				  	   function(x){
				  	   	if (x %in% gpsCount$ID){
				  	   		return(gpsCount[x, -1, with=F])
				  	   	} else {
				  	   		res = setNames(c(0,0,0), colnames(gpsCount)[-1])
				  	   		return(as.data.table(as.list(res)))
				  	   	}
				  	   }))

## finally, put together the full dataset
## TODO: hardcoded for human, make it generic
## TODO: rerun
iSumoData = cbind(proteome[, .(uniprotKb, symbol, Length, Mass)],
				  goMat[,-1,with=F], ## GO annotation
				  ptmMat[,-1,with=F],
				  ppiDegree, humanComplex, 
				  gpsSumo, sumo[, .(isSumo)])
write.table(iSumoData, paste("data/",taxId,".iSumo.txt",sep=""),
			quote = F, row.names = F, sep = "\t")

## Part 2: fitting RF


## part 3: visualize results
