## TODO: in the final release, remove these two lines
setwd("~/gitLcl/iSumo/")
load(".RData")
## loading library
print("Loading required R packages...")
tryCatch({
    library(data.table)
    library(reshape2)
    library(RMySQL)
    library(ROCR)
    library(ggplot2)
    library(UniProt.ws) ## bioc
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
orgMap = setNames(c("hsapiens","scerevisiae"), c("9606","559292"))

proteomeFn =
    paste("./data/", paste(taxId, "proteome", "rds", sep="."), sep = "")
if (!file.exists(proteomeFn)){
    ## download ref proteome from Uniprot
    uniprotTab = paste("data/", taxId, ".tab", sep="")
    proteome = fread(uniprotTab)

    if (file.exists(uniprotTab)){
        up = UniProt.ws(taxId = as.numeric(taxId))


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
        } else if (taxId == "559292") {
            print("analyzing S cerevisiae (S288c) data.")

            ##
            names(proteome) = gsub(pattern = "Entry",
                                   replacement = "uniprotKb", x = names(proteome))
            names(proteome) = gsub(pattern = "Protein names",
                                   replacement = "protein", x = names(proteome))
            names(proteome) = gsub(pattern = "Cross-reference",
                                   replacement = "", x = names(proteome))
            names(proteome) = gsub(pattern = " (GeneID)",
                                   replacement = "geneId",
                                   x = names(proteome), fixed = T)
            names(proteome) = gsub(pattern = " (SGD)",
                                   replacement = "sgd",
                                   x = names(proteome), fixed = T)

            ## use primary gene name, if empty use first item in gene names
            proteome[, gene := `Gene names  (primary )`]
            proteome[, gene := ifelse(nchar(gene)!=0 & !grepl(";", gene), gene,
                                      sapply(strsplit(proteome$`Gene names`[1:10], " "), head, 1))]

            ## remove semi colons from sgd and geneId
            proteome[, geneId := gsub(";", "", x = geneId)]
            proteome[, sgd := gsub(";", "", x = sgd)]

            ## convert Mass and Annotation to numeric
            proteome[, Mass := as.numeric(gsub(",","",Mass))]
            proteome[, Annotation :=
                         as.numeric( sapply(strsplit(Annotation, split = " "), head, 1) )]

            ## get rid of the duplicated sgd
            proteome = proteome[!duplicated(sgd), .(uniprotKb, sgd, gene,
                                                    protein, geneNames=`Gene names`,
                                                    Length, Mass,
                                                    Status, Annotation)]
            setkey(proteome, "uniprotKb")
        }
    } else {
        stop("Download Uniprot ref proteome in tab-delimited file first!")
    }
    saveRDS(proteome, proteomFn)
} else {
    proteome = readRDS(proteomeFn)
}


## task 2: assemble training labels of SUMO substrates
sumo = fread(paste("data/", taxId, ".sumo.txt", sep = ""))
setkey(sumo, uniprotKb)
sumo[, hits := rowSums(sumo[, -("uniprotKb"), with=F])]
sumo[, isSumo := hits>0]


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
	sigGo = enrich[domain %in% c("MF", "CC", "BP") &
				   	term.size>5 & term.size<1000 &
				   	significant==T & !grepl("sumo", term.name),
				   .(term.id, term.name, term.size, p.value)]
}
setkey(sigGo, term.name)
## Cleaning: don't include any annotation directly show SUMO status
sigGo = sigGo[!grepl("sumo", term.name, ignore.case = T)]

## Table 1: significantly enrich GO terms in SUMO set
write.table(sigGo[order(p.value)], paste(taxId,".sigGo.txt",sep = ""),
			sep = "\t", row.names = F, quote = F)

## task 4: retrieve GO-gene association for all genes and all selected terms
## set up connection
getUniprotKbByTermId = function(term.id, taxId="9606"){
    ## Note: this function is vectorized in respect to term.id
    ## TODO: check if taxId is of length 1

    ## set up GO MySQL connection
    mysql = dbDriver("MySQL")
#    goConn2 = dbConnect(mysql, user='go_select', password='',
 #                       host='spitz.lbl.gov', dbname='go_latest')
	goConn2 = dbConnect(mysql, user='go_select', password='amigo',
					   host='mysql-amigo.ebi.ac.uk', port=4085,
					   dbname='go_latest')

    ## set up query
	if (taxId=="9606"){
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
	} else if (taxId=="559292"){
		stmtUniprotKbByTermId = sprintf(
			"SELECT DISTINCT dbxref.xref_key AS sgd
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
			dbxref.xref_dbname = 'SGD'",
			term.id, taxId)
	} else {
		stop("Not yet implemented for this organism!!!")
	}

    ## execute the query and return a list
    ## of the same length of term.id, each element is a char vec of UniprotKb
    res = lapply(stmtUniprotKbByTermId,
             function(x){
             	if (taxId=="9606"){
             		dbGetQuery(goConn2, x)$uniprotKb
             	} else if (taxId == "559292") {
             		dbGetQuery(goConn2, x)$sgd
             	}
             })
    res = setNames(res, term.id)
    dbDisconnect(goConn2)
    return(res)
}

## expand proteome with GO assocaition
goMat = as.data.table(lapply(sigGo$term.id, function(x){
    ids = unlist(getUniprotKbByTermId(x, taxId = taxId), use.names = F)
    if (taxId == "9606"){
        proteome$uniprotKb %in% ids
    } else if (taxId == "559292"){
        proteome$sgd %in% ids
    }
}))
colnames(goMat) = sigGo[, term.name]

write.table(goMat, paste("data/", taxId, ".goMat.txt", sep = ""),
            quote = F, row.names = F, sep = "\t")

## task 5: retrieve STRING database
## download protein links file from website
if (taxId=="9606"){
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

} else if (taxId=="559292") {
	## NOTE: STRING db doesn't have data for 559292 (S288c), but only 4932 (S.
	## cerevisiae). We will use that instead.
	protInteract = fread("data/4932.protein.actions.v10.txt")
	protInteract = protInteract[mode=="binding",
								.(p1=item_id_a, p2=item_id_b)]
	if (any(grepl("\\.", protInteract$p1))){
		protInteract$p1 = sapply(strsplit(protInteract$p1, split = "\\."),
								 function(x) x[2])
		protInteract$p2 = sapply(strsplit(protInteract$p2, split = "\\."),
								 function(x) x[2])
	}

	## mapping of ORF name and sgd
	sysName = read.table("data/sysNameSgdMapping.csv",
						 sep = " ", quote = '"', header = T)
	sysName = as.data.table(sysName)[reason=="MATCH"]
	setkey(sysName, "secondaryIdentifier")

	protInteract$s1 = sysName[protInteract$p1, as.character(primaryIdentifier)]
	protInteract$s2 = sysName[protInteract$p2, as.character(primaryIdentifier)]
	sgd2uniprotKb = do.call("rbind",
							lapply(unique(c(protInteract$s1, protInteract$s2)),
				function(x) proteome[grepl(x, sgd), .(sgd, uniprotKb)]))
	sgd2uniprotKb = sgd2uniprotKb[!duplicated(sgd)]
	setkey(sgd2uniprotKb, "sgd")
	protInteract$u1 = sgd2uniprotKb[protInteract$s1, uniprotKb]
	protInteract$u2 = sgd2uniprotKb[protInteract$s2, uniprotKb]
	protInteract = protInteract[!is.na(u1) & !is.na(u2), .(u1, u2)]
	write.table(protInteract, paste("data/", taxId, ".stringInt.txt", sep=""),
				quote = F, sep = "\t", row.names = F)
}
## calculate degrees
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
} else if (taxId == "559292"){
	## use the same dataset in A Baryshnikova 2016 yeast genetic interaction ppr
	complexes = fread("data/559292.complexes.txt", sep = "\t")
	## a list of complex compositions
	cpxSubunits = setNames(
		sapply(complexes$`ORFs annotated to complex`, strsplit, "; "),
		complexes$`Protein Complex Name`)

	## change key for sysName
	sysName$primaryIdentifier = as.character(sysName$primaryIdentifier)
	sysName$secondaryIdentifier = as.character(sysName$secondaryIdentifier)
	setkey(sysName, "primaryIdentifier")

	## convert into complex memberships
	yeastComplex =
		t(sapply(proteome$uniprotKb,
				function(x){
					xOrf = sysName[proteome[x, sgd], secondaryIdentifier]
					inComp = which(sapply(cpxSubunits, function(y) xOrf %in% y))
					nComp = length(inComp) # n comp has the prot
					avgCompSz = 0
					if (nComp != 0){
						avgCompSz = # avg n subunit per comp
						mean(sapply(cpxSubunits[inComp], length))
					}
					return(c(nComp, avgCompSz))
				}))
	yeastComplex = as.data.table(yeastComplex)
	colnames(yeastComplex) = c("nComp", "avgCompSz")
}

## Table 2/Figure 2. analyzing RNA-binding, SUMO in protein complexes
## get proteins that are RNA-binding
humanRnaBinding = unlist(getUniprotKbByTermId(term.id = "GO:0003723"),
						 use.names = F)
yeastRnaBinding = unlist(getUniprotKbByTermId(term.id = "GO:0003723",
											  taxId = "559292"),
						 use.name = F)
yeastRnaBinding = sysName[primaryIdentifier %in% yeastRnaBinding,
						  secondaryIdentifier]
humanSumo = Reduce(union, readRDS("data/9606.sumoUniprotKb.rds"))
yeastSumo = readRDS("data/559292.sumoSgd.rds")

## assemble human and yeast complex centric data separately
hDt = do.call("rbind", lapply(names(corumSubunits), function(x){
	prots = corumSubunits[[x]]
	size = length(prots)
	nRna = sum(prots %in% humanRnaBinding)
	nSumo = sum(prots %in% humanSumo)
	return(c(x, size, "human", nRna, nSumo))
}))
yDt = do.call("rbind", lapply(names(cpxSubunits), function(x){
	prots = cpxSubunits[[x]]
	size = length(prots)
	nRna = sum(prots %in% yeastRnaBinding)
	nSumo = sum(prots %in% yeastSumoSys)
	return(c(x, size, "yeast", nRna, nSumo))
}))
## put them together into one dt, rename, convert classes
complexComp = data.table(rbind(hDt, yDt))
colnames(complexComp) = c("Complex name", "Complex size", "Organism",
						  "Number of RNA-binding subunits",
						  "Number of SUMOylated subunits")
class(complexComp$`Complex size`)="numeric"
class(complexComp$`Number of RNA-binding subunits`)="numeric"
class(complexComp$`Number of SUMOylated subunits`)="numeric"
complexComp$Organism = as.factor(complexComp$Organism)
complexComp[, ":="("RNA binding" = `Number of RNA-binding subunits`>0,
				   "SUMOylated" = `Number of SUMOylated subunits`>0)][,
				   "RNA.SUMO" := interaction(`RNA binding`, `SUMOylated`)]
## save it.
write.table(complexComp, "tableS3.complexes.txt",
			sep = "\t", row.names = F, quote = F)


## analysis of complex size correlation with RNA-binding or SUMO
## TODO: annotate the graph with the following!!!
complexComp[Organism=="human", table(`RNA binding`, SUMOylated)]
complexComp[Organism=="yeast", table(`RNA binding`, SUMOylated)]
## pairwise wilcox rank test
pairwise.wilcox.test(x = complexComp[Organism=="human", `Complex size`],
					 g = complexComp[Organism=="human", `RNA.SUMO`],
					 p.adjust.method = "fdr")
pairwise.wilcox.test(x = complexComp[Organism=="yeast", `Complex size`],
					 g = complexComp[Organism=="yeast", `RNA.SUMO`],
					 p.adjust.method = "fdr")

## first make boxplot to show complex size related to RNA/SUMO
## (violin and jitter too busy since majority of points are at the bottom)
fig2a = ggplot(data = complexComp,
			   mapping = aes(x = RNA.SUMO, y = `Complex size`))
svg(width = 12, height = 8, filename = "fig2a.compSzRnaSumo.svg")
plot.new()
fig2a + geom_boxplot() +
#	geom_violin(draw_quantiles = T) +
#	geom_jitter(width = 0.2) +
	facet_wrap(facets = ~Organism) +
#	annotate("text", )
	theme(plot.background = element_blank(),
		  panel.grid = element_blank(),
		  panel.background = element_blank(),
		  plot.margin = margin(t = 20, r = 60, b = 20, l = 20, unit = "pt"),
		  strip.placement = "indside",
		  text = element_text(size = 16),
		  axis.line =
		  	element_line(size=1, colour = "black",
		  				 arrow = arrow(angle = 30,
		  				 			  length = unit(x = 0.1, units = "inch"),
		  				 			  type = "open", ends = "last")),
#		  axis.title.x = element_text(size = 16),
#		  axis.title.y = element_text(size = 16),
		  axis.text.x = element_text(size = 16, angle = 335, hjust = 0),
#		  axis.text.y = element_text(size = 16),
		  axis.ticks = element_blank())
dev.off()

## second make scatter plot of nSUMO-nRNA, with size correspond to complex size
fig2b = ggplot(data = complexComp,
			   mapping = aes(x = `Number of RNA-binding subunits`,
			   			  y = `Number of SUMOylated subunits`,
			   			  color = Organism))
svg(width = 8, height = 8, filename = "fig2b.sumoRna.svg")
plot.new()
fig2b + geom_point(aes(size = complexComp$`Complex size`), alpha=0.5) +
	scale_size(breaks = c(5, 10, 20, 60, 100), range = c(1,8),
			   guide = guide_legend(title = "Complex size")) +
	## labeling human big SUMOs
	geom_text(mapping = aes(label = `Complex name`, size = 24),
		data = complexComp[Organism=="human"][
			order(`Number of SUMOylated subunits`,decreasing = T)][1:6],
		check_overlap = T, nudge_y = 3, nudge_x = 8, show.legend = F) +
	geom_text(mapping = aes(label = `Complex name`, size = 24),
			  data = complexComp[Organism=="human"][
			  	order(`Number of SUMOylated subunits`,decreasing = T)][7],
			  check_overlap = T, nudge_y = 3, nudge_x = 2, show.legend = F) +
	## labeling yeast big SUMOs
	geom_text(mapping = aes(label = `Complex name`, size = 24),
			  data = complexComp[Organism=="yeast"][
			  	order(`Number of SUMOylated subunits`,decreasing = T)][1],
			  check_overlap = T, nudge_y = 3, nudge_x = 16, show.legend = F) +
	geom_text(mapping = aes(label = `Complex name`, size = 24),
			  data = complexComp[Organism=="yeast"][
			  	order(`Number of SUMOylated subunits`,decreasing = T)][3],
			  check_overlap = T, nudge_y = 3, nudge_x = 2, show.legend = F) +
	geom_text(mapping = aes(label = `Complex name`, size = 32),
			  data = complexComp[Organism=="yeast"][
			  	order(`Number of SUMOylated subunits`,decreasing = T)][c(2,4)],
			  check_overlap = T, nudge_y = -3, nudge_x = 16, show.legend = F) +
	## labeling human big RNA binding low SUMOs
	geom_text(mapping = aes(label = `Complex name`, size = 32),
			  data = complexComp[Organism=="human" &
			  				   	`Number of SUMOylated subunits`<10][
			  	order(`Number of RNA-binding subunits`,decreasing = T)][1:3],
			  check_overlap = T, nudge_y = 4, nudge_x = 8, show.legend = F) +
	xlim(-7,133) +
	theme(plot.background = element_blank(),
		  panel.grid = element_blank(),
		  panel.background = element_blank(),
		  plot.margin = margin(t = 20, r = 20, b = 20, l = 20, unit = "pt"),
		  text = element_text(size = 16),
		  axis.line =
		  	element_line(size=1, colour = "black",
		  				 arrow = arrow(angle = 30,
		  				 			  length = unit(x = 0.1, units = "inch"),
		  				 			  type = "open", ends = "last")),
		  axis.ticks = element_blank(),
		  legend.background = element_blank(),
		  legend.key = element_blank(),
		  legend.key.size = unit(16, "pt"),
		  legend.position = c(0.2, 0.7),
		  legend.box = "vertical",
		  legend.text = element_text(size = 16))
dev.off()

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

} else if (taxId == "559292"){
	## no phosphosite data, use dbPTM data instead
	dbPTM = fread("data/dbPTM3.txt")
	yeastPtm = dbPTM[V2 %in% proteome$uniprotKb]

	phosCount = yeastPtm[V8=="Phosphorylation", table(V2)]
	methCount = yeastPtm[V8=="Methylation", table(V2)]
	acetCount = yeastPtm[V8=="Acetylation", table(V2)]
	ubiqCount = yeastPtm[V8=="Ubiquitylation", table(V2)]

}
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

## task 8: pre-compute GPS-SUMO
seq = select(x = up, columns = c("UNIPROTKB", "SEQUENCE"),
             keytype = "UNIPROTKB", keys = proteome[, uniprotKb])
seq = data.table(seq)
writeLines(seq[, paste(">", UNIPROTKB, "\n", SEQUENCE, "\n", sep="")],
           con = paste("data/", taxId, ".seq.fa", sep=""), sep = "\n")

## Manual step: feed the sequences into GPS-SUMO 2.0 web server
## download the result as text and rename it ./data/[taxId].gps.txt
gps = fread(paste("data/",taxId,".gps.txt", sep=""))
gpsCount = data.table(dcast(gps, ID ~ Type, value.var = "Type"))[,1:4,with=F]
setkey(gpsCount, ID)
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
if (taxId=="9606"){
	iSumoData = cbind(proteome[, .(uniprotKb, symbol, Length, Mass)],
					  goMat, ## GO annotation
					  ptmMat[,-1,with=F],
					  ppiDegree, humanComplex,
					  gpsSumo, sumo[, .(isSumo)])
} else if (taxId=="559292"){
	iSumoData = cbind(proteome[, .(uniprotKb, sgd, Length, Mass)],
					  goMat[, , with=F],
					  ptmMat[, -1, with=F],
					  ppiDegree = ppiDegree[proteome$uniprotKb],
					  yeastComplex,
					  gpsSumo, sumo[, .(isSumo)])
}
write.table(iSumoData, paste("data/",taxId,".iSumo.txt",sep=""),
			quote = F, row.names = F, sep = "\t")

## Part 2: fitting RF
source(paste(taxId, ".modelTuning.R", sep=""))
fn = dir(paste("models/", taxId, ".rf.h2o", sep=""), full.names = T)
rf = h2o.loadModel(fn)
fn.null = dir(paste("models/", taxId, ".rf.null.h2o", sep=""))
rf.null = h2o.loadModel(fn.null)

## part 3: visualize results
####################
## Figure 1. Stacked bar of number of SUMO proteins found in each paper,
## decomposed into confirmed by X studies
if (taxId == "9606"){
    m1 = melt(data = sumo, id.vars = "hits", measure.vars = 1:18)
    m2 = m1[hits>0 & value==TRUE, table(hits, variable)]
    m3 = melt(m2, value.name = "nSumo")
    m3$variable <- reorder(x = m3$variable, X = m3$nSumo, FUN = sum)
    m3 = data.table(m3)
    m3[, hits := ifelse(hits>=9, "9+", as.character(hits))]
} else if (taxId == "559292"){
	m1 = melt(data = sumo, id.vars = "hits", measure.vars = 1:6)
	m2 = m1[hits>0 & value==TRUE, table(hits, variable)]
	m3 = melt(m2, value.name = "nSumo")
	m3$variable <- reorder(x = m3$variable, X = m3$nSumo, FUN = sum)
	m3$hits = as.factor(m3$hits)
}

## some species specific pars
if (taxId == "9606"){
    pl.mar = margin(t = 0, r = 120, b = 12, l = 4, unit = "pt")
    width = 11
} else if (taxId == "559292"){
    pl.mar = margin(t = 0, r = 100, b = 12, l = 4, unit = "pt")
    width = 10
}

svg(filename = paste(taxId, ".sumoStudies.svg", sep = ""),
    width = width, height = 8)
plot.new()
fig1 = ggplot(data = m3) +
	geom_bar(mapping = aes(x = variable, y = nSumo, fill = hits),
			 stat = "identity") +
	ylab("Number of SUMOylated proteins") +
	scale_fill_brewer(name="Number of studies\nin consensus", type = "seq",
					  palette = "BuGn", direction = -1)+
	theme(axis.line =
		  	element_line(size=1, colour = "black",
		  				 arrow = arrow(angle = 30,
		  				 			  length = unit(x = 0.1, units = "inch"),
		  				 			  type = "open", ends = "last")),
		  panel.grid.major = element_blank(),
		  panel.grid.minor = element_blank(),
		  panel.border = element_blank(),
		  panel.background = element_blank()) +
	theme(axis.title.x = element_blank(),
		  axis.title.y = element_text(size = 16),
		  axis.text.x = element_text(size = 16, angle = 335, hjust = 0),
		  axis.text.y = element_text(size = 16),
		  axis.ticks = element_blank()) +
	theme(legend.background = element_blank(),
	      legend.position = c(0.2, 0.5),
		  legend.text = element_text(size = 16),
		  legend.key.size = unit(0.5, "inch"),
		  legend.title = element_text(size = 16)) +
    theme(plot.margin = pl.mar)
print(fig1)
dev.off()

###################
## Figure2.
## ROC curve, for either model, plot 10-CV in grey, validation in dotted
## red/blue, and test in solid red/blue. Add guideline y=x.
plotRoc = function(predLabelList, col=grey(0.25), roc=T, cv=F, dashed=F, lwd=1){
	## make performance measure
	pred = prediction(predLabelList[[1]], predLabelList[[2]])
	if (roc) perf = performance(pred, 'tpr', 'fpr') else
		perf = performance(Pred, 'prec', 'rec')

	## set line type
	lType = ifelse(dashed, 2, 1)

	## plot
	if (!cv){
		if (dev.cur()==1) {
			plot(perf, col=col, lty=lType, lwd = lwd)
		} else {
			plot(perf, col=col, lty=lType, lwd = lwd, add=T)
		}
	} else {
		if (dev.cur()==1) plot(perf, col='grey', lwd=0.2) else
			plot(perf, col='grey', lwd=0.2, add=T)
		##plot(Perf, avg='vertical', col=grey(0.3), add=T)
	}
}

## plotting
## rf: cv --> valid --> test
getModelPred = function(x, train, valid, test){
	pred = list()

	## get CV models predictions
	cvFold = as.data.frame(
		h2o.cross_validation_fold_assignment(x))$fold_assignment
	cvPred = lapply(h2o.cross_validation_predictions(x),
					function(y) as.data.frame(y))
	for (i in 1:length(cvPred)){
		thisName = paste("cv",i,sep="_")
		thisPred = as.data.frame(cvPred[[i]])[cvFold==(i-1), 3]
		thisLabel = as.data.frame(train)[cvFold==(i-1), "isSumo"]
		pred[[thisName]] = list(thisPred, thisLabel)
	}

	## get valid predictions
	validPred = as.data.frame(h2o.predict(x, valid))[, 3]
	validLabel = as.data.frame(valid)[, "isSumo"]
	pred$valid = list(validPred, validLabel)

	## get test predictions
	testPred = as.data.frame(h2o.predict(x, test))[, 3]
	testLabel = as.data.frame(test)[, "isSumo"]
	pred$test = list(testPred, testLabel)

	return(pred)
}

## make ROC, and P-R curves
pred = getModelPred(rf, train, valid, test)
pred.null = getModelPred(rf.null, train, valid, test)

svg(filename=paste(taxId, ".ROC.svg", sep=""),
	width=8, height=8,
	pointsize=12)
plot.new()
par(cex=1, cex.lab=1, cex.main=1, tcl=-0.1)
for (pr in names(pred)) {
	thisPred = pred[[pr]]
	thisPred.null = pred.null[[pr]]
	if (grepl("cv",pr)){
		## grey line: CV
		plotRoc(thisPred, roc = T, cv = T)
		plotRoc(thisPred.null, roc = T, cv = T)
	} else if (pr=="valid"){
		## dashed line: valid
		plotRoc(thisPred, roc = T, cv = F, dashed = T, col = "red", lwd=2)
		plotRoc(thisPred.null, roc = T, cv = F, dashed = T, col = "blue", lwd=2)
	} else {
		## solid line: test
		plotRoc(thisPred, roc = T, cv = F, dashed = F, col = "red", lwd=3)
		plotRoc(thisPred.null, roc = T, cv = F, dashed = F, col = "blue", lwd=3)
	}
}
## some garnish
axis(side = 1, at = seq(0,5)/5, labels = seq(0,5)/5)
axis(side = 2, at = seq(0,5)/5, labels = seq(0,5)/5)
title(##main = "Receiver operating characteristics (ROC)
##of iSUMO model versus with only predicted sequence motif",
	  xlab = "False positive rate", ylab = "True positive rate")
abline(a = 0, b = 1, lty=3, lwd=0.5, col=grey(0.25))
legend(x = 0.7, y = 0.3,
	   fill = c("red","blue",rgb(0, 0, 0, 0),rgb(0, 0, 0, 0)),lty = c(0,0,2,1),
	   border = F, bty = "n", xjust = 0.5, y.intersp = 1.25,
	   legend = c("iSUMO","seq motif only","validation set","test set"))
dev.off()

#################################
## Figure 3.
## collect relative importance of features
varImp = as.data.table(rf@model$variable_importances)
## sort variable levels by importance
varImp$variable = reorder(varImp$variable, varImp$relative_importance)

svg(filename=paste(taxId, ".varImp.svg", sep=""),
	width=10, height=8,
	pointsize=12)
plot.new()

fig3 = ggplot(data = varImp[1:30]) +
	geom_bar(mapping = aes(x = variable, y = scaled_importance),
			 stat = "identity", width = 0.6) +
	ylab(label = "Relative importance") +
	xlab(label = "Variable name") +
	theme(axis.line =
		  	element_line(size=1, colour = "black",
		  				 arrow = arrow(angle = 30,
		  				 			  length = unit(x = 0.1, units = "inch"),
		  				 			  type = "open", ends = "last")),
		  panel.grid.major = element_blank(),
		  panel.grid.minor = element_blank(),
		  panel.border = element_blank(),
		  panel.background = element_blank()) +
	theme(axis.title.x = element_blank(),
		  axis.title.y = element_text(size = 16),
		  axis.text.x = element_text(size = 16),
		  axis.text.y = element_text(size = 16),
		  axis.ticks = element_blank()) +
	coord_flip() +
	theme(legend.position = c(0.2, 0.8),
		  legend.text = element_text(size = 16),
		  legend.key.size = unit(0.5, "inch"),
		  legend.title = element_text(size = 16))
print(fig3)
dev.off()

#################################
## Table
finalPerformance = h2o.performance(rf, dt)
## find the threshold maximizing F2 value
threshold = as.data.table(
	finalPerformance@metrics$max_criteria_and_metric_scores)[
		metric=="max f2", threshold]

## dt is the h2o frame used for model tuning
finalPrediction = as.data.table(h2o.predict(rf, newdata = dt))[, -1, with=F]
finalPrediction[, ":="(uniprotKb = proteome$uniprotKb,
					   protein = proteome$protein,
					   gene = proteome$gene,
					   isSumo = iSumoData$isSumo)]
finalPrediction$predict = finalPrediction$TRUE.>threshold
setkey(finalPrediction, "uniprotKb")
saveRDS(finalPrediction, paste(taxId, ".finalPrediction.rds", sep=""))
write.table(finalPrediction, paste(taxId, ".finalPrediction.csv", sep=""),
			sep = "\t", row.names = F, quote = F)

##### SAVE WORKSPCE
save.image(file = paste(taxId, "RData", sep = "."))
