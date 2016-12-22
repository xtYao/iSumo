## munging human
taxId = "9606"
sumo = fread(paste("data/", taxId, ".SumoData.txt"))
sumoEnsg = lapply(sumo[,-(1:2),with=F], function(x) sumo$ENSG[which(!is.na(x))])

names(ensg) = c("uniprotKb", "ensg")
write.csv(ensg, "data/9606.uniprotKb2ensg.txt",
          quote = F, row.names = F, sep = "\t")
sumoUniprotKb = lapply(sumoEnsg, function(x) ensg[ensg %in% x, uniprotKb])

## newly added 4 studies
hendriks = fread("data/hendriks.txt")
hendriks = hendriks[, unique(uniprotKb)]
hendriks2 = fread("data/hendriks2.txt", sep = "\t")
hendriks2 = sort(unique(do.call(c, strsplit(hendriks2$uniprotKb, split = ";"))))
hendriks2 = unique(gsub(hendriks2, pattern = "-*", replacement = ""))
impens = fread("data/Impens.txt")
impens = unique(impens$uniprotKb)
lamoliatte = fread("data/lamoliatte.txt")
lamoliatte = unique(lamoliatte$uniprotKb)
tammsalu = fread("data/tammsalu.txt")
tammsalu = unique(tammsalu$uniprotKb)

sumoUniprotKb$hendriks = hendriks
sumoUniprotKb$impens = impens
sumoUniprotKb$lamoliatte = lamoliatte
sumoUniprotKb$tammsalu = tammsalu

humanStudies = c("Lamoliatte et al., 2013",
                 "Grant, 2010",
                 "Matic et al., 2010",
                 "Blomster et al., 2010",
                 "Blomster et al., 2009",
                 "Manza et al., 2004",
                 "Vertegaal et al., 2006",
                 "Golebiowski et al., 2009",
                 "Bruderer et al., 2011",
                 "Galisson et al., 2011",
                 "Rosas-Acosta et al., 2005",
                 "Vertegaal et al., 2004",
                 "Tatham et al., 2011",
                 "Schimmel et al., 2008",
                 "Hendriks et al., 2014",
                 "Impens et al., 2014",
                 "Lamoliatte et al., 2014",
                 "Tammsalu et al., 2014")
sumoUniprotKb = setNames(sumoUniprotKb, humanStudies)

saveRDS(sumoUniprotKb, "9606.sumoUniprotKb.rds")

sumo = as.data.table(lapply(sumoUniprotKb, function(x) proteome$uniprotKb %in% x))
sumo[, uniprotKb := proteome$uniprotKb]
write.table(sumo, "data/9606.sumo.txt", quote = F, row.names = F, sep = "\t")
