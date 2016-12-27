taxId = "559292"
proteome = readRDS("data/559292.proteome.rds")
yeastSumoRaw = fread("data/559292.sumo.raw.txt")

write.table(proteome$sgd, "data/559292.sgd.txt", quote = F, row.names = F)

## systematic name -- sgd mapping
sysName = read.table("data/sysNameSgdMapping.csv",
					 sep = " ", quote = '"', header = T)
sysName = as.data.table(sysName)[reason=="MATCH"]
setkey(sysName, "secondaryIdentifier")

dup = yeastSumoRaw[duplicated(systematic_name),systematic_name]
yeastSumoRaw[systematic_name %in% dup]
yeastSumoRaw[, sgd := sysName[systematic_name, primaryIdentifier]]

yeastSumoSys = yeastSumoRaw[systematic_name != "not_found",systematic_name]

write.table(yeastSumoRaw$systematic_name, "data/559292.sumo.raw.sysName.txt",
			quote=F, sep="\t", row.names = F)
saveRDS(unique(yeastSumoRaw$sgd), "data/559292.sumoSgd.rds")

yeastSumo = proteome[, 1:5, with=F]
yeastSumo[, ":="(
	`Wykoff and O'Shea, 2005`=
		ifelse(sgd %in% yeastSumoRaw[Wykoff==1, sgd], TRUE, FALSE),
	`Hannich et al., 2005`=
		ifelse(sgd %in% yeastSumoRaw[Hannich_MS==1, sgd], TRUE, FALSE),
	`Denison et al., 2005`=
		ifelse(sgd %in% yeastSumoRaw[Denison==1, sgd], TRUE, FALSE),
	`Zhou et al., 2004`=
		ifelse(sgd %in% yeastSumoRaw[Zhou==1, sgd], TRUE, FALSE),
	`Wolschlegel et al., 2004`=
		ifelse(sgd %in% yeastSumoRaw[Wohlschlegel==1, sgd], TRUE, FALSE),
	`Panse et al., 2004`=
		ifelse(sgd %in% yeastSumoRaw[Panse==1, sgd], TRUE, FALSE)
	)]

write.table(yeastSumo[,c(6:11,1),with=F], "data/559292.sumo.txt",
			quote = F, row.names = F, sep = "\t")
