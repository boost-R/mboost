
###############################################
###
###	Data preprocessing for data published in 
###
###     Use of Gene-Expression Profiling to Identify Prognostic
###     Subclasses in Adult Acute Myloid Leukemia
###
###     Bullinger et al. (2004), New England Journal of Medicine 350 (16),
###     1605-1616.
###
###############################################

if (!require("pamr"))
    stop("cannot attach package ", sQuote("pamr"))

### file as downloaded from http://www.ncbi.nlm.nih.gov/geo, 
### accession number GSE425
### use Format: SOFT to download this file
txt <- readLines(url("http://www.imbe.med.uni-erlangen.de/~hothorn/GSE425.txt"))

### separate clinical and expression data and selected cDNA's
writeLines(txt[121:237], con = "clinical.txt")
writeLines(txt[265:6548], con = "expressions.txt")
writeLines(txt[6788:6937], con = "selections.txt")

### clinical data
clinical <- read.table("clinical.txt", header = TRUE, sep = "\t",
                       strip.white = TRUE)
file.remove("clinical.txt")
rownames(clinical) <- gsub(" ", ".", as.character(clinical$Sample))
clinical$Sample <- NULL

### Study: Treatment Protocol 
### (HD98A between 16 and 60 years, HD98B > 60 years)
summary(clinical$Study)

### Sample.source. (BM: Bone marrow, PB: peripheral blood)
summary(clinical$Sample.source.)

### Cytogenetic.group (Chromosomal abnormalities)
summary(clinical$Cytogenetic.group)

### Karyotype ?
summary(clinical$Karyotype)

### FISH ?
summary(clinical$FISH)

### FLT3.aberration. -> mutation
summary(clinical$FLT3.aberration.)

### MLL.PTD  -> mutation
summary(clinical$MLL.PTD) 

### Sex
summary(clinical$Sex)

### Age (in years)
names(clinical) <- gsub("Age..years..", "Age", names(clinical))
summary(clinical$Age)

### Preceding.Malignancy.
summary(clinical$Preceding.Malignancy.)

### WBC: white blood cell count
names(clinical) <- gsub("WBC..x1000.ul..", "WBC", names(clinical))
summary(clinical$WBC)

### PB.Blasts (blood?)
names(clinical) <- gsub("PB.Blasts......", "PB.Blasts", names(clinical))
summary(clinical$PB.Blasts)

### BM.Blasts (bone marrow?)
names(clinical) <- gsub("BM.Blasts.....", "BM.Blasts", names(clinical))
summary(clinical$BM.Blasts)

### LDH
names(clinical) <- gsub("LDH..U.l..", "LDH", names(clinical))
summary(clinical$LDH)

### FAB.subtype. (French-American-British classification)
summary(clinical$FAB.subtype.)

### Tx.group: Treatment
### HSCT hematopoietic stem-cell transplantation, 
### IC intensive chemotherapy, AUTO autologous HSCT, and ALLO allogeneic HSCT.
summary(clinical$Tx.Group.)

### Remission.status.
### CR complete response, PR partial response, RD refractory disease, R randomization
summary(clinical$Remission.status.)

### Censoring Status
clinical$event <- clinical$Status == "dead"
clinical$Status <- NULL
summary(clinical$event)

### Survival Time (in days)
clinical$time <- clinical[["Overall.Survival..days.."]]
clinical[["Overall.Survival..days.."]] <- NULL
### hm, time = 0 induces problems
clinical$time[clinical$time == 0] <- 0.1
summary(clinical$time)

### = 333
median(clinical$time)
### = 611
median(clinical$time[!clinical$event])

### Training.Test.Set. (for Bullinger paper only)
summary(clinical$Training.Test.Set.)

clinical_predictors <- names(clinical)[-c(1, 2, 17:20)]
dim(clinical)

### expression data
expressions <- read.table("expressions.txt", header = TRUE, sep = "\t", 
                           quote = "", comment.char = "", 
                           colClasses = c(rep("character", 2), rep("numeric", 119)))
file.remove("expressions.txt")
CLID <- expressions$CLID

expressions <- t(expressions[,-(1:2)])
colnames(expressions) <- CLID

expressions <- expressions[match(rownames(clinical), rownames(expressions)),]

dim(expressions)

### impute missing values, use package `pamr'
iexpressions <- t(pamr.knnimpute(list(x = t(expressions)))$x)

### Selected genes
selgenes <- read.table("selections.txt", header = TRUE, sep = "\t",
                       quote = "", comment.char = "", na.string = "", 
                       colClasses = rep("character", 3))[,-1]
file.remove("selections.txt")
selgenes <- selgenes[!is.na(selgenes$Clone.ID),]
dim(selgenes)
summary(selgenes)
sum(!is.na(selgenes$Gene.Symbol)) ### = 133

### no expression levels?
selgenes[!selgenes$Clone.ID %in% colnames(expressions),]

### save data
save(clinical, clinical_predictors, expressions, iexpressions, selgenes,
     file = "AML_Bullinger.Rda")
