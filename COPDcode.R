
# Reads in data
exp.copd <- read.table("ECLIPSE_Blood_Exp_COPD.txt")
exp.smok <- read.table("ECLIPSE_Blood_Exp_SmokerControl.txt")

# Stick the tables together side by side into one megatable
full.exp <- cbind(exp.smok,exp.copd)

# Load limma
library(limma)
library(tidyverse)


# Makes a big factor list that looks like "CTRL CTRL CTRL.... COPD COPD COPD..."
lev <- c("CTRL","COPD")
f <- factor(c(rep("CTRL",ncol(exp.smok)),rep("COPD",ncol(exp.copd))),levels=lev)

# Creates a matrix that denotes which patient is COPD and which is CTRL
# I'm still foggy on the syntax of formulae and model.matrix
# I *think* that ~0+f means "Match 0 to all the factors, with all duplicates removed"
# That is, 0 gets tacked on to both COPD and CTRL. How the ordering works out, IDK 
# Finally, I know it's a one-sided formula operator. I have no idea what that actually means.
design <- model.matrix(~0+f)
colnames(design) <- lev

# With the knowledge of who's who, do a linear fit on the data
fit <- lmFit(full.exp,design)

cont.diff <- makeContrasts(COPD-CTRL, levels=design)
fit2 <- contrasts.fit(fit,cont.diff)
fit3 <- eBayes(fit2)

top.copd <- topTable(fit3,adjust="BH",number=100)

copd.deg <- rownames(top.copd)[top.copd[,5]<0.05]

copd.drivers <- t(read.table("COPD_drivers.txt"))[1,]

copd.drivers.onarray <- intersect(copd.drivers,rownames(exp.copd))

a <- length(intersect(copd.deg,copd.drivers.onarray))
b <- length(copd.deg)-a
c <- length(copd.drivers.onarray)-a
d <- nrow(exp.copd)-a-b-c

ftest <- fisher.test(array(c(a,b,c,d),dim=c(2,2)),alternative="greater") ## P = 0.71

#Find the means of all the rows:
copdMeans <- as.tibble(rowMeans(exp.copd))
copdMeans[,2] <- copdMeans[,1]
copdMeans[,1] <- rownames(exp.copd)
smokMeans <- as.tibble(rowMeans(exp.smok))
smokMeans[,2] <- smokMeans[,1]
smokMeans[,1] <- rownames(exp.smok)
topCopdRownames <- as.tibble(rownames(top.copd))
for (i in 1:nrow(topCopdRownames)){
  grepSearch[i] <- which(copdMeans[,1] == as.list(topCopdRownames[i,]))
}
topCopdMeans <- copdMeans[grepSearch,]
topSmokMeans <- smokMeans[grepSearch,] #Assumes the rows are in the same order

differenceBetween <- topCopdMeans[,2] - topSmokMeans[,2]
differenceBetween[,2] <- topCopdMeans[,1]
print(differenceBetween)

ggplot(differenceBetween, aes(x = value, y = value.1))+
  geom_tile(aes(fill = value.1), color = "white") +
  scale_fill_gradient(low = "white", high = "steelblue")
