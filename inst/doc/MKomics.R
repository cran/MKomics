## ---- eval = FALSE------------------------------------------------------------
#  ## Install package BiocManager
#  install.packages("BiocManager")
#  ## Use BiocManager to install limma
#  BiocManager::install("limma")

## -----------------------------------------------------------------------------
library(MKomics)

## -----------------------------------------------------------------------------
## One-sample test
X <- matrix(rnorm(10*20, mean = 1), nrow = 10, ncol = 20)
mod.t.test(X)

## Two-sample test
set.seed(123)
X <- rbind(matrix(rnorm(5*20), nrow = 5, ncol = 20),
           matrix(rnorm(5*20, mean = 1), nrow = 5, ncol = 20))
g2 <- factor(c(rep("group 1", 10), rep("group 2", 10)))
mod.t.test(X, group = g2)

## Paired two-sample test
subjID <- factor(rep(1:10, 2))
mod.t.test(X, group = g2, paired = TRUE, subject = subjID)

## -----------------------------------------------------------------------------
set.seed(123)
X <- rbind(matrix(rnorm(5*20), nrow = 5, ncol = 20),
           matrix(rnorm(5*20, mean = 1), nrow = 5, ncol = 20))
gr <- factor(c(rep("A1", 5), rep("B2", 5), rep("C3", 5), rep("D4", 5)))
mod.oneway.test(X, gr)

## -----------------------------------------------------------------------------
X <- rbind(matrix(rnorm(6*18), nrow = 6, ncol = 18),
           matrix(rnorm(6*18, mean = 1), nrow = 6, ncol = 18))
gr <- factor(c(rep("T1", 6), rep("T2", 6), rep("T3", 6)))
subjectID <- factor(c(rep(1:6, 3)))
mod.oneway.test(X, gr, repeated = TRUE, subject = subjectID)

## -----------------------------------------------------------------------------
pairwise.mod.t.test(X, gr)

## -----------------------------------------------------------------------------
x <- rnorm(100) ## assumed as log-data
g <- factor(sample(1:4, 100, replace = TRUE))
levels(g) <- c("a", "b", "c", "d")
## modified FC
pairwise.fc(x, g)
## "true" FC
pairwise.fc(x, g, mod.fc = FALSE)
## without any transformation
pairwise.logfc(x, g)

## -----------------------------------------------------------------------------
pairwise.logfc(x, g, ave = median)

## -----------------------------------------------------------------------------
M <- matrix(rnorm(100), ncol = 5)
FL <- matrix(rpois(100, lambda = 10), ncol = 5) # only for this example
repMeans(x = M, flags = FL, use.flags = "max", ndups = 5, spacing = 4)

## -----------------------------------------------------------------------------
af <- oneWayAnova(c(rep(1,5),rep(2,5)))
## p value
af(rnorm(10))
x <- matrix(rnorm(12*10), nrow = 10)
## 2-way ANOVA with interaction
af1 <- twoWayAnova(c(rep(1,6),rep(2,6)), rep(c(rep(1,3), rep(2,3)), 2))
## p values
apply(x, 1, af1)
## 2-way ANOVA without interaction
af2 <- twoWayAnova(c(rep(1,6),rep(2,6)), rep(c(rep(1,3), rep(2,3)), 2), 
                   interaction = FALSE)
## p values
apply(x, 1, af2)

## -----------------------------------------------------------------------------
M <- matrix(rcauchy(1000), nrow = 5)
## Pearson
corDist(M)
## Spearman
corDist(M, method = "spearman")
## Minimum Covariance Determinant
corDist(M, method = "mcd")

## -----------------------------------------------------------------------------
madMatrix(t(M))

## ---- fig.width=8, fig.height=7-----------------------------------------------
M <- matrix(rnorm(1000), ncol = 20)
colnames(M) <- paste("Sample", 1:20)
M.cor <- cor(M)
corPlot(M.cor, minCor = min(M.cor), labels = colnames(M))

## ---- fig.width=8, fig.height=7-----------------------------------------------
corPlot2(M.cor, minCor = min(M.cor), labels = colnames(M))

## ---- fig.width=8, fig.height=7-----------------------------------------------
## random data
x <- matrix(rnorm(1000), ncol = 10)
## outliers
x[1:20,5] <- x[1:20,5] + 10
madPlot(x, new = TRUE, maxMAD = 2.5, labels = TRUE,
        title = "MAD: Outlier visible")
madPlot2(x, new = TRUE, maxMAD = 2.5, labels = TRUE,
        title = "MAD: Outlier visible")
## in contrast
corPlot2(x, new = TRUE, minCor = -0.5, labels = TRUE,
        title = "Correlation: Outlier masked")

## ---- fig.width=7, fig.height=9-----------------------------------------------
## generate some random data
data.plot <- matrix(rnorm(100*50, sd = 1), ncol = 50)
colnames(data.plot) <- paste("patient", 1:50)
rownames(data.plot) <- paste("gene", 1:100)
data.plot[1:70, 1:30] <- data.plot[1:70, 1:30] + 3
data.plot[71:100, 31:50] <- data.plot[71:100, 31:50] - 1.4
data.plot[1:70, 31:50] <- rnorm(1400, sd = 1.2)
data.plot[71:100, 1:30] <- rnorm(900, sd = 1.2)
nrcol <- 128
## Load required packages
library(RColorBrewer)
myCol <- rev(colorRampPalette(brewer.pal(10, "RdBu"))(nrcol))
heatmap(data.plot, col =  myCol, main = "standard colors")
myCol2 <- heatmapCol(data = data.plot, col = myCol, 
                     lim = min(abs(range(data.plot)))-1)
heatmap(data.plot, col = myCol2, main = "heatmapCol colors")

## -----------------------------------------------------------------------------
x <- "GACGGATTATG"
y <- "GATCGGAATAG"
## Hamming distance
stringDist(x, y)
## Levenshtein distance
d <- stringDist(x, y)
d

## -----------------------------------------------------------------------------
attr(d, "ScoringMatrix")
attr(d, "TraceBackMatrix")

## -----------------------------------------------------------------------------
## optimal global alignment score
d <- stringSim(x, y)
d
attr(d, "ScoringMatrix")
attr(d, "TraceBackMatrix")

## optimal local alignment score
d <- stringSim(x, y, global = FALSE)
d
attr(d, "ScoringMatrix")
attr(d, "TraceBackMatrix")

## -----------------------------------------------------------------------------
x <- "GACGGATTATG"
y <- "GATCGGAATAG"
## Levenshtein distance
d <- stringDist(x, y)
## optimal global alignment
traceBack(d)

## Optimal global alignment score
d <- stringSim(x, y)
## optimal global alignment
traceBack(d)

## Optimal local alignment score
d <- stringSim(x, y, global = FALSE)
## optimal local alignment
traceBack(d, global = FALSE)

## -----------------------------------------------------------------------------
sessionInfo()

