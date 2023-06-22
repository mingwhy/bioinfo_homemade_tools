#http://www.bioconductor.org/packages/devel/bioc/vignettes/BASiCS/inst/doc/BASiCS.html

#BiocManager::install("BASiCS")
library(BASiCS)

# quick start 
set.seed(1)
Data <- makeExampleBASiCS_Data()
Data #50 gene x 30 cell
colnames(colData(Data)) #BatchInfo is included
altExp(Data) #spike-in, 20gene x 30cell

Chain <- BASiCS_MCMC(
  Data = Data,
  N = 20000, Thin = 20, Burn = 10000,
  PrintProgress = FALSE, Regression = TRUE
)
#Please ensure the acceptance rates displayed in the console output of BASiCS_MCMC are around 0.44. 
#If they are too far from this value, you should increase N and Burn.
Chain


# complete workflow
set.seed(1)
Counts <- matrix(rpois(50 * 40, 2), ncol = 40)
rownames(Counts) <- c(paste0("Gene", 1:40), paste0("Spike", 1:10))
colnames(Counts) <- paste0("Cell", 1:40)

Tech <- c(rep(FALSE, 40), rep(TRUE, 10)) #If Tech[i] = FALSE the gene i is biological; otherwise the gene is spike-in. This vector must be specified in the same order of genes as in the Counts matrix.
set.seed(2)
SpikeInput <- rgamma(10, 1, 1)
SpikeInfo <- data.frame(
  "SpikeID" = paste0("Spike", 1:10),
  "SpikeInput" = SpikeInput
)

# No batch structure
DataExample <- newBASiCS_Data(Counts, Tech, SpikeInfo)

# With batch structure
DataExample <- newBASiCS_Data(
  Counts, Tech, SpikeInfo,
  BatchInfo = rep(c(1, 2), each = 20)
)

# without spike-in genes
set.seed(1)
CountsNoSpikes <- matrix(rpois(50 * 40, 2), ncol = 40)
rownames(CountsNoSpikes) <- paste0("Gene", 1:50)
colnames(CountsNoSpikes) <- paste0("Cell", 1:40)

# With batch structure
DataExampleNoSpikes <- SingleCellExperiment(
  assays = list(counts = CountsNoSpikes),
  colData = data.frame(BatchInfo = rep(c(1, 2), each = 20))
)

# Analysis for two groups of cells
data(ChainSC)
data(ChainRNA)
ChainSC;
ChainRNA
Test <- BASiCS_TestDE(
  Chain1 = ChainSC, Chain2 = ChainRNA,
  GroupLabel1 = "SC", GroupLabel2 = "PaS",
  EpsilonM = log2(1.5), EpsilonD = log2(1.5),
  EFDR_M = 0.10, EFDR_D = 0.10,
  Offset = TRUE, PlotOffset = TRUE, Plot = TRUE
)
Test
head(as.data.frame(Test, Parameter = "Mean"))
head(as.data.frame(Test, Parameter = "Disp"))
rowData(Test)

rowData(Test) <- cbind(rowData(Test), Index = 1:nrow(rowData(Test)))
as.data.frame(Test, Parameter = "Mean")

BASiCS_PlotDE(Test)
BASiCS_PlotDE(Test, Plots = c("MA", "Volcano"))
BASiCS_PlotDE(Test, Plots = "MA", Parameters = "Mean")

Test <- BASiCS_TestDE(
  Chain1 = ChainSC, Chain2 = ChainRNA,
  GroupLabel1 = "SC", GroupLabel2 = "PaS",
  EpsilonM = 0, EpsilonD = log2(1.5),
  EFDR_M = 0.10, EFDR_D = 0.10,
  Offset = TRUE, PlotOffset = FALSE, Plot = FALSE
)

# If WithSpikes = FALSE
DataNoSpikes <- SingleCellExperiment(
  assays = list(counts = Counts),
  colData = data.frame(BatchInfo = rep(c(1, 2), each = 20))
)

ChainNoSpikes <- BASiCS_MCMC(
  Data = DataNoSpikes, N = 1000,
  Thin = 10, Burn = 500,
  WithSpikes = FALSE,  Regression = TRUE,
  PrintProgress = FALSE
)

# If Regression = TRUE
DataRegression <- newBASiCS_Data(
  Counts, Tech, SpikeInfo,
  BatchInfo = rep(c(1, 2), each = 20)
)
ChainRegression <- BASiCS_MCMC(
  Data = DataRegression, N = 1000,
  Thin = 10, Burn = 500,
  Regression = TRUE,
  PrintProgress = FALSE
)
data("ChainRNAReg")
BASiCS_ShowFit(ChainRNAReg)

