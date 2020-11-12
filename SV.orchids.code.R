library(ape)
library(ade4)
library(seqinr)

setwd("~/Documents/SV.Orchids")

fig.1=read.csv("Fig.1.csv", header=T)
fig.1.a=fig.1[1:2,1:4]
rownames(fig.1.a)=fig.1.a[,1]
fig.1.a=fig.1.a[,-1]
fig.1.a=as.matrix(fig.1.a)
fig.1.b=fig.1[1:2,6:8]
rownames(fig.1.b)=fig.1.b[,1]
fig.1.b=fig.1.b[,-1]
fig.1.b=as.matrix(fig.1.b)

pdf("rplot.pdf") 
fig.1.a.plot=barplot(fig.1.a, xlab = "Match Types for Species on GenBank", legend=rownames(fig.1.a), ylab= "Number of Matches",
                     names=c("Species Match", "Genus Match", "No Match"),beside = T, ylim=c(0,80))
dev.off()
pdf("rplot.pdf") 
fig.1.b.plot=barplot(fig.1.b, xlab = "Match Types for Genera on GenBank", legend=rownames(fig.1.b), ylab= "Number of Matches",
                     names=c("Genus Match", "No Match"),beside = T, ylim=c(0,80))
dev.off()

fig.2=read.csv("Fig.2.csv", header=T, row.names = 1)
fig.2=as.matrix(fig.2)
pdf("rplot.pdf") 
fig.2.plot=barplot(fig.2, xlab = "All-to-All Match Types", legend=rownames(fig.2), ylab= "Number of Matches",
                     names=c("Species Match", "Genus Match", "Combo Match"),beside = T, ylim=c(0,120))
dev.off()


# Building the phylogenies

# export .fasta files from Geneious

# run MAFFT alignment in Terminal
# place .fasta files into home/username folder
# Terminal code = mafft rbcl.only.fasta > rbcl.aln.fasta
# Terminal code = mafft matk.only.fasta > matk.aln.fasta

library(ape)
library(picante)
library(geiger)
library(vegan)
library(seqinr)
library(phangorn)
library(phytools)
# For R crashing
library(uwot)
library(RcppParallel)

setwd("~/Documents/SV.Orchids/PhyloFiles")

rbcl.aln=read.dna("rbcl.aln.fasta", format="fasta")
matk.aln=read.dna("matk.aln.fasta", format="fasta")

# neighbor joining tree

nj.tree.rbcl=njs(dist.dna(rbcl.aln, model = "K80"))
plot(nj.tree.rbcl, cex=0.5)

# Modeltest

rbcl.phydat=read.phyDat("rbcl.aln.fasta", format = "fasta", type = "DNA")
matk.phydat=read.phyDat("matk.aln.fasta", format = "fasta", type = "DNA")

modelTest(rbcl.phydat) # GTR + G + I is the best
modelTest(matk.phydat)# GTR + G + I is the best

# Maximum likelihood tree

dist.rbcl=dist.logDet(rbcl.phydat)
dist.matk=dist.logDet(matk.phydat)

rbcl.nj.tree=NJ(dist.rbcl)
matk.nj.tree=NJ(dist.matk)

rbcl.ml.model=pml(rbcl.nj.tree, rbcl.phydat)
matk.ml.model=pml(matk.nj.tree, matk.phydat)

rbcl.ml.tree=optim.pml(rbcl.ml.model, model="GTR", optGamma = TRUE, optInv = TRUE)
matk.ml.tree=optim.pml(matk.ml.model, model="GTR", optGamma = TRUE, optInv = TRUE)

rbcl.optim.ml.tree=optim.pml(rbcl.ml.model, model="GTR", optGamma = TRUE, optInv = TRUE, optNni = TRUE)
matk.optim.ml.tree=optim.pml(matk.ml.model, model="GTR", optGamma = TRUE, optInv = TRUE, optNni = TRUE)

plot(midpoint(rbcl.optim.ml.tree$tree), cex=0.3)
plot(midpoint(matk.optim.ml.tree$tree), cex=0.3)

plot.phylo(rbcl.optim.ml.tree$tree, cex=0.3)

# Bootstrap trees
# Remove NA node labels in text editor

rbcl.boot=bootstrap.pml(rbcl.optim.ml.tree, bs=1000, optNni=TRUE)
tree.rbcl=plotBS(rbcl.optim.ml.tree$tree, rbcl.boot, p=75, type="cladogram")
tree.rbcl$tip.label=rbcl.tips
write.tree(tree.rbcl, file="tree.rbcl.bs.txt")

matk.boot=bootstrap.pml(matk.optim.ml.tree, bs=1000, optNni=TRUE)
tree.matk=plotBS(matk.optim.ml.tree$tree, matk.boot, p=75, type="cladogram")
tree.matk$tip.label=matk.tips
write.tree(tree.matk, file="tree.matk.bs.txt")
