library(ape)
library(ade4)
library(seqinr)

setwd("~/Documents/SV.Orchids")

fig.1=read.csv("Fig.1.csv", header=T)
fig.1.a=fig.1[1:2,1:4]
rownames(fig.1.a)=fig.1.a[,1]
fig.1.a=fig.1.a[,-1]
fig.1.a=as.matrix(fig.1.a)
fig.1.b=fig.1[1:2,6:9]
rownames(fig.1.b)=fig.1.b[,1]
fig.1.b=fig.1.b[,-1]
fig.1.b=as.matrix(fig.1.b)

pdf("rplot.pdf") 
fig.1.a.plot=barplot(fig.1.a, xlab = "Match Types for Species on GenBank", legend=rownames(fig.1.a), ylab= "Number of Matches",
                     names=c("Species Match", "Genus Match", "No Match"),beside = T, ylim=c(0,55), cex.lab=1.5, cex.axis = 1.5, cex.names=1.5)
dev.off()
pdf("rplot.pdf") 
fig.1.b.plot=barplot(fig.1.b, xlab = "Match Types for Genera on GenBank", legend=rownames(fig.1.b), ylab= "Number of Matches",
                     names=c("Genus Match", "Multiple Matches", "No Match"),beside = T, ylim=c(0,55), cex.axis = 1.5, cex.names = 1.5, cex.lab=1.5)
dev.off()

fig.2=read.csv("Fig.2.csv", header=T, row.names = 1)
fig.2=as.matrix(fig.2)
pdf("rplot.pdf") 
fig.2.plot=barplot(fig.2, xlab = "All-to-All Match Types", legend=rownames(fig.2), ylab= "Number of Matches",cex.lab=1.5,
                     names=c("Species Match", "Genus Match", "Unresolved"),beside = T, ylim=c(0,120), cex.axis = 1.5, cex.names = 1.5)
dev.off()


# Building the phylogenies

# export .fasta files from Geneious

# run MAFFT alignment in Terminal
# place .fasta files into home/username folder
# Terminal code = mafft rbcl.only.fasta > rbcl.aln.fasta
# Terminal code = mafft matk.only.fasta > matk.aln.fasta
# Outgroup code = mafft rbcl.outgroup.fasta > rbcl.outgp.aln.fasta
# Outgroup code = mafft matk.outgroup.fasta > matk.outgp.aln.fasta

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

rbcl.og.aln=read.dna("rbcl.outgp.aln.fasta", format="fasta")
matk.og.aln=read.dna("matk.outgp.aln.fasta", format="fasta")

# Distance between sequences in alignment
rbcl.aln.2=read.alignment("rbcl.aln.fasta", format="fasta")

rbcl.aln.dist=as.matrix(dist.alignment(rbcl.aln.2, matrix = "identity"))

write.csv(rbcl.aln.dist, file="rbcl.aln.dist.mat.csv")

matk.aln.2=read.alignment("matk.aln.fasta", format="fasta")

matk.aln.dist=as.matrix(dist.alignment(matk.aln.2, matrix = "identity"))

write.csv(matk.aln.dist, file="matk.aln.dist.mat.csv")


# neighbor joining tree

nj.tree.rbcl=njs(dist.dna(rbcl.aln, model = "K80"))
plot(nj.tree.rbcl, cex=0.5)

nj.tree.rbcl.og=njs(dist.dna(rbcl.og.aln, model = "K80"))
plot(root(nj.tree.rbcl.og,133),cex = 0.3)

nj.tree.matk.og=njs(dist.dna(matk.og.aln, model = "K80"))
plot(root(nj.tree.matk.og,132),cex = 0.3)

# Modeltest

rbcl.phydat=read.phyDat("rbcl.aln.fasta", format = "fasta", type = "DNA")
matk.phydat=read.phyDat("matk.aln.fasta", format = "fasta", type = "DNA")

rbcl.phydat.og=read.phyDat("rbcl.outgp.aln.fasta", format = "fasta", type = "DNA")
matk.phydat.og=read.phyDat("matk.outgp.aln.fasta", format = "fasta", type = "DNA")

modelTest(rbcl.phydat) # GTR + G + I is the best
modelTest(matk.phydat)# GTR + G + I is the best

modelTest(rbcl.phydat.og) # GTR + G + I is the best
modelTest(matk.phydat.og)# GTR + G + I is the best

# Maximum likelihood tree

dist.rbcl=dist.logDet(rbcl.phydat)
dist.matk=dist.logDet(matk.phydat)

dist.rbcl.og=dist.logDet(rbcl.phydat.og)
dist.matk.og=dist.logDet(matk.phydat.og)

rbcl.nj.tree=NJ(dist.rbcl)
matk.nj.tree=NJ(dist.matk)

rbcl.nj.tree.og=NJ(dist.rbcl.og)
matk.nj.tree.og=NJ(dist.matk.og)

rbcl.ml.model=pml(rbcl.nj.tree, rbcl.phydat)
matk.ml.model=pml(matk.nj.tree, matk.phydat)

rbcl.ml.model.og=pml(rbcl.nj.tree.og, rbcl.phydat.og)
matk.ml.model.og=pml(matk.nj.tree.og, matk.phydat.og)

rbcl.ml.tree=optim.pml(rbcl.ml.model, model="GTR", optGamma = TRUE, optInv = TRUE)
matk.ml.tree=optim.pml(matk.ml.model, model="GTR", optGamma = TRUE, optInv = TRUE)

rbcl.ml.tree.og=optim.pml(rbcl.ml.model.og, model="GTR", optGamma = TRUE, optInv = TRUE)
matk.ml.tree.og=optim.pml(matk.ml.model.og, model="GTR", optGamma = TRUE, optInv = TRUE)

rbcl.optim.ml.tree=optim.pml(rbcl.ml.model, model="GTR", optGamma = TRUE, optInv = TRUE, optNni = TRUE)
matk.optim.ml.tree=optim.pml(matk.ml.model, model="GTR", optGamma = TRUE, optInv = TRUE, optNni = TRUE)

rbcl.optim.ml.tree.og=optim.pml(rbcl.ml.model.og, model="GTR", optGamma = TRUE, optInv = TRUE, optNni = TRUE)
matk.optim.ml.tree.og=optim.pml(matk.ml.model.og, model="GTR", optGamma = TRUE, optInv = TRUE, optNni = TRUE)

rbcl.optim.ml.tree.og.rooted=root(rbcl.optim.ml.tree.og$tree, outgroup="FBPL464-12|Arabidopsis thaliana|rbcLa", resolve.root=T)
matk.optim.ml.tree.og.rooted=root(matk.optim.ml.tree.og$tree, outgroup="FBPL464-12|Arabidopsis_2 thaliana|matK", resolve.root=T)

plot(midpoint(rbcl.optim.ml.tree.og.rooted), cex=0.3)
plot(midpoint(matk.optim.ml.tree.og.rooted), cex=0.3)

# Bootstrap trees
# Remove NA node labels in text editor

rbcl.boot.og=bootstrap.pml(rbcl.optim.ml.tree.og, bs=1000, optNni=TRUE)
rbclbootsrooted <- lapply(rbcl.boot.og, function(x) root(x, resolve.root=TRUE, outgroup="FBPL464-12|Arabidopsis thaliana|rbcLa"))
class(rbclbootsrooted) <- "multiPhylo"
tree.rbcl.og=plotBS(rbcl.optim.ml.tree.og$tree, rbclbootsrooted, p=75, type="phylo")
tree.rbcl.og$tip.label=rbcl.tips
write.tree(tree.rbcl.og, file="tree.rbcl.bs.og.txt")

matk.boot.og=bootstrap.pml(matk.optim.ml.tree.og, bs=1000, optNni=TRUE)
matkbootsrooted <- lapply(matk.boot.og, function(x) root(x, resolve.root=TRUE, outgroup="FBPL464-12|Arabidopsis_2 thaliana|matK"))
class(matkbootsrooted) <- "multiPhylo"
tree.matk.og=plotBS(matk.optim.ml.tree.og$tree, matkbootsrooted, p=75, type="phylo")
tree.matk.og$tip.label=matk.tips
write.tree(tree.matk.og, file="tree.matk.bs.og.txt")


# Genus tree

# Building the phylogenies

# export .fasta files from Geneious

# run MAFFT alignment in Terminal
# place .fasta files into home/username folder
# Terminal code = mafft rbcl.genus.fasta > rbcl.genus.aln.fasta
# Terminal code = mafft matk.genus.fasta > matk.genus.aln.fasta

setwd("~/Documents/SV.Orchids/PhyloFiles")

rbcl.genus.aln=read.dna("rbcl.genus.aln.fasta", format="fasta")
matk.genus.aln=read.dna("matk.genus.aln.fasta", format="fasta")
concat=cbind(rbcl.genus.aln, matk.genus.aln, fill.with.gaps=TRUE, check.names=FALSE)

write.dna(concat, "concat.fasta", format="fasta")

# Add two sequences that only had matk to the concatenated alignmentment in Terminal
# mafft --auto --addfragments new.seq.fasta concat.fasta>new.aln.fasta

all.genus.aln=read.dna("new.aln.fasta", format="fasta")
all.genus.phydat=read.phyDat("new.aln.fasta", format = "fasta", type="DNA")

# neighbor joining tree

nj.tree.rbcl=njs(dist.dna(rbcl.aln, model = "K80"))
plot(nj.tree.rbcl, cex=0.5)

nj.tree.rbcl.og=njs(dist.dna(rbcl.og.aln, model = "K80"))
plot(root(nj.tree.rbcl.og,133),cex = 0.3)

nj.tree.matk.og=njs(dist.dna(matk.og.aln, model = "K80"))
plot(root(nj.tree.matk.og,132),cex = 0.3)

# Modeltest
modelTest(all.genus.phydat) # GTR + G + I is the best

# Maximum likelihood tree

dist.genus=dist.logDet(all.genus.phydat)
genus.nj.tree=NJ(dist.genus)
genus.ml.model=pml(genus.nj.tree, all.genus.phydat)
genus.ml.tree=optim.pml(genus.ml.model, model="GTR", optGamma = TRUE, optInv = TRUE)
genus.optim.ml.tree=optim.pml(genus.ml.model, model="GTR", optGamma = TRUE, optInv = TRUE, optNni = TRUE)
genus.optim.ml.tree.rooted=root(genus.optim.ml.tree$tree, outgroup="FBPL464-12|Arabidopsis thaliana|rbcLa", resolve.root=T)
plot(midpoint(genus.optim.ml.tree.rooted), cex=0.3)

# Bootstrap trees
# Remove NA node labels in text editor

genus.boot=bootstrap.pml(genus.optim.ml.tree, bs=1000, optNni=TRUE)
genusbootsrooted <- lapply(genus.boot, function(x) root(x, resolve.root=TRUE, outgroup="FBPL464-12|Arabidopsis thaliana|rbcLa"))
class(genusbootsrooted) <- "multiPhylo"
tree.genus=plotBS(genus.optim.ml.tree$tree, genusbootsrooted, p=75, type="phylo")
tree.genus$tip.label=genus.tips
write.tree(tree.genus, file="tree.genus.bs.txt")


