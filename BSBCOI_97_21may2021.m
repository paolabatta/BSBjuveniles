    set.dir(input=/home/FCAM/mwojcicki/BSB/COIFGBSB)
set.current(processors=36)

make.file(inputdir=/home/FCAM/mwojcicki/BSB/COIFGBSB, type=fastq, prefix=COIFGBSB_97)
make.contigs(file=current, insert=30)
summary.seqs(fasta=current)
screen.seqs(fasta=current, group=current, summary=current, minlength=155, maxambig=0)
summary.seqs(fasta=current)
unique.seqs(fasta=current)
summary.seqs(fasta=current, name=current)
count.seqs(name=current, group=current)
summary.seqs(fasta=current, count=current)

##using MetaZooGene Global database as reference
align.seqs(fasta=current, reference=/labs/Bucklin/Databases/MZGCombinedGlobal7March2021MAFFT.fasta, flip=t)
summary.seqs(fasta=current, count=current)
screen.seqs(fasta=current, count=current, summary=current, start=36086, end=23506, minlength=155, maxhomop=16)
summary.seqs(fasta=current, count=current)
filter.seqs(fasta=current, vertical=T)
unique.seqs(fasta=current, count=current)
summary.seqs(fasta=current, count=current)
get.current()

##diffs=2 is 0.005 for COI with 344 bp average length
pre.cluster(fasta=current, count=current, diffs=2, method=unoise)
chimera.vsearch(fasta=current, count=current, dereplicate=t)
remove.seqs(fasta=current, accnos=current, dups=f)
summary.seqs(fasta=current, count=current)
classify.seqs(fasta=current, count=current, reference=/labs/Bucklin/Databases/MZGdbCombined_7March2021_NAtl_A.fasta,taxonomy=/labs/Bucklin/Databases/MZGdbCombined_7March2021_NAtl_A.mothur, cutoff=97)

##need version 1.44 or newer to do make.shared for asv
##For ASVs, convert fasta and count table from pre.cluster to generate shared and list asv files
##adding80 to label=asv to indicate the cutoff used
make.shared(count=current, label=asv97)
classify.otu(list=current, count=current, taxonomy=current, label=asv97)
summary.tax(taxonomy=current, count=current)
summary.seqs(fasta=current, count=current)

##PLOT RAREFACTION OUTSIDE OF MOTHUR USING iNEXT https://chao.shinyapps.io/iNEXTOnline/
rarefaction.single(shared=current, calc=sobs)

##alpha diversity measures
summary.single(shared=current, calc=nseqs-sobs-coverage-shannon-shannoneven-invsimpson, label=asv97)

##beta diversity measures
dist.shared(shared=current, calc=braycurtis-thetayc-jest, subsample=T, label=asv97)
split.groups(fasta=current, count=current)
