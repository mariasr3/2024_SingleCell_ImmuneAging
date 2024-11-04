
# Data exploratory plots -------

library(ggplot2);library(dplyr)

#paths 
if (getwd()== "/Users/mariasopenar"){
  basepath <- "/Users/mariasopenar/cluster/"
  data_path <- "/Users/mariasopenar/cluster/Projects/scRNAseq/"
  
}else if(getwd()== "/home/mariasr"){
  basepath <- "/home/mariasr/cluster/"
  data_path <- "/home/mariasr/cluster/Projects/scRNAseq/"
}else{
  basepath <- "/gpfs/projects/bsc83/"
  data_path <- "/gpfs/projects/bsc83/Projects/scRNAseq/"
}
fig_path <- paste0(data_path, "/msopena/02_OneK1K_Age/figures/")

metadata_all <- readRDS(paste0(basepath, "Data/scRNAseq/Yazar2022/metadata_processed.rds"))
order_cells<- readRDS(paste0(basepath, "Data/scRNAseq/Yazar2022/order_cells.rds"))


metadata <- metadata_all[!duplicated(metadata_all$assignment),]

#Age distributioncat 
age <- ggplot(metadata, aes(x=Age))+geom_histogram(fill= blue)+theme+ylab("Number of donors")+geom_vline(xintercept = 63, color=blue, linetype="dashed")+theme(aspect.ratio=1, axis.text=element_text(size=10))

#Sex distribution
sex <- ggplot(metadata,  aes(x=Gender, fill=Gender))+geom_bar(stat="count")+theme+ylab("Number of donors")+xlab("Sex")+scale_fill_manual(values=c(red, blue), name="Sex")+theme(aspect.ratio=1, legend.position="none", 
                                                                                                                                                                                     axis.text=element_text(size=10))

#Number of cells l1 
ncells <- metadata_all %>% dplyr::group_by(predicted.celltype.l1) %>% dplyr::count()
colnames(ncells)[1] <- "celltype"

ncells_l1 <- ggplot(ncells, aes(y=celltype, x=n)) + geom_bar(stat="identity", fill=alpha(blue, 0.7))+  scale_x_continuous(labels = function(x) format(x / 1e6, scientific = FALSE)) +
  theme  + xlab("Number of cells l1 (million)") + ylab("") +scale_y_discrete(limits = rev(levels(ncells$celltype)))+
  theme(axis.text.x=element_text(size=10), axis.text.y=element_text(size=10), axis.title=element_text(size=11), aspect.ratio=1)

#Number of cells l2
metadata_all <- metadata_all[!is.na(metadata_all$cell_type),]
ncells <- metadata_all %>% dplyr::group_by(cell_type) %>% dplyr::count() %>% dplyr::arrange(desc(n))
colnames(ncells)[1] <- "celltype"
ncells <-ncells[ncells$n > 1000,]

#saveRDS(ncells$celltype, data_path, "/msopena/02_OneK1K_Age/robjects/cells_to_keep.rds"))


ncells_l2 <-ggplot(ncells, aes(y=reorder(celltype,-n), x=n)) + geom_bar(stat="identity", fill=alpha(blue, 0.7))+  scale_x_continuous(labels = function(x) format(x / 1e6, scientific = FALSE)) +
  theme  + xlab("Number of cells (million)") + ylab("") +scale_y_discrete(limits = levels(ncells$celltype))+coord_flip()+
  theme(axis.text.y=element_text(size=10), axis.title=element_text(size=10), aspect.ratio=0.6, axis.text.x=element_text(angle = 90, hjust = 0.95, vjust = 0.6, size=10))


#Number of donors per cell type 

ncells_donor <- metadata_all[metadata_all$cell_type %in% ncells$celltype,] %>%   distinct(cell_type, assignment) %>% dplyr::group_by(cell_type) %>% dplyr::count() %>% dplyr::arrange(desc(n))


ncells_donor_l2 <-ggplot(ncells_donor, aes(y=reorder(cell_type,-n), x=n)) + geom_bar(stat="identity", fill=alpha(blue, 0.7)) +
  theme  + xlab("Number of donors per cell type") + ylab("") +scale_y_discrete(limits = levels(ncells$cell_type))+coord_flip()+
  theme(axis.text.y=element_text(size=10), axis.title=element_text(size=10), aspect.ratio=0.6, axis.text.x=element_text(angle = 90, hjust = 0.95, vjust = 0.6, size=10))



# Number of cells per donor   
  
ncells_donor <- metadata_all %>% dplyr::group_by(assignment) %>% dplyr::count()
ncells_donor$cells_donor <- "NCells_donor"

ncells_donor_p <-ggplot(ncells_donor, aes(x=cells_donor, y=n)) + geom_boxplot(fill=alpha(blue, 0.7), outlier.shape = NA)+geom_jitter(size=0.5, color=blue)+
  theme  + xlab("Number of cells (million)") + ylab("") +
theme(axis.text.x=element_text(size=10), axis.text.y=element_text(size=12), axis.title=element_text(size=10), aspect.ratio=1)

median(ncells_donor$n)
ncells_donor_p <-ggplot(ncells_donor, aes(x=n )) + geom_density(fill=alpha(blue, 0.5), color=blue)+
  theme  + xlab("Number of cells per donor") + ylab("Density") +geom_vline(xintercept = median(ncells_donor$n), color=blue, linetype="dashed")+
  theme(axis.text.x=element_text(size=10), axis.text.y=element_text(size=10), axis.title=element_text(size=11), aspect.ratio=1)


#expressed genes 

# Figure S1A
pdf(paste0(fig_path, "/FigS1/S1A_AgeDistribution.pdf"), width = 3, height = 3)
age 
dev.off()

# Figure S1B
pdf(paste0(fig_path, "/FigS1/S1B_SexDistribution.pdf"), width = 3, height = 3)
sex 
dev.off()

# Figure S1C
pdf(paste0(fig_path, "/FigS1/S1C_NumCells_l1.pdf"), width = 3, height = 3)
ncells_l1 
dev.off()
# Figure S1D
pdf(paste0(fig_path, "/FigS1/S1D_NumCells_l2.pdf"), width = 6, height = 3)
ncells_l2
dev.off()

pdf(paste0(fig_path, "/FigS1/S1E_NumCells_Donor.pdf"), width = 3, height = 3)
ncells_donor_p 
dev.off()

pdf(paste0(fig_path, "/FigS1/S1E_NumDonor_Celltype.pdf"), width = 6, height = 3)
ncells_donor_l2 
dev.off()








