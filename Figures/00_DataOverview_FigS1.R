
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
all_genes <- readRDS(paste0(data_path, "msopena/02_OneK1K_Age/robjects/01_DEG_pseudobulk/DEG_Age_cell_type_02.rds"))
cells_to_keep <- unique(all_genes$celltype)

metadata <- metadata_all[!duplicated(metadata_all$assignment),]


metadata <- metadata_all[!duplicated(metadata_all$assignment),]

#Age distribution
age <- ggplot(metadata, aes(x=Age))+geom_histogram(fill= blue)+theme+ylab("Number of donors")+geom_vline(xintercept = 63, color=blue, linetype="dashed")+theme(aspect.ratio=1, axis.text=element_text(size=10))

## (add) stratified by sex  

metadata.mod <- metadata
metadata.mod$Gender <- factor(metadata.mod$Gender,
                              levels = c('F', 'M')) #F=blue, M=red

metadata_by_sex.stats <- metadata.mod %>%
  group_by(Gender) %>%
  summarize(median = median(Age))

age_by_sex <- ggplot(metadata.mod, aes(x=Age, fill=Gender))+
  geom_histogram(alpha = 0.7)+
  theme+
  ylab("Number of donors")+
  geom_vline(data = metadata_by_sex.stats, 
             aes(xintercept = median, color = Gender),
             linetype="dashed")+
  scale_color_manual(values=c(blue, red), name="Sex")+
  scale_fill_manual(values=c(blue, red), name="Sex")+
  theme(aspect.ratio=1, axis.text=element_text(size=10))

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
ncells <- metadata_all[metadata_all$cell_type %in% cells_to_keep,] %>% dplyr::group_by(cell_type) %>% dplyr::count() %>% dplyr::arrange(desc(n))
colnames(ncells)[1] <- "celltype"
l2_order <- unique(ncells$celltype)

#saveRDS(ncells$celltype, data_path, "/msopena/02_OneK1K_Age/robjects/cells_to_keep.rds"))


ncells_l2 <-ggplot(ncells, aes(y=reorder(celltype,-n), x=n)) + geom_bar(stat="identity", fill=alpha(blue, 0.7))+  scale_x_continuous(labels = function(x) format(x / 1e6, scientific = FALSE)) +
  theme  + xlab("Number of cells (million)") + ylab("") +scale_y_discrete(limits = levels(ncells$celltype))+coord_flip()+
  theme(axis.text.y=element_text(size=10), axis.title=element_text(size=10), aspect.ratio=0.6, axis.text.x=element_text(angle = 90, hjust = 0.95, vjust = 0.6, size=9))

## (add) stratified by sex  

ncells_by_sex <- metadata_all[metadata_all$cell_type %in% cells_to_keep,] %>% 
  dplyr::group_by(cell_type, Gender) %>% 
  dplyr::count() %>% dplyr::arrange(desc(n))
colnames(ncells_by_sex)[1] <- "celltype"

ncells_by_sex$Gender <- factor(ncells_by_sex$Gender,
                               levels = c('F', 'M'))
l2_order.in <- l2_order[l2_order%in%unique(ncells_by_sex$celltype)]
ncells_by_sex$celltype <- factor(ncells_by_sex$celltype,
                                 levels = l2_order.in)

ncells_l2_by_sex <-ggplot(ncells_by_sex, aes(y=celltype, x=n, fill=Gender)) +
  geom_bar(position="stack", stat="identity", alpha = 0.7)+ 
  scale_x_continuous(labels = function(x) format(x / 1e6, scientific = FALSE)) +
  theme  + 
  xlab("Number of cells (million)") + 
  ylab("") +
  scale_y_discrete(limits = levels(ncells_by_sex$celltype))+
  coord_flip()+
  scale_fill_manual(values=c(blue, red), name="Sex")+
  theme(axis.text.y=element_text(size=10), axis.title=element_text(size=10), aspect.ratio=0.5, axis.text.x=element_text(angle = 90, hjust = 0.95, vjust = 0.6, size=9))


#Number of donors per cell type 

ncells_donor <- metadata_all[metadata_all$cell_type %in% ncells$celltype,] %>%   distinct(cell_type, assignment) %>% dplyr::group_by(cell_type) %>% dplyr::count() %>% dplyr::arrange(desc(n))


ncells_donor_l2 <-ggplot(ncells_donor, aes(y=reorder(cell_type,-n), x=n)) + geom_bar(stat="identity", fill=alpha(blue, 0.7)) +
  theme  + xlab("Number of donors per cell type") + ylab("") +scale_y_discrete(limits = levels(ncells$cell_type))+coord_flip()+
  theme(axis.text.y=element_text(size=10), axis.title=element_text(size=10), aspect.ratio=0.6, axis.text.x=element_text(angle = 90, hjust = 0.95, vjust = 0.6, size=9))

## (add) stratified by sex  

ncells_donor_by_sex <- metadata_all[metadata_all$cell_type %in% ncells$celltype,] %>% 
  distinct(cell_type, assignment, Gender) %>% 
  dplyr::group_by(cell_type, Gender) %>% 
  dplyr::count() %>% dplyr::arrange(desc(n))

ncells_donor_by_sex$Gender <- factor(ncells_donor_by_sex$Gender,
                                     levels = c('F', 'M'))
l2_donors_order <- unique(ncells_donor_by_sex$cell_type)

ncells_donor_l2_by_sex <-ggplot(ncells_donor_by_sex, aes(y=reorder(cell_type,-n), x=n, fill=Gender)) + 
  geom_bar(position="stack", stat="identity", alpha = 0.7)+ 
  theme  + 
  xlab("Number of donors per cell type") + 
  ylab("") +
  scale_y_discrete(limits = levels(ncells_donor_by_sex$cell_type))+
  coord_flip()+
  scale_fill_manual(values=c(blue, red), name="Sex")+
  theme(axis.text.y=element_text(size=10), axis.title=element_text(size=10), aspect.ratio=0.5, axis.text.x=element_text(angle = 90, hjust = 0.95, vjust = 0.6, size=9))


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

## (add) stratified by sex  

ncells_donor_by_sex <- metadata_all %>% 
  dplyr::group_by(assignment, Gender) %>% 
  dplyr::count()
ncells_donor_by_sex$cells_donor <- "NCells_donor"

ncells_donor_by_sex.stats <- ncells_donor_by_sex %>%
  group_by(Gender) %>%
  summarize(median = median(n))
ncells_donor_by_sex$Gender <- factor(ncells_donor_by_sex$Gender,
                                     levels = c('F', 'M'))
ncells_donor_by_sex_p <- ggplot(ncells_donor_by_sex, aes(x=n, color = Gender)) + 
  geom_density(alpha= 0.5)+
  theme + 
  xlab("Number of cells per donor") + 
  ylab("Density") +
  geom_vline(data = ncells_donor_by_sex.stats, 
             aes(xintercept = median, color = Gender),
             linetype="dashed")+
  scale_color_manual(values=c(blue, red), name="Sex")+
  theme(axis.text.x=element_text(size=10), axis.text.y=element_text(size=10), axis.title=element_text(size=11), aspect.ratio=1)


#expressed genes 

# Figure S1A
pdf(paste0(fig_path, "/FigS1/S1A_AgeDistribution.pdf"), width = 3, height = 3)
age 
dev.off()

# Figure S1B
pdf(paste0(fig_path, "/FigS1/S1B_NumCells_l2.pdf"), width = 6, height = 3)
ncells_l2
dev.off()

# Figure S1C
pdf(paste0(fig_path, "/FigS1/S1C_NumDonor_Celltype.pdf"), width = 6, height = 3)
ncells_donor_l2 
dev.off()

# Figure S1D
pdf(paste0(fig_path, "/FigS1/S1D_SexDistribution.pdf"), width = 3, height = 3)
sex
dev.off()

# # Figure S1C
# pdf(paste0(fig_path, "/FigS1/S1C_NumCells_l1.pdf"), width = 3, height = 3)
# ncells_l1 
# dev.off()

# pdf(paste0(fig_path, "/FigS1/S1E_NumCells_Donor.pdf"), width = 3, height = 3)
# ncells_donor_p 
# dev.off()

# Figure S1E
pdf(paste0(fig_path, "/FigS1/S1E_AgeDistribution_by_sex.pdf"), width = 3, height = 3)
age_by_sex
dev.off()

# Figure S1F
pdf(paste0(fig_path, "/FigS1/S1F_NumCells_Donor_by_Sex.pdf"), width = 3, height = 3)
ncells_donor_by_sex_p 
dev.off()

# Figure S1G
pdf(paste0(fig_path, "/FigS1/S1G_NumCells_l2_by_Sex.pdf"), width = 6, height = 3)
ncells_l2_by_sex 
dev.off()

# Figure S1H
pdf(paste0(fig_path, "/FigS1/S1H_NumDonor_Celltype_by_sex.pdf"), width = 6, height = 3)
ncells_donor_l2_by_sex 
dev.off()



