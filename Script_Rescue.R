### Transcriptional predictors of rescue behaviour in ants

### Jaimes-Nino, Luisa Maria1, Bar, Adi2, Scharf, Inon2, Foitzik, Susanne1

### 1 Institute of Organismic and Molecular Evolution, Johannes Gutenberg University Mainz, Germany
### 2 School of Zoology, George S Wise Faculty of Life Sciences, Tel Aviv University, Tel Aviv, Israel


#### Brain size rescuers vs non-rescuers ####
library(glmmTMB)
library(emmeans)
library(DHARMa)
library(openxlsx)

size_rescue<- read.xlsx("~/Documents/04_Catag_rescue/MS_rescue/New_subm/JExpBio/Revision_JExpB_rescue/Supl_material_File1.xlsx", sheet = 2)


ggplot(data=size_rescue, aes(x=LengthB, y =WidthB, color = Rescue))+
geom_point()+
geom_smooth(method=lm)+
ylab("Width (mm)")+
xlab("Length (mm)")+
scale_color_manual(values=c("blue3", "red3"))+
theme_minimal()
sizeL.model <- glmmTMB(data = size_rescue, LengthB ~ Rescue + (1|Individual_id))
summary(sizeL.model)

##### Sensitivity analysis#####
HR_track<- read.xlsx("~/Documents/04_Catag_rescue/MS_rescue/New_subm/JExpBio/Revision_JExpB_rescue/Supl_material_File1.xlsx", sheet = 6)

HR_track$Average_Speed <- as.numeric(HR_track$Average_Speed)
HR_track$Average_Speed_i<- as.numeric(HR_track$Average_Speed_i)
HR_track$Average_Speed_o<- as.numeric(HR_track$Average_Speed_o)
table(HR_track$Resolution, HR_track$Rescue)
HR_track$ID <- paste(HR_track$Video, HR_track$Individual, sep = "_")
speed_o_lm <- glmmTMB(data= HR_track, Average_Speed_o ~ Resolution +  Rescue + (1|ID), family = Gamma)

simulateResiduals(speed_o_lm, plot = T)
summary(speed_o_lm)
emm.resc <- emmeans(speed_o_lm, "Rescue")
str(emm.resc)
summary(emm.resc, infer = TRUE, type = "response")


speed_i_lm <- glmmTMB(data= HR_track, Average_Speed_i ~ Resolution + Rescue + (1|ID), family =  Gamma)
simulateResiduals(speed_i_lm, plot = T)

summary(speed_i_lm)


emm.resc <- emmeans(speed_i_lm, "Rescue")

summary(emm.resc, infer = TRUE, type = "response")

# Average speed

speed_lm <- glmmTMB(data= HR_track, log(Average_Speed) ~ Resolution + Rescue + (1|ID))
simulateResiduals(speed_lm, plot = T)
summary(speed_lm)

#### General tracking analysis ####
library(patchwork)
library(ggplot2)
library(ggbeeswarm)
library(dplyr)

general_track<- read.xlsx("~/Documents/04_Catag_rescue/MS_rescue/New_subm/JExpBio/Revision_JExpB_rescue/Supl_material_File1.xlsx", sheet = 4)

gen_rescuers <- subset(general_track, Rescue == "Yes")
gen_nonrescuers <- subset(general_track, Rescue == "No")
gen_rescue_track <- rbind(gen_rescuers, gen_nonrescuers)
gen_rescue_track$Average_Speed <- as.numeric(gen_rescue_track$Average_Speed)

gen_plotA <- ggplot(data=gen_rescue_track, aes(x=Rescue, y=Average_Speed, color = Rescue))+
  geom_beeswarm(priority = "none")+
  stat_summary(fun.y= mean, fun.ymin=mean, fun.ymax=mean, geom="crossbar", width=0.3, linewidth =0.3)+
  annotate("text", color = "black",label = "n.s", x = 1.45, y = 1.65) +
  scale_color_manual(values=c("blue3", "red3"))+
  theme_classic()+
  ylab("Speed (cm/sec)")+
  xlab("Rescue")+
  theme(legend.position="none")

speed_gen_lm <- glmmTMB(data= gen_rescue_track, Average_Speed ~ Rescue + (1|Video))
simulateResiduals(speed_gen_lm, plot = T)
summary(speed_gen_lm)
# Family: gaussian  ( identity )
# Formula:          Average_Speed ~ Rescue + (1 | Video)
# Data: gen_rescue_track
# 
# AIC       BIC    logLik -2*log(L)  df.resid 
# 26.6      32.9      -9.3      18.6        31 
# 
# Random effects:
#   
#   Conditional model:
#   Groups   Name        Variance Std.Dev.
# Video    (Intercept) 0.03561  0.1887  
# Residual             0.07880  0.2807  
# Number of obs: 35, groups:  Video, 7
# 
# Dispersion estimate for gaussian family (sigma^2): 0.0788 
# 
# Conditional model:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)  0.93005    0.10361   8.977   <2e-16 ***
#   RescueYes    0.01454    0.09967   0.146    0.884     


circle_track<- read.xlsx("~/Documents/04_Catag_rescue/MS_rescue/New_subm/JExpBio/Revision_JExpB_rescue/Supl_material_File1.xlsx", sheet = 5)
colnames(circle_track)
rescuers <- subset(circle_track, Rescue == "Yes")
nonrescuers <- subset(circle_track, Rescue == "No")

rescue_track <- rbind(rescuers, nonrescuers)
rescue_track$Latency_i<- as.numeric(rescue_track$Latency_i)
rescue_track$Latency_o<- as.numeric(rescue_track$Latency_o)
rescue_track$Average_Speed_i<- as.numeric(rescue_track$Average_Speed_i)
rescue_track$Average_Speed_o<- as.numeric(rescue_track$Average_Speed_o)

table( rescue_track$Rescue)
#rescue_track$Rescue <- factor(rescue_track$Rescue, levels = c("Yes", "No"))


inner_plot <- ggplot(data=rescue_track, aes(x=Rescue, y=Average_Speed_i, colour = Rescue))+
  geom_beeswarm(priority = "none")+
  stat_summary(fun.y= mean, fun.ymin=mean, fun.ymax=mean, geom="crossbar", width=0.3, linewidth =0.3)+
  scale_color_manual(values=c("blue3", "red3"))+
  theme_classic()+
  annotate("text", color = "black",label = "**", x = 1.45, y = 3, size = 7) +
  ylab("Speed inner circle (cm/sec)")+
  xlab("Rescue")+
  theme(legend.position="none")

speed_i_lm <- glmmTMB(data= rescue_track, log(Average_Speed_i) ~ Rescue + (1|Video))
simulateResiduals(speed_i_lm, plot = T)
summary(speed_i_lm)
# Family: gaussian  ( identity )
# Formula:          log(Average_Speed_i) ~ Rescue + (1 | Video)
# Data: rescue_track
# 
# AIC       BIC    logLik -2*log(L)  df.resid 
# 48.0      53.2     -20.0      40.0        23 

# Random effects:
#   
#   Conditional model:
#   Groups   Name        Variance Std.Dev.
# Video    (Intercept) 0.3268   0.5717  
# Residual             0.1430   0.3782  
# Number of obs: 27, groups:  Video, 7
# 
# Dispersion estimate for gaussian family (sigma^2): 0.143 
# 
# Conditional model:
#   Estimate Std. Error z value Pr(>|z|)   
# (Intercept)  -0.5094     0.2649  -1.923  0.05446 . 
# RescueYes    -0.5776     0.1776  -3.252  0.00114 **

outer_plot <- ggplot(data=rescue_track, aes(x=Rescue, y=Average_Speed_o, colour = Rescue))+
  geom_beeswarm(priority = "none")+
  stat_summary(fun.y= mean, fun.ymin=mean, fun.ymax=mean, geom="crossbar", width=0.3, linewidth =0.3)+
  annotate("text", color = "black",label = "***", x = 1.45, y = 3, size = 7) +
  scale_color_manual(values=c("blue3", "red3"))+
  theme_classic()+
  ylab("Speed outer circle (cm/sec)")+
  xlab("Rescue")+
  theme(legend.position="none")

speed_o_lm <- glmmTMB(data= rescue_track, log(Average_Speed_o) ~ Rescue + (1|Video))
simulateResiduals(speed_o_lm, plot = T)
summary(speed_o_lm)
# Family: gaussian  ( identity )
# Formula:          log(Average_Speed_o) ~ Rescue + (1 | Video)
# Data: rescue_track
# 
# AIC       BIC    logLik -2*log(L)  df.resid 
# 52.5      58.2     -22.3      44.5        27 
# 
# Random effects:
#   
#   Conditional model:
#   Groups   Name        Variance Std.Dev.
# Video    (Intercept) 0.06873  0.2622  
# Residual             0.19971  0.4469  
# Number of obs: 31, groups:  Video, 7
# 
# Dispersion estimate for gaussian family (sigma^2):  0.2 
# 
# Conditional model:
          #   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)   0.1223     0.1726   0.709    0.478    
# RescueYes    -0.6902     0.1773  -3.893 9.88e-05 ***
coefficients(speed_o_lm)
exp(0.2329)
exp(-0.6902)

speed_plotC <- gen_plotA+ outer_plot + inner_plot + plot_annotation(tag_levels = 'A')

inner_plot <- ggplot(data=rescue_track, aes(x=Rescue, y=Exploration_relative_value_i, colour = Rescue))+
  geom_beeswarm(priority = "none")+
  stat_summary(fun.y= mean, fun.ymin=mean, fun.ymax=mean, geom="crossbar", width=0.3, linewidth =0.3)+
  scale_color_manual(values=c("blue3", "red3"))+
  annotate("text", color = "black",label = "***", x = 1.45, y = 1, size = 7) +
  theme_classic()+
  ylab("Relative exploration inner circle")+
  xlab("Rescue")+
  theme(legend.position="none")

explor_i_lm <- glmmTMB(data= rescue_track, Exploration_relative_value_i ~ Rescue + (1|Video))
testDispersion(explor_i_lm ) # ok

simulateResiduals(explor_i_lm, plot = T)
summary(explor_i_lm)
# Conditional model:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)  0.30077    0.07447   4.039 5.38e-05 ***
#   RescueYes    0.48710    0.08848   5.505 3.68e-08 ***

table(rescue_track$Rescue)
outer_plot <- ggplot(data=rescue_track, aes(x=Rescue, y=Exploration_relative_value_o, colour = Rescue))+
  geom_beeswarm(priority = "none")+
  stat_summary(fun.y= mean, fun.ymin=mean, fun.ymax=mean, geom="crossbar", width=0.3, linewidth =0.3)+
  scale_color_manual(values=c("blue3", "red3"))+
  theme_classic()+
  annotate("text", color = "black",label = "**", x = 1.45, y = 1, size = 7) +
  ylab("Relative exploration outer circle")+
  xlab("Rescue")+
  theme(legend.position="none")
explor_o_lm <- glmmTMB(data= rescue_track, Exploration_relative_value_o ~ Rescue + (1|Video))
testDispersion(explor_o_lm ) # ok
simulateResiduals(explor_o_lm, plot = T)
summary(explor_o_lm)

# Conditional model:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)  0.31087    0.08519   3.649 0.000263 ***
#   RescueYes    0.23932    0.08581   2.789 0.005285 ** 
generalExpl_plot <- ggplot(data=gen_rescue_track, aes(x=Rescue, y=Exploration_relative_value, colour = Rescue))+
  geom_beeswarm(priority = "none")+
  stat_summary(fun.y= mean, fun.ymin=mean, fun.ymax=mean, geom="crossbar", width=0.3, linewidth =0.3)+
  scale_color_manual(values=c("blue3", "red3"))+
  theme_classic()+
  annotate("text", color = "black",label = "n.s", x = 1.45, y = 1) +
  ylab("Relative exploration arena")+
  xlab("Rescue")+
  theme(legend.position="none")


explor_lm <- glmmTMB(data= gen_rescue_track, Exploration_relative_value ~ Rescue + (1|Video))
testDispersion(explor_lm ) # ok
simulateResiduals(explor_lm, plot = T)
summary(explor_lm)
# Conditional model:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)  0.41451    0.06695   6.191 5.97e-10 ***
#  RescueYes    0.03151    0.06690   0.471    0.638   
Exploration_plot_D <- generalExpl_plot + outer_plot + inner_plot  + plot_annotation(tag_levels = 'A')



inner_lat_plot <- ggplot(data=rescue_track, aes(x=Rescue, y=Latency_i, colour = Rescue))+
  geom_beeswarm(priority = "none")+
  stat_summary(fun.y= mean, fun.ymin=mean, fun.ymax=mean, geom="crossbar", width=0.3, linewidth =0.3)+
  scale_color_manual(values=c("blue3", "red3"))+
  annotate("text", color = "black",label = "**", x = 1.45, y = 850, size = 7) +
  theme_classic()+
  ylab("Latency to reach inner circle (s)")+
  xlab("Rescue")+
  theme(legend.position="none")

latenc_i_lm <- glmmTMB(data= rescue_track, Latency_i ~ Rescue + (1|Video))
simulateResiduals(latenc_i_lm, plot = T)
summary(latenc_i_lm)
# Conditional model:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)   594.66      89.15   6.671 2.55e-11 ***
#   RescueYes    -254.03      90.10  -2.819  0.00481 ** 
outer_lat_plot <- ggplot(data=rescue_track, aes(x=Rescue, y=Latency_o, colour = Rescue))+
  geom_beeswarm(priority = "none")+
  stat_summary(fun.y= mean, fun.ymin=mean, fun.ymax=mean, geom="crossbar", width=0.3, linewidth =0.3)+
  scale_color_manual(values=c("blue3", "red3"))+
  annotate("text", color = "black",label = "**", x = 1.45, y = 850, size = 7) +
  theme_classic()+
  ylab("Latency to reach outer circle (s)")+
  xlab("Rescue")+
  theme(legend.position="none")

latenc_o_lm <- glmmTMB(data= rescue_track, Latency_o ~ Rescue + (1|Video))
simulateResiduals(latenc_o_lm, plot = T)
summary(latenc_o_lm)
# Conditional model:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)   484.45      75.43   6.422 1.34e-10 ***
#   RescueYes    -189.42      72.78  -2.603  0.00925 ** 
inner_lat_plot + outer_lat_plot

gen_rescue_track<- gen_rescue_track %>%
  mutate(across(c(3:14), as.numeric))
travDist_plot <- ggplot(data=gen_rescue_track, aes(x=Rescue, y=Traveled_Dist, colour = Rescue))+
  geom_beeswarm(priority = "none")+
  stat_summary(fun.y= mean, fun.ymin=mean, fun.ymax=mean, geom="crossbar", width=0.3, linewidth =0.3)+
  scale_color_manual(values=c("blue3", "red3"))+
  annotate("text", color = "black",label = "n.s", x = 1.45, y = 1150) +
  theme_classic()+
  ylab("Travelled distance (cm)")+
  xlab("Rescue")+
  theme(legend.position="none")
dist_lm <- glmmTMB(data= gen_rescue_track, log(Traveled_Dist) ~ Rescue + (1|Video))
simulateResiduals(dist_lm, plot = T)
summary(dist_lm)
# Conditional model:
# Intercept)   4.9692     0.2533  19.619   <2e-16 ***
#   RescueYes     0.1112     0.3078   0.361    0.718  
((generalExpl_plot + outer_plot + inner_plot)  /
  (travDist_plot + outer_lat_plot + inner_lat_plot ) )+ plot_annotation(tag_levels = 'A')

######################## GOterms #######################
library(openxlsx)
library(ggplot2)
library(Rgraphviz)
library(topGO)
library("RColorBrewer")
library(tagcloud)


#### MBs rescuers vs non-rescuers ####
MB_up_genes <- as.data.frame(c("gene_3661",  "gene_3662",  "gene_3426",  "gene_1607",  "gene_2555", "gene_1243" , "gene_930" ,  "gene_5298" , "gene_6089",  "gene_5245",  "gene_13149", "gene_9895",  "gene_9900",  "gene_8499",  "gene_7997"))

colnames(MB_up_genes) <- c("genes") 

#### GOterms Upregulated when prone to rescue
#we start with an enrichment analysis
#load the universe reference file you created above using bash
universe <- readMappings("~/Documents/Catagl_cognition/From_Eyal/Cnig_gn3.1.fasta.transcripts.fasta.trinotate.report.topGOinput.txt")
# add newly found annotations:
universe$gene_3426 <- c("GO:0008061", "GO:0005576")
universe$gene_5298 <- c("GO:0005549")
universe$gene_1607 <- c("GO:0035176","GO:0005549", "GO:0005615")
universe$gene_9900 <- c("GO:0006955","GO:0005164", "GO:0005515", "GO:0016020", "GO:0048018", "GO:0005615")
universe$gene_5301 <- c("GO:0035176","GO:0005549", "GO:0005615")

# generate list of the protein IDs
universe_genelist <- names(universe)

# generate a reference list of the protein IDs
#testset_genelist_idown <- readLines(LD_down_genes)
#determine matches between both datasets
gene_list_iup <- factor(as.integer(universe_genelist %in% MB_up_genes$genes))
names(gene_list_iup) <- universe_genelist
#set the ontology type you want to investigate, we first choose "Biological Processes"=BP, '93Cellular Components'94=CC, '93Molecular Function'94=MF
ontology_type = "BP" # or CC or MF
GO_data_iup <- new("topGOdata", description="GO_Enrichment", ontology=ontology_type, allGenes=gene_list_iup, annot = annFUN.gene2GO, gene2GO=universe)

result_topGO_w01_iup <- runTest(GO_data_iup, algorithm = "weight01", statistic = "fisher")
result_table_w01_iup_BP <- GenTable(GO_data_iup, Fisher = result_topGO_w01_iup, orderBy = "Fisher", ranksOf = "Fisher", topNodes = numSigGenes(GO_data_iup))
result_table_w01_iup_BP
#### Graphics ####
result_table_w01_iup_BP_sig <- result_table_w01_iup_BP[which(result_table_w01_iup_BP$Fisher<0.05),]
setwd("~/Documents/Catag_rescue/")
#write.csv(result_table_w01_idown_BP_sig, file=paste0("topGO_result_weight01_significant_OL_Corrprop_BP_down_cooks.csv")
colors <- colorRampPalette( brewer.pal( 12, "Paired" ) )( length(result_table_w01_iup_BP_sig$Term))
#Create the tag cloud
tagcloud(strmultline(result_table_w01_iup_BP_sig$Term), weights = - log(as.numeric(result_table_w01_iup_BP_sig$Fisher)),col=colors)



## DESeq2 - mRNA sequencing data from Cat niger. 7 samples of rescuers and non-rescuers. The brain tissue was sequenced by separating the optic lobes (OL), the mushroom bodies (MB), the rest of the brain tissue including antennal lobes and central complex (CC), and the antennae (AA). 

tissues <- read.xlsx("~/Documents/Catag_rescue/Sequencing_rescue/Gene_counts/Rescue_rawcounts_perGene.xlsx", colNames = T)
tissue_counts <- tissues[-c(1:4),]
colnames(tissue_counts)[1] <- c("Gene_ID")
rownames(tissue_counts) <- tissue_counts$Gene_ID
tissue_counts <- tissue_counts[,-c(1)]
colnames(tissue_counts)
# [1] "Ant_A_727" "Ant_B_724" "Ant_C_723" "Ant_D_727" "Ant_E_718" "Ant_F_718" "Ant_G_720" "Ant_H_720"
# [9] "Ant_I_723" "Ant_J_724" "Ant_K_725" "Ant_L_725" "Ant_M_726" "Ant_N_726" "CC_A_727"  "CC_B_724" 
# [17] "CC_C_723"  "CC_D_727"  "CC_F_718"  "CC_G_720"  "CC_H_720"  "CC_I_723"  "CC_J_724"  "CC_K_725" 
# [25] "CC_L_725"  "CC_M_726"  "CC_N_726"  "MB_A_727"  "MB_C_723"  "MB_D_727"  "MB_E_718"  "MB_F_718" 
# [33] "MB_G_720"  "MB_H_720"  "MB_I_723"  "MB_J_724"  "MB_K_725"  "MB_L_725"  "MB_M_726"  "MB_N_726" 
# [41] "OL_A_727"  "OL_B_724"  "OL_C_723"  "OL_E_718"  "OL_F_718"  "OL_G_720"  "OL_H_720"  "OL_I_723" 
# [49] "OL_J_724"  "OL_K_725"  "OL_L_725"  "OL_M_726"  "OL_N_726" 

# First check differences among tissues to highlight any outliers, or if the dissection failed. 

design_rescue<- read.xlsx("~/Documents/Catag_rescue/Samples_rescue.xlsx", sheet = 6)
head(design_rescue)
# Sample_ID Tissue Pool Colony_ID Rescue
# 1 Ant_A_727    Ant    A       727     No
# 2 Ant_B_724    Ant    B       724    Yes
# 3 Ant_C_723    Ant    C       723    Yes
# 4 Ant_D_727    Ant    D       727    Yes
# 5 Ant_E_718    Ant    E       718     No
# 6 Ant_F_718    Ant    F       718    Yes
design_rescue[,1]
# [1] "Ant_A_727" "Ant_B_724" "Ant_C_723" "Ant_D_727" "Ant_E_718" "Ant_F_718" "Ant_G_720" "Ant_H_720"
# [9] "Ant_I_723" "Ant_J_724" "Ant_K_725" "Ant_L_725" "Ant_M_726" "Ant_N_726" "CC_A_727"  "CC_B_724" 
# [17] "CC_C_723"  "CC_D_727"  "CC_F_718"  "CC_G_720"  "CC_H_720"  "CC_I_723"  "CC_J_724"  "CC_K_725" 
# [25] "CC_L_725"  "CC_M_726"  "CC_N_726"  "MB_A_727"  "MB_C_723"  "MB_D_727"  "MB_E_718"  "MB_F_718" 
# [33] "MB_G_720"  "MB_H_720"  "MB_I_723"  "MB_J_724"  "MB_K_725"  "MB_L_725"  "MB_M_726"  "MB_N_726" 
# [41] "OL_A_727"  "OL_B_724"  "OL_C_723"  "OL_E_718"  "OL_F_718"  "OL_G_720"  "OL_H_720"  "OL_I_723" 
# [49] "OL_J_724"  "OL_K_725"  "OL_L_725"  "OL_M_726"  "OL_N_726" 
design_rescue$Colony_ID <- as.factor(design_rescue$Colony_ID)
design_rescue$Rescue <- as.factor(design_rescue$Rescue)
design_rescue$Tissue <- as.factor(design_rescue$Tissue)

dds_alltissues <- DESeqDataSetFromMatrix(countData = tissue_counts, colData = design_rescue , design = ~ Tissue + Rescue)
# converting counts to integer mode
keep <- rowSums(counts(dds_alltissues) >= 10) >= 6
dds_alltissues <- dds_alltissues[keep,]

# Visualization - let's first have a look at the data in a PCA plot. For this the data need to be transformed, for which several methods are described within the manual. Choose one, e.g. varianceStabilizingTransformation and carry out the command on the dds object. Then run the command plotPCA on the new transformed object, naming the grouping variables with intgroup
dds_vst_allTissues <- varianceStabilizingTransformation(dds_alltissues)
pcaData_allTissues <-  plotPCA(dds_vst_allTissues, intgroup = c("Tissue", "Colony_ID"), returnData = T)
# using ntop=500 top features by variance

pca_allTissues <- ggplot(pcaData_allTissues, aes(x = PC1, y = PC2, label =name, color = Tissue,shape =Colony_ID )) + 
  geom_point(size=3) + 
  # geom_label_repel(label.size = 0)+
  # geom_text(aes(label=name), vjust=3, size=4,nudge_x = 3, nudge_y = 3)+
  xlab("PC1: 85% variance")+
  ylab("PC2 : 9% variance")+
  ggtitle("PCA all Tissues Rescue")+
  scale_shape_manual(values = c(0, 5, 10, 14, 1,8,17))+
  theme_minimal()+
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour="black"),
        legend.position = "bottom",
        legend.text=element_text(size=12),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 12))
pca_allTissues

# Wald test in all tissues, but then test one to one each tissue, better all tissues that have DEGs with an LRT test
# dds_alltissues <- DESeq(dds_alltissues)
# Wald_rescue_results <- results(dds_alltissues, contrast = c("Rescue", "Yes", "No"), alpha = 0.05)
# Wald_rescue_results_sig <- subset(Wald_rescue_results , padj < 0.05)
# Wald_rescue_results_sig # DataFrame with 
# A reduced model without the rescue variable
dds_allTissue_LTR_rescue <- DESeq(dds_alltissues, test="LRT", full= ~  Tissue + Rescue , reduced = ~Tissue)

LRT_results_allTissue_rescue <- results(dds_allTissue_LTR_rescue)
#write.csv(LRT_results_allTissue_rescue, "DEGs_alltissues.csv")
LRT_results_allTissue_rescue_sig <- subset(LRT_results_allTissue_rescue , padj < 0.05)
LRT_results_allTissue_rescue_sig # DataFrame with 9 rows and 6 columns
as.data.frame(LRT_results_allTissue_rescue_sig)
# baseMean log2FoldChange      lfcSE     stat       pvalue        padj
# gene_930    278.50568      2.8281217 0.47015841 27.54054 1.538359e-07 0.001583740
# gene_2555   117.04428      2.7623696 0.53417219 23.27896 1.401249e-06 0.004808619
# gene_5245  2322.00176      1.0977149 0.22555155 22.50527 2.095677e-06 0.005393748
#non-shared with MBs genes
# gene_1997    34.03176      0.7625534 0.17878296 17.96929 2.244977e-05 0.029849548
# gene_13318  789.81013      0.7336380 0.17270608 17.68300 2.609480e-05 0.029849548
# gene_13632  631.05258     -0.2551454 0.05794553 19.29975 1.117209e-05 0.023003325
# gene_14429  450.17155      0.3543307 0.08393271 17.79513 2.460122e-05 0.029849548
# gene_14458 2677.47760      0.5531638 0.10852223 25.92235 3.554297e-07 0.001829574
# gene_7603    95.96581      0.5167914 0.11947306 18.71925 1.514459e-05 0.025985597
LRT_results_allTissue_rescue_sig@rownames
alltissues_genes <- c("gene_1997",  "gene_2555",  "gene_930",   "gene_5245",  "gene_13318", "gene_13632", "gene_14429", "gene_14458", "gene_7603") 

LRT_results_allTissue_rescue_sig$gene <- rownames(LRT_results_allTissue_rescue_sig)
df_allT <- lapply(LRT_results_allTissue_rescue_sig$gene, (x) {
  y <- plotCounts(dds_allTissue_LTR_rescue,  returnData=TRUE, x, intgroup = "Rescue")
  y$feature <- x
  return(y)
})

df_allTs <- do.call(rbind, df_allT)
df_allTs<- df_allTs %>%
  group_by(feature) %>%
  mutate(y.position=max(count))


padj_allTs <- as.data.frame(LRT_results_allTissue_rescue_sig[[6]])
padj_allTs$feature <- LRT_results_allTissue_rescue_sig[[7]]
padj_allTs$`LRT_results_allTissue_rescue_sig[[6]]` <- lapply(padj_allTs$`LRT_results_allTissue_rescue_sig[[6]]`, signif, digits=3)

large_padj_allTs <- merge(padj_allTs, df_allTs, by = "feature", all= FALSE)
padj_allTs_small <- large_padj_allTs[!duplicated(large_padj_allTs$feature), ]

plot_a_list <- function(master_list_with_plots, no_of_rows, no_of_cols) {
  
  patchwork::wrap_plots(master_list_with_plots,
                        nrow = no_of_rows, ncol = no_of_cols)
}


# apply seperately for each gene showing the mean
p_allT <- lapply(LRT_results_allTissue_rescue_sig$gene, function(gene) {
  ggplot(df_allTs[df_allTs[, "feature"]==gene,]) +
    geom_beeswarm(aes(x=Rescue, y=count, color= Rescue),priority = "none")+
    stat_summary(aes(x=Rescue, y=count), fun.y= mean, fun.ymin=mean, fun.ymax=mean, geom="crossbar", width=0.3, linewidth =0.3)+
    #  annotate("text", label= "Padj=", x= 1,y=(padj_MBs_small[match( gene, padj_MBs_small$feature),5]))+
    geom_text(data = as.data.frame(padj_allTs_small),
              aes(label = (padj_allTs_small[match( gene,feature),2]),
                  x = 1.45, y = (padj_allTs_small[match( gene, feature),5]))) +
    ylab("Normalized count")+
    scale_color_manual(values=c("blue3", "red3"))+
    theme_classic()+
    ggtitle(paste(gene))+
    theme(axis.title.y=element_blank(),
          legend.position="none")
})


# finally print your plots
plot_a_list(p_allT, 3, 3)

##### 15 genes in MB, how are they in other tissues? ###
MB_genes

df_allT_15MBs <- lapply(MB_genes, (x) {
  y <- plotCounts(dds_allTissue_LTR_rescue,  returnData=TRUE, x, intgroup =c("Rescue","Tissue"))
  y$feature <- x
  return(y)
})

p_allT_15 <- lapply(df_allT_15MBs, function(df) {
  ggplot(df) +
    geom_beeswarm(aes(x=Tissue, y=count, color = Rescue),priority = "none")+
    stat_summary(aes(x=Tissue, y=count, color = Rescue), fun.y= mean, fun.ymin=mean, fun.ymax=mean, geom="crossbar", width=0.3, linewidth =0.3)+
    ylab("Normalized count")+
    scale_color_manual(values=c("blue3", "red3"))+
    theme_classic()+
    ggtitle(paste(df$feature))+
    theme(axis.title.y=element_blank(),
          legend.position="none")
})
plot_a_list(p_allT_15, 5, 3)

gene_14458 <- plotCounts(dds_allTissue_LTR_rescue,  returnData=TRUE, "gene_14458", intgroup =c("Rescue","Tissue")) 
ggplot(data = gene_14458) +
  geom_beeswarm(aes(x=Tissue, y=count, color = Rescue),priority = "none")+
  stat_summary(aes(x=Tissue, y=count, color = Rescue), fun.y= mean, fun.ymin=mean, fun.ymax=mean, geom="crossbar", width=0.3, linewidth =0.3)+
  ylab("Normalized count")+
  scale_color_manual(values=c("blue3", "red3"))+
  theme_classic()+
  #  ggtitle(paste(df$feature))+
  theme(legend.position="right")

# A reduced model without the tissue variable
dds_allTissue_LTR_tissue <- DESeq(dds_alltissues, test="LRT", full= ~ Tissue +Rescue, reduced = ~Rescue)

LRT_results_allTissue_tissue <- results(dds_allTissue_LTR_tissue)
LRT_results_allTissue_tissue_sig <- subset(LRT_results_allTissue_tissue , padj < 0.05)
LRT_results_allTissue_tissue_sig # DataFrame with 9913 rows and 6 columns
LRT_results_allTissue_tissue_sig$gene <- rownames(LRT_results_allTissue_tissue_sig)
LRT_results_allTissue_frame <- as.data.frame(LRT_results_allTissue_tissue_sig)
#write.xlsx(LRT_results_allTissue_frame, "DEGs_alltissues.xlsx")

# DESEq per tissue 

##### Optic lobes ####
design_OL <- subset(design_rescue, Tissue =="OL")
OL_counts <- tissue_counts[,colnames(tissue_counts) %in% design_OL$Sample_ID]
dds_OL <- DESeqDataSetFromMatrix(countData = OL_counts, colData = design_OL , design = ~ Colony_ID + Rescue)
# converting counts to integer mode
keep_OL <- rowSums(counts(dds_OL) >= 10) >= 6
dds_OL <- dds_OL[keep_OL,]

# A reduced model without the rescue factor
dds_OL_ltr_rescue <- DESeq(dds_OL, test="LRT", full= ~ Colony_ID + Rescue , reduced = ~ Colony_ID)

LRT_results_rescue_OL <- results(dds_OL_ltr_rescue)
#write.csv(LRT_results_rescue_OL, "DEGs_OL.csv")
LRT_results_rescue_OL_sig <- subset(LRT_results_rescue_OL, padj < 0.05)
LRT_results_rescue_OL_sig # 1 gene "affected"
# log2 fold change (MLE): Rescue Yes vs No 
# LRT p-value: '~ Colony_ID + Rescue' vs '~ Colony_ID' 
# DataFrame with 1 row and 6 columns
# baseMean log2FoldChange     lfcSE      stat     pvalue        padj
# <numeric>      <numeric> <numeric> <numeric>  <numeric>   <numeric>
#   gene_5301   399.255         1.3244  0.240348    28.644 8.6982e-08 0.000842073
opticlobes_genes <- c("gene_5301")

LRT_results_rescue_OL_sig$gene <- rownames(LRT_results_rescue_OL_sig)

y <- plotCounts(dds_OL_ltr_rescue ,  returnData=TRUE,opticlobes_genes , intgroup = "Rescue")
y$feature <- c("gene_5301")
y.position=max(y$count)

padj_OLs <- as.data.frame(LRT_results_rescue_OL_sig[[6]])
padj_OLs$feature <- LRT_results_rescue_OL_sig[[7]]
padj_OLs$`LRT_results_rescue_OL_sig[[6]]` <- lapply(padj_OLs$`LRT_results_rescue_OL_sig[[6]]`, signif, digits=3)

large_padj_OLs <- merge(padj_OLs, y, by = "feature", all= FALSE)
padj_OLs_small <- large_padj_OLs[!duplicated(large_padj_OLs$feature), ]

plot_a_list <- function(master_list_with_plots, no_of_rows, no_of_cols) {
  
  patchwork::wrap_plots(master_list_with_plots,
                        nrow = no_of_rows, ncol = no_of_cols)
}

# apply seperately for each gene showing the mean
p_OL <- ggplot(y) +
  geom_beeswarm(aes(x=Rescue, y=count, color= Rescue),priority = "none")+
  stat_summary(aes(x=Rescue, y=count), fun.y= mean, fun.ymin=mean, fun.ymax=mean, geom="crossbar", width=0.3, linewidth =0.3)+
  #  annotate("text", label= "Padj=", x= 1,y=(padj_MBs_small[match( gene, padj_MBs_small$feature),5]))+
  geom_text(data = as.data.frame(padj_OLs_small),
            aes(label = (padj_OLs_small[1,2]),
                x = 1.45, y = y.position)) +
  ylab("Normalized count")+
  scale_color_manual(values=c("blue3", "red3"))+
  theme_classic()+
  ggtitle(paste("gene_5301"))+
  theme(axis.title.y=element_blank(),
        legend.position="none")


# finally print your plots
plot_a_list(p_MB, 4, 4)

##### Mushroom Bodies ####
design_MB <- subset(design_rescue, Tissue =="MB")
MB_counts <- tissue_counts[,colnames(tissue_counts) %in% design_MB$Sample_ID]
dds_MB <- DESeqDataSetFromMatrix(countData = MB_counts, colData = design_MB , design = ~ Colony_ID + Rescue)
# converting counts to integer mode
keep_MB <- rowSums(counts(dds_MB) >= 10) >= 6
dds_MB <- dds_MB[keep_MB,]
# A reduced model without the rescue factor
dds_MB_ltr_rescue <- DESeq(dds_MB, test="LRT", full= ~ Colony_ID + Rescue , reduced = ~ Colony_ID)

LRT_results_rescue_MB <- results(dds_MB_ltr_rescue)
#write.csv(LRT_results_rescue_MB, "DEGs_MB.csv")
LRT_results_rescue_MB_sig <- subset(LRT_results_rescue_MB, padj < 0.05)
LRT_results_rescue_MB_sig 
# log2 fold change (MLE): Rescue Yes vs No 
# LRT p-value: '~ Colony_ID + Rescue' vs '~ Colony_ID' 
# DataFrame with 15 rows and 6 columns
# baseMean log2FoldChange     lfcSE      stat      pvalue        padj
# <numeric>      <numeric> <numeric> <numeric>   <numeric>   <numeric>
#   gene_3661   156.4731       2.496080  0.426528   28.9972 7.24820e-08 9.90312e-05
# gene_3662    98.3451       1.977841  0.435669   17.4028 3.02374e-05 2.10276e-02
# gene_3426   250.8087       0.962072  0.225606   17.3690 3.07807e-05 2.10276e-02
# gene_1607   167.2192       1.373933  0.242181   30.4875 3.36027e-08 6.42752e-05
# gene_2555    56.0721       2.647648  0.446123   30.8862 2.73609e-08 6.42752e-05
# ...              ...            ...       ...       ...         ...         ...
# gene_13149 1019.5838       0.505851  0.114291   19.4021 1.05890e-05 9.43983e-03
# gene_9895    64.7391       2.145827  0.369639   31.1385 2.40262e-08 6.42752e-05
# gene_9900   169.3375       0.948923  0.234393   15.7326 7.29542e-05 4.65156e-02
# gene_8499    46.5866       1.933660  0.399757   21.4861 3.56396e-06 4.26071e-03
# gene_7997   137.1182       1.403238  0.301844   19.9517 7.94241e-06 8.44014e-03

LRT_results_rescue_MB_sig@rownames
MB_genes <-c( "gene_3661",  "gene_3662",  "gene_3426",  "gene_1607",  "gene_2555", "gene_1243" , "gene_930" ,  "gene_5298" , "gene_6089",  "gene_5245",  "gene_13149", "gene_9895",  "gene_9900",  "gene_8499",  "gene_7997") 
library(ggVennDiagram)
list_genes <- list(alltissues_genes, opticlobes_genes, MB_genes)
Venn_genes <-   ggVennDiagram(list_genes, label = c("count"), label_alpha = 0, category.names = c("All tissues", "OL", "MB"))+
  scale_fill_gradient(low = "#FFFFFF", high = "#FF0033")


dds_vst_MB <- varianceStabilizingTransformation(dds_MB)
pcaData_MB <-  plotPCA(dds_vst_MB, intgroup = c("Rescue", "Colony_ID", "Pool"), returnData = T)

# using ntop=500 top features by variance

pca_MB <- ggplot(pcaData_MB, aes(x = PC1, y = PC2, color = Rescue, shape =  Colony_ID)) + 
  geom_point(size=3) + 
  # geom_label_repel(label.size = 0)+
  xlab("PC1: 33% variance")+
  ylab("PC2 : 19% variance")+
  ggtitle("PCA MB Rescue")+
  #   scale_color_manual(values = c("#C77CFF", "#00BFC4", "#7CAE00", "#F8766D"))+
  theme_minimal()+
  scale_shape_manual(values = c(0, 5, 10, 14, 1,8,17))+
  scale_color_manual(values=c("blue3", "red3"))+
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour="black"),
        legend.position = "bottom",
        legend.text=element_text(size=12),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 12))
pca_MB

Venn_genes + pca_MB + plot_layout(widths = c(1,2))


plotMA(LRT_results_rescue_MB, ylim=c(-2,2), alpha = 0.05) # logarithmic fold change log2(yes/no).
# resLFC <- lfcShrink(dds_MB_ltr_rescue, coef="Rescue_Yes_vs_No", type="apeglm")
# plotMA(resLFC, ylim=c(-2,2), alpha = 0.05) # shrunken log2 fold changes, which remove the noise associated with log2 fold changes from low count genes without requiring arbitrary filtering thresholds.

#### Plot MB DEGs count ###

#The counts are which normalizes counts by the estimated size factors (or normalization factors if these were used) and adds a pseudocount of 1/2 to allow for log scale plotting. 

#plotCounts(dds_MB_ltr_rescue, gene=which.min(LRT_results_rescue_MB$padj), intgroup="Rescue")
LRT_results_rescue_MB_sig$gene <- rownames(LRT_results_rescue_MB_sig)
df_MB <- lapply(LRT_results_rescue_MB_sig$gene, (x) {
  y <- plotCounts(dds_MB_ltr_rescue,  returnData=TRUE, x, intgroup = "Rescue")
  y$feature <- x
  return(y)
})



df_MBs <- do.call(rbind, df_MB)
df_MBs<- df_MBs %>%
  group_by(feature) %>%
  mutate(y.position=max(count))

df_MBs$sample <- rep(sort_samples_data19_compl$shortID,length(df_3))
df_threeruns$number_runs <- rep(sort_samples_data19_compl$number_runs,length(df_3))

padj_MBs <- as.data.frame(LRT_results_rescue_MB_sig[[6]])
padj_MBs$feature <- LRT_results_rescue_MB_sig[[7]]
padj_MBs$`LRT_results_rescue_MB_sig[[6]]` <- lapply(padj_MBs$`LRT_results_rescue_MB_sig[[6]]`, signif, digits=3)

large_padj_MBs <- merge(padj_MBs, df_MBs, by = "feature", all= FALSE)
padj_MBs_small <- large_padj_MBs[!duplicated(large_padj_MBs$feature), ]

plot_a_list <- function(master_list_with_plots, no_of_rows, no_of_cols) {
  
  patchwork::wrap_plots(master_list_with_plots,
                        nrow = no_of_rows, ncol = no_of_cols)
}

# apply seperately for each gene showing the mean
p_MB <- lapply(LRT_results_rescue_MB_sig$gene, function(gene) {
  ggplot(df_MBs[df_MBs[, "feature"]==gene,]) +
    geom_beeswarm(aes(x=Rescue, y=count, color= Rescue),priority = "none")+
    stat_summary(aes(x=Rescue, y=count), fun.y= mean, fun.ymin=mean, fun.ymax=mean, geom="crossbar", width=0.3, linewidth =0.3)+
    #  annotate("text", label= "Padj=", x= 1,y=(padj_MBs_small[match( gene, padj_MBs_small$feature),5]))+
    geom_text(data = as.data.frame(padj_MBs_small),
              aes(label = (padj_MBs_small[match( gene,feature),2]),
                  x = 1.45, y = (padj_MBs_small[match( gene, feature),5]))) +
    ylab("Normalized count")+
    scale_color_manual(values=c("blue3", "red3"))+
    theme_classic()+
    ggtitle(paste(gene))+
    theme(axis.title.y=element_blank(),
          legend.position="none")
})


# finally print your plots
plot_a_list(p_MB, 4, 4)


##### Central Complex+ ####
design_CC <- subset(design_rescue, Tissue =="CC")
CC_counts <- tissue_counts[,colnames(tissue_counts) %in% design_CC$Sample_ID]
dds_CC <- DESeqDataSetFromMatrix(countData = CC_counts, colData = design_CC , design = ~ Colony_ID + Rescue)
# converting counts to integer mode
keep_CC <- rowSums(counts(dds_CC) >= 10) >= 6
dds_CC <- dds_CC[keep_CC,]
# A reduced model without the rescue factor
dds_CC_ltr_rescue <- DESeq(dds_CC, test="LRT", full= ~ Colony_ID + Rescue , reduced = ~ Colony_ID)

LRT_results_rescue_CC <- results(dds_CC_ltr_rescue)
#write.csv(LRT_results_rescue_CC, "DEGs_CC.csv")
LRT_results_rescue_CC_sig <- subset(LRT_results_rescue_CC, padj < 0.05)
LRT_results_rescue_CC_sig 
# DataFrame with 0 rows and 6 columns

##### Antennae ####
design_Ant <- subset(design_rescue, Tissue =="Ant")
Ant_counts <- tissue_counts[,colnames(tissue_counts) %in% design_Ant$Sample_ID]
dds_Ant <- DESeqDataSetFromMatrix(countData = Ant_counts, colData = design_Ant , design = ~ Colony_ID + Rescue)
# converting counts to integer mode
keep_Ant <- rowSums(counts(dds_Ant) >= 10) >= 6
dds_Ant <- dds_Ant[keep_Ant,]
# A reduced model without the rescue factor
dds_Ant_ltr_rescue <- DESeq(dds_Ant, test="LRT", full= ~ Colony_ID + Rescue , reduced = ~ Colony_ID)

LRT_results_rescue_Ant <- results(dds_Ant_ltr_rescue)
#write.csv(LRT_results_rescue_Ant, "DEGs_antennae.csv")
LRT_results_rescue_Ant_sig <- subset(LRT_results_rescue_Ant, padj < 0.05)
LRT_results_rescue_Ant_sig
# 0
LRT_results_rescue_Ant_df<- as.data.frame(LRT_results_rescue_Ant)



#### GOterms upregulated in MB

ontology_type = "CC" # or CC or MF
GO_data_iup <- new("topGOdata", description="GO_Enrichment", ontology=ontology_type, allGenes=gene_list_iup, annot = annFUN.gene2GO, gene2GO=universe)
result_topGO_w01_iup <- runTest(GO_data_iup, algorithm = "weight01", statistic = "fisher")
result_table_w01_iup_CC <- GenTable(GO_data_iup, Fisher = result_topGO_w01_iup, orderBy = "Fisher", ranksOf = "Fisher", topNodes = numSigGenes(GO_data_iup))
result_table_w01_iup_CC

ontology_type = "MF" # or CC or MF
GO_data_iup <- new("topGOdata", description="GO_Enrichment", ontology=ontology_type, allGenes=gene_list_iup, annot = annFUN.gene2GO, gene2GO=universe)
result_topGO_w01_iup <- runTest(GO_data_iup, algorithm = "weight01", statistic = "fisher")
result_table_w01_iup_MF <- GenTable(GO_data_iup, Fisher = result_topGO_w01_iup, orderBy = "Fisher", ranksOf = "Fisher", topNodes = numSigGenes(GO_data_iup))
result_table_w01_iup_MF


####### With the new annotation ####
setwd("~/Documents/Catagl_cognition/From_Eyal")

# Load InterProScan GO mappings 
go_raw <- read.delim("go_raw.tsv", header = FALSE, stringsAsFactors = FALSE)
colnames(go_raw) <- c("transcript_id", "go_terms")
head(go_raw)

# Load gene-transcript mapping file
mapping <- read.csv("Cnig_gn3.1.fasta.transcripts.fasta.trinotate.csv", header = TRUE, stringsAsFactors = FALSE)

# Preview to confirm columns
head(mapping[, 1:2])

# Join GO annotations with gene mapping
go_merged <- merge(go_raw, mapping[, c("X.gene_id", "transcript_id")], by = "transcript_id")
head(go_merged)

library(dplyr)
library(tidyr)

# Clean and collapse
gene_go <- go_merged %>%
filter(go_terms != "") %>%
mutate(go_split = strsplit(go_terms, "[,|]")) %>%
unnest(go_split) %>%
group_by(X.gene_id) %>%
summarise(go_terms = paste(unique(go_split), collapse = ",")) %>%
ungroup()

# Remove any suffixes like (InterPro), (Pfam), etc.
gene_go$go_terms <- gsub("(.*?)", "", gene_go$go_terms)

# Build geneID2GO Object for topGO
geneID2GO <- setNames(strsplit(gene_go$go_terms, ","), gene_go$X.gene_id)
write.table(geneID2GO, file = "geneID2GO_universe.tsv", sep = "t", quote = FALSE, row.names = FALSE)
# generate list of the protein IDs
universe_genelist <- names(geneID2GO)

#### MBs rescuers vs non-rescuers ####
MB_up_genes <- as.data.frame(c("gene_3661",  "gene_3662",  "gene_3426",  "gene_1607",  "gene_2555", "gene_1243" , "gene_930" ,  "gene_5298" , "gene_6089",  "gene_5245",  "gene_13149", "gene_9895",  "gene_9900",  "gene_8499",  "gene_7997"))

colnames(MB_up_genes) <- c("genes") 

# generate a reference list of the protein IDs
#testset_genelist_idown <- readLines(LD_down_genes)
#determine matches between both datasets
gene_list_iup <- factor(as.integer(universe_genelist %in% MB_up_genes$genes))
names(gene_list_iup) <- universe_genelist
#set the ontology type you want to investigate, we first choose "Biological Processes"=BP, '93Cellular Components'94=CC, '93Molecular Function'94=MF
ontology_type = "BP" # or CC or MF
GO_data_iup <- new("topGOdata", description="GO_Enrichment", ontology=ontology_type, allGenes=gene_list_iup, annot = annFUN.gene2GO, gene2GO= geneID2GO)

result_topGO_w01_iup <- runTest(GO_data_iup, algorithm = "weight01", statistic = "fisher")
result_table_w01_iup_BP <- GenTable(GO_data_iup, Fisher = result_topGO_w01_iup, orderBy = "Fisher", ranksOf = "Fisher", topNodes = numSigGenes(GO_data_iup))
result_table_w01_iup_BP
#### Graphics ####
result_table_w01_iup_BP_sig <- result_table_w01_iup_BP[which(result_table_w01_iup_BP$Fisher<0.05),]
setwd("~/Documents/Catag_rescue/")
#write.csv(result_table_w01_idown_BP_sig, file=paste0("topGO_result_weight01_significant_OL_Corrprop_BP_down_cooks.csv")
colors <- colorRampPalette( brewer.pal( 12, "Paired" ) )( length(result_table_w01_iup_BP_sig$Term))
#Create the tag cloud
tagcloud(strmultline(result_table_w01_iup_BP_sig$Term), weights = - log(as.numeric(result_table_w01_iup_BP_sig$Fisher)),col=colors)



#### GOterms upregulated in MB

ontology_type = "CC" # or CC or MF
GO_data_iup <- new("topGOdata", description="GO_Enrichment", ontology=ontology_type, allGenes=gene_list_iup, annot = annFUN.gene2GO, gene2GO=geneID2GO)
result_topGO_w01_iup <- runTest(GO_data_iup, algorithm = "weight01", statistic = "fisher")
result_table_w01_iup_CC <- GenTable(GO_data_iup, Fisher = result_topGO_w01_iup, orderBy = "Fisher", ranksOf = "Fisher", topNodes = numSigGenes(GO_data_iup))
result_table_w01_iup_CC

ontology_type = "MF" # or CC or MF
GO_data_iup <- new("topGOdata", description="GO_Enrichment", ontology=ontology_type, allGenes=gene_list_iup, annot = annFUN.gene2GO, gene2GO=geneID2GO)
result_topGO_w01_iup <- runTest(GO_data_iup, algorithm = "weight01", statistic = "fisher")
result_table_w01_iup_MF <- GenTable(GO_data_iup, Fisher = result_topGO_w01_iup, orderBy = "Fisher", ranksOf = "Fisher", topNodes = numSigGenes(GO_data_iup))
result_table_w01_iup_MF

############ reactome ##########
# Load Reactome mappings 
reactome_raw <- read.delim("reactome_raw.tsv", header = FALSE, stringsAsFactors = FALSE)
colnames(reactome_raw) <- c("transcript_id", "reactome_id")
head(reactome_raw)

# Load gene-transcript mapping file
mapping <- read.csv("Cnig_gn3.1.fasta.transcripts.fasta.trinotate.csv", header = TRUE, stringsAsFactors = FALSE)

# Preview to confirm columns
head(mapping[, 1:2])

# Join GO annotations with gene mapping
reactome_merged <- merge(reactome_raw, mapping[, c("X.gene_id", "transcript_id")], by = "transcript_id")
head(reactome_merged)

# Clean and collapse
gene_reactome <- reactome_merged %>%
filter(reactome_id != "") %>%
mutate(react_split = strsplit(reactome_id, "[,|]")) %>%
unnest(react_split) %>%
group_by(X.gene_id) %>%
summarise(reactome_id = paste(unique(react_split), collapse = ",")) %>%
ungroup()

MB_up_reactome <- merge(MB_up_genes, gene_reactome, by.x = "genes", by.y="X.gene_id", all.x=T, all.y = F)
MB_up_reactome_list <- setNames(strsplit(MB_up_reactome$reactome_id, ","), MB_up_reactome$genes)
write.csv(MB_up_reactome, "MB_upDEGs_reactome_pathways.csv")
library(ReactomeContentService4R)

reactome_ids_gene9900 <- MB_up_reactome_list$gene_9900
# Remove "Reactome:" prefix
reactome_ids_gene9900clean <- sub("^Reactome:", "", reactome_ids_gene9900)

# Get pathway info
pathway_info <- lapply(reactome_ids_gene9900clean, function(id) {
  tryCatch(
           ReactomeContentService4R::getPathwaySummation(id),
           error = function(e) NA
  )
  })
# Combine results
names(pathway_info) <- reactome_ids_gene9900clean
pathway_info


