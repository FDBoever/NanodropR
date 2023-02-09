#A260/280 260nm yo 280nm purity ratio ~1.8 is generally accepted as "pure" for DNA (~2.0 for RNA)
#A260/230 260nm to 280nm purity ratio 1.8-2.2 is generally accepted as "pure" for both DNA and RNA
#
#A low A260/A230 ratio may be the result of:
  # Carbohydrate carryover (often a problem with plants).
  # Residual phenol from nucleic acid extraction.
  # Residual guanidine (often used in column-based kits).
  # Glycogen used for precipitation.


#=======================================================================================#
# -- salt contamination?
# has been suggested to wash the pellete twice with EthOH - isopronol ~ salts?
# better homogenization, load them all up!
# https://www.researchgate.net/post/How_to_improve_my_260_230_ratio_for_DNA_high_purity
#=======================================================================================#
library(dplyr)
library(ggplot2)


read_nanodrop <- function(file){
  df <- read.delim(file,skip=2)
  df$analysis <- stringr::str_split(basename(file),pattern=' ')[[1]][1]
  #df$date <- stringr::str_split(basename(file),pattern=' ')[[1]][2]
  return(df)
}

nanodrop_ratios <- function(df){
  df.summary <- df[,c('Date.and.Time','Sample.ID', 'Username','analysis','X230.0','X260.0','X280.0')]
  df.summary$r260_280 <- df.summary$X260.0/df.summary$X280.0
  df.summary$r260_230 <-  df.summary$X260.0/df.summary$X230.0
  return(df.summary)
}

nanodrop_long <- function(df){
  df.long <- reshape2::melt(df,ids=c('Date.and.Time','Sample.ID', 'Username','analysis'))
  df.long$variable <- as.numeric(gsub('X','',df.long$variable))
  return(df.long)
}


#' Adds custom theme to ggplot
#'
#' @param aspect.ratio sets aspect.ratio in ggplot2::theme()
#'
#' @keywords fdb_style
#' @export
#' @examples
#' ggplot(data=mtcars, aes(cyl,mpg)) +
#' geom_point() +
#' geom_smooth(method='lm')+
#' ggtitle('mtcars') +
#' fdb_style()

fdb_style <- function(aspect.ratio=1) {
  #font <- "Helvetica"
  theme_classic() + ggplot2::theme(
    #Set aspect ratio, this is changable
    aspect.ratio = aspect.ratio,

    #Format title
    plot.title = ggplot2::element_text(hjust = 0.5),

    #Axis format
    axis.title = ggplot2::element_text(size = 11, colour = '#000000'),
    axis.text = ggplot2::element_text(size = 10, colour = '#000000'),

    #Legend format
    legend.justification = c(1, 1),
    legend.key.width = unit(0.25, 'cm'),
    legend.key.height = unit(0.55, 'cm'),
    legend.text = ggplot2::element_text(size = 10),
    legend.title = ggplot2::element_text(size = 11),
    #legend.position='none',

    #Blank background
    #This sets the panel background as blank, removing the standard grey ggplot background colour from the plot
    panel.background = ggplot2::element_blank(),

    #Strip background (#This sets the panel background for facet-wrapped plots to white, removing the standard grey ggplot background colour and sets the title size of the facet-wrap title to font size 22)
    strip.background = ggplot2::element_blank(),
    strip.text = ggplot2::element_text(size  = 11,  hjust = 0)
    # theme(strip.background = element_blank())

    #Grid lines
    # panel.grid.minor = ggplot2::element_blank(),
    # panel.grid.major.y = ggplot2::element_line(color="#cbcbcb"),
    # panel.grid.major.x = ggplot2::element_blank(),
  )
}




#readig nanodrop csv, it has wierd characters, R runs into error
read_nanodrop_csv <- function(file){
  rl <- readLines(file, skipNul=T, warn=F)
  rl <- rl[rl!=""]
  sl <- stringr::str_split(rl,pattern='\t')

  df <- plyr::ldply(sl, rbind)  %>% data.frame()
  df[] <- lapply(df, as.character)

  cl.names <- df %>% slice(1) %>%
    c(., recursive=TRUE) %>%
    unname

  cl.names <- gsub('\xff\xfe','',
                   gsub('_Cor','Cor',
                        gsub(" ","_",
                             gsub("\\(","",
                                  gsub("\\)","",
                                       gsub("\\/","",
                                            gsub("\\%","",
                                                 cl.names)))))))
  df <- df[2:nrow(df),]
  colnames(df) <-cl.names

  df$Nucleic_AcidnguL <- as.numeric(df$Nucleic_AcidnguL)
  df$Corrected_nguL <- as.numeric(df$Corrected_nguL)
  df$Impurity_1_A260 <- as.numeric(df$Impurity_1_A260)
  df$Impurity_1_CV <- as.numeric(df$Impurity_1_CV)
  df$A260A280 <- as.numeric(df$A260A280)
  df$A260A230 <- as.numeric(df$A260A230)
  return(df)

}



#multiread functions

multiread_nanodrop <- function(input,dir){
  out <- c()
  for(i in input){
    tmp <- read_nanodrop(file = paste0(dir, i, " AM_table.tsv"))
    tmp$Sample.ID <- paste0( tmp$Sample.ID ,' ', strsplit(i,split=' ')[[1]][2])
    out <- rbind(out,tmp)
  }
  return(out)
}

multiread_nanodrop_csv <- function(input,dir){
  out <- c()
  for(i in input){
    tmp <- read_nanodrop_csv(file = paste0(dir, i, " AM.csv"))
    tmp$Sample_Name <- paste0( tmp$Sample_Name ,' ', strsplit(i,split=' ')[[1]][2])
    out <- rbind(out,tmp)
  }
  return(out)
}


multiread_nanodrop_meta <- function(input,dir){
  out <- c()
  for(i in input){
    tmp <- read.delim(file = paste0(dir, i, " metadata.txt"))
    tmp$Sample_Name <- paste0( tmp$Sample_Name ,' ', strsplit(i,split=' ')[[1]][2])
    out <- rbind(out,tmp)
  }
  return(out)
}
#===







nanodrop <- read.delim('/Users/sa01fd/DATA/NanodropOne_AZY2121604/dsDNA 5_10_2021 8_12_03 AM_table.tsv',skip=2)
nanodrop2 <- read.delim('/Users/sa01fd/DATA/NanodropOne_AZY2121604/dsDNA 5_11_2021 10_44_03 AM_table.tsv',skip=2)
nanodrop <- read.delim('/Users/sa01fd/DATA/NanodropOne_AZY2121604/dsDNA 5_12_2021 7_49_09 AM_table.tsv',skip=2)


read.delim('/Users/sa01fd/DATA/NanodropOne_AZY2121604/dsDNA 5_10_2021 8_12_03 AM_table.tsv',skip=2)


#=======================================================================================#
#
#nanodrop <- read_nanodrop(file='/Users/sa01fd/DATA/NanodropOne_AZY2121604/dsDNA 5_10_2021 8_12_03 AM_table.tsv')
#nanodrop.csv <- read_nanodrop_csv(file = '/Users/sa01fd/DATA/NanodropOne_AZY2121604/dsDNA 5_10_2021 8_12_03 AM.csv')
#nanodrop.meta <- read.delim(file = '/Users/sa01fd/DATA/NanodropOne_AZY2121604/dsDNA 5_10_2021 8_12_03 metadata.txt')

nanodrop <- read_nanodrop(file='/Users/sa01fd/DATA/NanodropOne_AZY2121604/dsDNA 5_11_2021 10_44_03 AM_table.tsv')
nanodrop <- read_nanodrop(file='/Users/sa01fd/DATA/NanodropOne_AZY2121604/dsDNA 5_12_2021 7_49_09 AM_table.tsv')



dir <- '/Users/sa01fd/DATA/NanodropOne_AZY2121604/'
input <- "dsDNA 5_10_2021 8_12_03"
input <- "dsDNA 5_11_2021 10_44_03"



 #=====

input <- "dsDNA 6_18_2021 4_05_06"
input <- "dsDNA 6_22_2021 4_38_33"
input <- "dsDNA 6_23_2021 3_46_16"
input <- "dsDNA 6_24_2021 3_39_16"



#===================================================================================
#
#
#
#===================================================================================

# Initialisation
#--------------------------------------------------#

#Directory where the Nanodrop files are stored
dir <- '/Users/sa01fd/DATA/NanodropOne_AZY2121604/'

#Vector containing input prefixes, the corresponding files will be combined
inputvector <-   c("dsDNA 6_18_2021 4_05_06",
                  "dsDNA 6_22_2021 4_38_33",
                  "dsDNA 6_23_2021 3_46_16",
                  "dsDNA 6_24_2021 3_39_16",
                  "dsDNA 6_25_2021 4_44_49",
                  "dsDNA 6_29_2021 3_37_27")


# Read in multiple Nanodrop experiments using the multi_read functions
nanodrop <- multiread_nanodrop(input=inputvector,dir=dir)
nanodrop.csv <- multiread_nanodrop_csv(input=inputvector,dir=dir)
nanodrop.meta <- multiread_nanodrop_meta(input=inputvector,dir=dir)

#make 2-way analysis possible by splittng protocol column into HighSalt, and PEG
nanodrop.meta <- nanodrop.meta %>% mutate(PEG=ifelse(grepl('PEG',Protocol),'Y','N'))
nanodrop.meta <- nanodrop.meta %>% mutate(NaCl=ifelse(grepl('HS',Protocol),'Y','N'))

#Generate long formatted data
nanodrop.long <- nanodrop_long(df = nanodrop)

#Calculate ratios using nanodrop_ratios
nanodrop.summary <- nanodrop_ratios(df=nanodrop)

# Generating a clean matrix from the absobrance spectra
nanodrop.matrix <- nanodrop[,3:ncol(nanodrop)-1] #removes unwanted columns
rownames(nanodrop.matrix) <- nanodrop.matrix[,1] #Set rownames to sample id column
nanodrop.matrix <- nanodrop.matrix[,3:ncol(nanodrop.matrix)] #remove the unwanted columns


#remove columns with zero variance
nanodrop.matrix <- nanodrop.matrix[,sapply(nanodrop.matrix, var)!=0]


# normalise Nanodrop matrix
nanodrop.norm <- t(nanodrop.matrix)
nanodrop.norm = apply(nanodrop.norm,2,function(x){x/sum(x)})
nanodrop.norm <- t(nanodrop.norm) %>%
  data.frame() %>%
  tibble::rownames_to_column(var='Sample.ID') %>%
  left_join(nanodrop %>% dplyr::select(Date.and.Time,Sample.ID,Username,analysis))

nanodrop.norm.long <- nanodrop_long(df = nanodrop.norm)


nanodrop %>% mutate_at(X220.0:X350.0, ~(scale(.) %>% as.vector))
nanodrop %>% dplyr::select(across(X220.0:X350.0))

# Idenfification of Sample sets
#--------------------------------------------------#

#CCAP strains with Pigmented gDNA extracts
pigmented.sw <- c('CCAP1312/1',
                  'CCAP1311/1',
                  'CCAP1327/1',
                  'CCAP1319/1',
                  'CCAP1333/1')

#remove those samples that are questionable
pigmented.samples <- nanodrop.meta %>%
  dplyr::filter(Strain %in% pigmented.sw) %>%
  dplyr::select(Sample_Name) %>%
  dplyr::pull()


# samples to remove due to low concentration
conc.cutoff <- 5
low.conc.samples <- nanodrop.csv %>%
  dplyr::filter(Nucleic_AcidnguL<5)%>%
  dplyr::select(Sample_Name) %>%
  dplyr::pull()


#NanodropOne idenfified impurities
nanodrop.csv %>% filter(!is.na(Corrected_nguL))

nanodrop.csv %>%
  dplyr::filter(!is.na(Corrected_nguL)) %>%
  ggplot2::ggplot(aes(x=Nucleic_AcidnguL,y=Corrected_nguL)) +
  ggplot2::geom_point() +
  ggplot2::geom_smooth(method='lm') +
  fdb_style()

impure.samples <- nanodrop.csv %>%
  dplyr::filter(!is.na(Corrected_nguL)) %>%
  dplyr::select(Sample_Name) %>%
  dplyr::pull()


#Stringent sample filtering based on combined values
#remove low-cutoff samples, and pigmented samples
selected.samples <- nanodrop.csv %>%
  dplyr::filter(!Sample_Name %in% low.conc.samples) %>% #filters on low concentration samples (< threshold)
  dplyr::filter(!Sample_Name %in% pigmented.samples) %>% # filters out pigmented gDNA samples
  dplyr::filter(A260A230>0.4) %>%
  dplyr::select(Sample_Name) %>%
  dplyr::pull()


# Export
#--------------------------------------------------#

# export an output table
df.out <- nanodrop.csv %>% left_join(nanodrop.meta, by=c('Sample_Name'='Sample_Name'))
write.table(df.out,file='~/DATA/nanodrop.table.txt',sep='\t',quote = FALSE)







# Basic analysis
#--------------------------------------------------#


#Run PCA on scaled absorbance spectra
pca.res <- prcomp(nanodrop.matrix, scale = TRUE)


#Run PCA on restricted sample set
pca.res <- prcomp(nanodrop.matrix[rownames(nanodrop.matrix) %in% selected.samples,], scale = TRUE)


pca.df <- pca.res$x %>% data.frame() %>%
  tibble::rownames_to_column(var='Sample_Name') %>%
  left_join(nanodrop.meta, by=c('Sample_Name'='Sample_Name'))


pca.df %>%
#  filter(!Strain %in% pigmented.sw) %>%
#  filter(PC1>-10)%>%
  ggplot(aes(x=PC1,y=PC2, color = Protocol)) +
  geom_point() +
  scale_color_manual(values =  c('grey60','grey30','firebrick','skyblue')) +
  fdb_style() #+ geom_text(aes(label=Sample_Name))




nanodrop.long %>% left_join(nanodrop.meta, by=c('Sample.ID'='Sample_Name')) %>%
  #filter(Taxon == "Chlorophyta") %>%
  filter(Taxon == "Rhodophyta") %>%
  #filter(Taxon == "Ochrophyta") %>%

  ggplot2::ggplot()+
  ggplot2::geom_line(aes(x=variable,y=value, color=Protocol)) +
  #ggplot2::geom_area(aes(x=variable,y=value, fill=Description),alpha=0.5) +
  ggplot2::facet_wrap(~Organism+Strain,scales = 'free') +
  #geom_text(data= nanodrop.csv  %>% mutate(Sample.ID = Sample_Name) %>%
  #            left_join(nanodrop.meta, by=c('Sample.ID'='Sample_Name')) %>%
  #            filter(Sample.ID %in% selected.samples),
  #         aes(x=255,y=0.3, label = paste0('260/280=',round(A260A280,2),'\n260/230=',round(A260A230,2),'\n',round(Nucleic_AcidnguL,2),'ng/uL')))+
  geom_vline(xintercept=230,linetype='dashed',color='grey')+
  ggplot2::xlab('Wavelength (nm)')+
  ggplot2::ylab('Absorbance')+
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0),breaks = c(220,240,260,280,300,320,340))+
  ggplot2::theme_classic()+fdb_style(aspect.ratio=0.5)+
  scale_color_manual(values =  c('grey60','grey30','firebrick','skyblue'))
#theme(legend.position = 'none')





nanodrop.norm.long %>% left_join(nanodrop.meta, by=c('Sample.ID'='Sample_Name')) %>%
  filter(Taxon == "Chlorophyta") %>%
  #filter(Taxon == "Rhodophyta") %>%
  #filter(Taxon == "Ochrophyta") %>%
  ggplot2::ggplot()+
  ggplot2::geom_line(aes(x=variable,y=value, color=Protocol)) +
  ggplot2::facet_wrap(~Organism+Strain,scales = 'free') +
  geom_vline(xintercept=230,linetype='dashed',color='grey')+
  ggplot2::xlab('Wavelength (nm)')+
  ggplot2::ylab('Absorbance')+
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0),breaks = c(220,240,260,280,300,320,340))+
  ggplot2::theme_classic()+fdb_style(aspect.ratio=0.5)+
  scale_color_manual(values =  c('grey60','grey30','firebrick','skyblue'))


#=
#normalised spectra

p.chloro <- nanodrop.norm.long %>%
  dplyr::left_join(nanodrop.meta, by=c('Sample.ID'='Sample_Name')) %>%
  dplyr::filter(Taxon == "Chlorophyta") %>%
  ggplot2::ggplot()+
  ggplot2::geom_line(aes(x=variable,y=value, color=Protocol)) +
  ggplot2::facet_wrap(~Organism+Strain,scales = 'free',ncol=8) +
  ggplot2::geom_vline(xintercept=230,linetype='dashed',color='grey')+
  ggplot2::xlab('Wavelength (nm)')+
  ggplot2::ylab('Absorbance')+
  ggplot2::scale_y_continuous(expand = c(0, 0)) +
  ggplot2::scale_x_continuous(expand = c(0, 0),breaks = c(220,240,260,280,300,320,340))+
  ggplot2::theme_classic()+fdb_style(aspect.ratio=0.5)+
  ggplot2::scale_color_manual(values =  c('grey60','grey30','firebrick','deepskyblue3'))

ggsave(plot = p.chloro, file='/Users/sa01fd/DATA/HMW_DNA_Seaweeds/Method_development/Chlorophyta_overlay_normalised_spectra.pdf',width=24,height = 4)

p.rhodo <- nanodrop.norm.long %>%
  dplyr::left_join(nanodrop.meta, by=c('Sample.ID'='Sample_Name')) %>%
  dplyr::filter(Taxon == "Rhodophyta") %>%
  ggplot2::ggplot()+
  ggplot2::geom_line(aes(x=variable,y=value, color=Protocol)) +
  ggplot2::facet_wrap(~Organism+Strain,scales = 'free',ncol=11) +
  ggplot2::geom_vline(xintercept=230,linetype='dashed',color='grey')+
  ggplot2::xlab('Wavelength (nm)')+
  ggplot2::ylab('Absorbance')+
  ggplot2::scale_y_continuous(expand = c(0, 0)) +
  ggplot2::scale_x_continuous(expand = c(0, 0),breaks = c(220,240,260,280,300,320,340))+
  ggplot2::theme_classic()+fdb_style(aspect.ratio=0.5)+
  ggplot2::scale_color_manual(values =  c('grey60','grey30','firebrick','deepskyblue3'))

ggsave(plot = p.rhodo, file='/Users/sa01fd/DATA/HMW_DNA_Seaweeds/Method_development/Rhodophyta_overlay_normalised_spectra.pdf',width=30,height = 4)

p.ochro<- nanodrop.norm.long %>%
  dplyr::left_join(nanodrop.meta, by=c('Sample.ID'='Sample_Name')) %>%
  dplyr::filter(Taxon == "Ochrophyta") %>%
  ggplot2::ggplot()+
  ggplot2::geom_line(aes(x=variable,y=value, color=Protocol)) +
  ggplot2::facet_wrap(~Organism+Strain,scales = 'free',ncol=11) +
  ggplot2::geom_vline(xintercept=230,linetype='dashed',color='grey')+
  ggplot2::xlab('Wavelength (nm)')+
  ggplot2::ylab('Absorbance')+
  ggplot2::scale_y_continuous(expand = c(0, 0)) +
  ggplot2::scale_x_continuous(expand = c(0, 0),breaks = c(220,240,260,280,300,320,340))+
  ggplot2::theme_classic()+fdb_style(aspect.ratio=0.5)+
  ggplot2::scale_color_manual(values =  c('grey60','grey30','firebrick','deepskyblue3'))

ggsave(plot = p.ochro, file='/Users/sa01fd/DATA/HMW_DNA_Seaweeds/Method_development/Ochrophyta_overlay_normalised_spectra.pdf',width=30,height = 4)

p.box.comb <- gridExtra::grid.arrange(p.chloro, p.rhodo, p.ochro,nrow=3)
ggsave(plot = p.box.comb, file='/Users/sa01fd/DATA/HMW_DNA_Seaweeds/Method_development/boxplot_collage_non_filtered_per_treatment.pdf',width=30,height = 9)

#=

p.chloro <- nanodrop.long %>%
  dplyr::left_join(nanodrop.meta, by=c('Sample.ID'='Sample_Name')) %>%
  dplyr::filter(Taxon == "Chlorophyta") %>%
  ggplot2::ggplot()+
  ggplot2::geom_line(aes(x=variable,y=value, color=Protocol)) +
  ggplot2::facet_wrap(~Organism+Strain,scales = 'free',ncol=8) +
  ggplot2::geom_vline(xintercept=230,linetype='dashed',color='grey')+
  ggplot2::xlab('Wavelength (nm)')+
  ggplot2::ylab('Absorbance')+
  ggplot2::scale_y_continuous(expand = c(0, 0)) +
  ggplot2::scale_x_continuous(expand = c(0, 0),breaks = c(220,240,260,280,300,320,340))+
  ggplot2::theme_classic()+fdb_style(aspect.ratio=0.5)+
  ggplot2::scale_color_manual(values =  c('grey60','grey30','firebrick','deepskyblue3'))

ggsave(plot = p.chloro, file='/Users/sa01fd/DATA/HMW_DNA_Seaweeds/Method_development/Chlorophyta_overlay_spectra.pdf',width=24,height = 4)

p.rhodo <- nanodrop.long %>%
  dplyr::left_join(nanodrop.meta, by=c('Sample.ID'='Sample_Name')) %>%
  dplyr::filter(Taxon == "Rhodophyta") %>%
  ggplot2::ggplot()+
  ggplot2::geom_line(aes(x=variable,y=value, color=Protocol)) +
  ggplot2::facet_wrap(~Organism+Strain,scales = 'free',ncol=11) +
  ggplot2::geom_vline(xintercept=230,linetype='dashed',color='grey')+
  ggplot2::xlab('Wavelength (nm)')+
  ggplot2::ylab('Absorbance')+
  ggplot2::scale_y_continuous(expand = c(0, 0)) +
  ggplot2::scale_x_continuous(expand = c(0, 0),breaks = c(220,240,260,280,300,320,340))+
  ggplot2::theme_classic()+fdb_style(aspect.ratio=0.5)+
  ggplot2::scale_color_manual(values =  c('grey60','grey30','firebrick','deepskyblue3'))

ggsave(plot = p.rhodo, file='/Users/sa01fd/DATA/HMW_DNA_Seaweeds/Method_development/Rhodophyta_overlay_spectra.pdf',width=30,height = 4)

p.ochro<- nanodrop.long %>%
  dplyr::left_join(nanodrop.meta, by=c('Sample.ID'='Sample_Name')) %>%
  dplyr::filter(Taxon == "Ochrophyta") %>%
  ggplot2::ggplot()+
  ggplot2::geom_line(aes(x=variable,y=value, color=Protocol)) +
  ggplot2::facet_wrap(~Organism+Strain,scales = 'free',ncol=11) +
  ggplot2::geom_vline(xintercept=230,linetype='dashed',color='grey')+
  ggplot2::xlab('Wavelength (nm)')+
  ggplot2::ylab('Absorbance')+
  ggplot2::scale_y_continuous(expand = c(0, 0)) +
  ggplot2::scale_x_continuous(expand = c(0, 0),breaks = c(220,240,260,280,300,320,340))+
  ggplot2::theme_classic()+fdb_style(aspect.ratio=0.5)+
  ggplot2::scale_color_manual(values =  c('grey60','grey30','firebrick','deepskyblue3'))

ggsave(plot = p.ochro, file='/Users/sa01fd/DATA/HMW_DNA_Seaweeds/Method_development/Ochrophyta_overlay_spectra.pdf',width=30,height = 4)


#combined plot
p.box.comb <- gridExtra::grid.arrange(p.chloro, p.rhodo, p.ochro,nrow=3)
ggsave(plot = p.box.comb, file='/Users/sa01fd/DATA/HMW_DNA_Seaweeds/Method_development/boxplot_collage_non_filtered_per_treatment.pdf',width=30,height = 9)





p.chloro <- nanodrop.long %>% left_join(nanodrop.meta, by=c('Sample.ID'='Sample_Name')) %>%
  filter(Taxon == "Chlorophyta") %>%
  ggplot2::ggplot()+
  ggplot2::geom_line(aes(x=variable,y=value, color=NaCl)) +
  ggplot2::facet_wrap(PEG~Organism+Strain,scales = 'free',ncol=8) +
   geom_vline(xintercept=230,linetype='dashed',color='grey')+
  ggplot2::xlab('Wavelength (nm)')+
  ggplot2::ylab('Absorbance')+
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0),breaks = c(220,240,260,280,300,320,340))+
  ggplot2::theme_classic()+fdb_style(aspect.ratio=0.5)+
  scale_color_manual(values = c('grey60','grey30'))+
  ggtitle('Chlorophyta')

p.chloro

ggsave(plot = p.chloro, file='/Users/sa01fd/DATA/HMW_DNA_Seaweeds/Method_development/Chlorophyta_spectra.pdf',width=24,height = 9)




p.rhodo <- nanodrop.long %>% left_join(nanodrop.meta, by=c('Sample.ID'='Sample_Name')) %>%
  filter(Taxon == "Rhodophyta") %>%
  ggplot2::ggplot()+
  ggplot2::geom_line(aes(x=variable,y=value, color=NaCl)) +
  ggplot2::facet_wrap(PEG~Organism+Strain,scales = 'free',ncol=11) +
  geom_vline(xintercept=230,linetype='dashed',color='grey')+
  ggplot2::xlab('Wavelength (nm)')+
  ggplot2::ylab('Absorbance')+
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0),breaks = c(220,240,260,280,300,320,340))+
  ggplot2::theme_classic()+fdb_style(aspect.ratio=0.5)+
  scale_color_manual(values = c('grey60','grey30'))+
  ggtitle('Rhodophyta')

p.rhodo

ggsave(plot = p.rhodo, file='/Users/sa01fd/DATA/HMW_DNA_Seaweeds/Method_development/Rhodophyta_spectra.pdf',width=24,height = 9)




p.ochro <- nanodrop.long %>% left_join(nanodrop.meta, by=c('Sample.ID'='Sample_Name')) %>%
  filter(Taxon == "Ochrophyta") %>%
  ggplot2::ggplot()+
  ggplot2::geom_line(aes(x=variable,y=value, color=NaCl)) +
  ggplot2::facet_wrap(PEG~Organism+Strain,scales = 'free',ncol=11) +
  geom_vline(xintercept=230,linetype='dashed',color='grey')+
  ggplot2::xlab('Wavelength (nm)')+
  ggplot2::ylab('Absorbance')+
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0),breaks = c(220,240,260,280,300,320,340))+
  ggplot2::theme_classic()+fdb_style(aspect.ratio=0.5)+
  scale_color_manual(values = c('grey60','grey30'))+
  ggtitle('Ochrophyta')

p.ochro

ggsave(plot = p.ochro, file='/Users/sa01fd/DATA/HMW_DNA_Seaweeds/Method_development/Ochrophyta_spectra.pdf',width=24,height = 9)




p.pigmented <- nanodrop.long %>% left_join(nanodrop.meta, by=c('Sample.ID'='Sample_Name')) %>%
  filter(Strain %in% pigmented.sw) %>%
  ggplot2::ggplot()+
  ggplot2::geom_line(aes(x=variable,y=value, color=NaCl)) +
  ggplot2::facet_wrap(PEG~Organism+Strain,scales = 'free',ncol=5) +
  geom_vline(xintercept=230,linetype='dashed',color='grey')+
  ggplot2::xlab('Wavelength (nm)')+
  ggplot2::ylab('Absorbance')+
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0),breaks = c(220,240,260,280,300,320,340))+
  ggplot2::theme_classic()+fdb_style(aspect.ratio=0.5)+
  scale_color_manual(values = c('grey60','grey30'))+
  ggtitle('Ochrophyta')

p.pigmented

ggsave(plot = p.pigmented, file='/Users/sa01fd/DATA/HMW_DNA_Seaweeds/Method_development/Ochrophyta_.pigmented.spectra.pdf',width=12,height = 9)





p.pigmented <- nanodrop.long %>% left_join(nanodrop.meta, by=c('Sample.ID'='Sample_Name')) %>%
  dplyr::filter(Strain %in% pigmented.sw) %>%
  ggplot2::ggplot()+
  ggplot2::geom_line(aes(x=variable,y=value, color=Protocol)) +
  ggplot2::facet_wrap(~Organism+Strain,scales = 'free',ncol=5) +
  geom_vline(xintercept=230,linetype='dashed',color='grey')+
  ggplot2::xlab('Wavelength (nm)')+
  ggplot2::ylab('Absorbance')+
  ggplot2::scale_y_continuous(expand = c(0, 0)) +
  ggplot2::scale_x_continuous(expand = c(0, 0),breaks = c(220,240,260,280,300,320,340))+
  ggplot2::theme_classic()+fdb_style(aspect.ratio=0.5)+
  ggplot2::scale_color_manual(values=c('grey70','grey50','grey35','black'))+
  ggplot2::theme(legend.position = 'none')

ggsave(plot = p.pigmented, file='/Users/sa01fd/DATA/HMW_DNA_Seaweeds/Method_development/Ochrophyta_.pigmented_overlay.spectra.pdf',width=12,height = 4)



p.pigmented <- nanodrop.long %>% left_join(nanodrop.meta, by=c('Sample.ID'='Sample_Name')) %>%
  dplyr::filter(Strain %in% pigmented.sw) %>%
  dplyr::filter(Protocol == 'HS-CTAB-PEG') %>%
  ggplot2::ggplot()+
  ggplot2::geom_line(aes(x=variable,y=value, color=Protocol)) +
  ggplot2::facet_wrap(~Organism+Strain,scales = 'free',ncol=5) +
  geom_vline(xintercept=230,linetype='dashed',color='grey')+
  ggplot2::xlab('Wavelength (nm)')+
  ggplot2::ylab('Absorbance')+
  ggplot2::scale_y_continuous(expand = c(0, 0)) +
  ggplot2::scale_x_continuous(expand = c(0, 0),breaks = c(220,240,260,280,300,320,340))+
  ggplot2::theme_classic()+fdb_style(aspect.ratio=0.5)+
  ggplot2::scale_color_manual(values=c('black'))+
  ggplot2::theme(legend.position = 'none')

ggsave(plot = p.pigmented, file='/Users/sa01fd/DATA/HMW_DNA_Seaweeds/Method_development/Ochrophyta_.pigmented_HS-CTAB-PEG.spectra.pdf',width=12,height = 4)



nanodrop.long %>% left_join(nanodrop.meta, by=c('Sample.ID'='Sample_Name')) %>%
  dplyr::filter(Strain %in% pigmented.sw) %>%
  ggplot2::ggplot()+
  ggplot2::geom_line(aes(x=variable,y=log(value), color=Protocol)) +
  ggplot2::facet_wrap(~Organism+Strain,scales = 'free',ncol=5) +
  geom_vline(xintercept=230,linetype='dashed',color='grey')+
  ggplot2::xlab('Wavelength (nm)')+
  ggplot2::ylab('Absorbance')+
  ggplot2::scale_y_continuous(expand = c(0, 0)) +
  ggplot2::scale_x_continuous(expand = c(0, 0),breaks = c(220,240,260,280,300,320,340))+
  ggplot2::theme_classic()+fdb_style(aspect.ratio=0.5)+
  ggplot2::scale_color_manual(values=c('grey70','grey50','grey35','black'))+
  ggplot2::theme(legend.position = 'none')


p.pigmented.decl <- nanodrop.csv %>% left_join(nanodrop.meta, by=c('Sample_Name'='Sample_Name')) %>%
  dplyr::filter(Strain %in% pigmented.sw) %>%
  ggplot2::ggplot(aes(x=Protocol,y=Nucleic_AcidnguL, color=Protocol))+
  ggplot2::geom_point() +
  ggplot2::geom_line(aes(group=Strain))+
  geom_vline(xintercept=230,linetype='dashed',color='grey')+
  ggplot2::ylab('Absorbance')+
  ggplot2::theme_classic()+fdb_style(aspect.ratio=2)+
  ggplot2::scale_color_manual(values=c('grey70','grey50','grey35','black'))+
  ggplot2::theme(legend.position = 'none')+
  ggplot2::theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

ggsave(plot = p.pigmented.decl, file='/Users/sa01fd/DATA/HMW_DNA_Seaweeds/Method_development/Ochrophyta_.pigmented_drop_in_abs.pdf',width=4,height = 4)



p_scatter <- df.nan %>%
  dplyr::filter(Sample_Name %in% selected.samples) %>%
  ggplot2::ggplot(aes(x=A260A280,y=A260A230, color=Protocol)) +
  ggplot2::geom_point() +
  fdb_style() +
  scale_color_manual(values =  c('grey60','grey30','firebrick','skyblue'))


ggsave(plot = p_scatter, file='/Users/sa01fd/DATA/HMW_DNA_Seaweeds/Method_development/scatter_A280_A230_filtered.pdf',width=4,height = 4)

p_scatter <- df.nan %>%
  dplyr::filter(Sample_Name %in% selected.samples) %>%
  ggplot2::ggplot(aes(x=Nucleic_AcidnguL,y=A260A230, color=Protocol)) +
  ggplot2::geom_point() +
  fdb_style() +
  scale_color_manual(values =  c('grey60','grey30','firebrick','skyblue'))


ggsave(plot = p_scatter, file='/Users/sa01fd/DATA/HMW_DNA_Seaweeds/Method_development/scatter_conc_A230_filtered.pdf',width=4,height = 4)




p_scatter <- df.nan %>% tidyr::pivot_longer(cols=c(A260A230,A260A280,Nucleic_AcidnguL),names_to='variable', values_to='values') %>%
  dplyr::filter(Sample_Name %in% selected.samples) %>%
  ggplot2::ggplot(aes(x=values, fill=Protocol)) +
  ggplot2::geom_density() +
  fdb_style(aspect.ratio = 0.5) +
  scale_fill_manual(values =  c('grey60','grey30','firebrick','skyblue')) +
  #facet_wrap(Protocol~variable,ncol=3,scales='free_y')
  facet_grid(Protocol~variable,scales='free')




p <- nanodrop.long %>% left_join(nanodrop.meta, by=c('Sample.ID'='Sample_Name')) %>%
  #filter(Sample.ID %in% selected.samples) %>%

  ggplot2::ggplot()+
  ggplot2::geom_line(aes(x=variable,y=value, color=Protocol)) +
  #ggplot2::geom_area(aes(x=variable,y=value, fill=Protocol),alpha=0.5) +
  ggplot2::facet_wrap(Protocol~Organism+Strain,scales = 'free',ncol = 6) +
  #geom_text(data= nanodrop.csv  %>% mutate(Sample.ID = Sample_Name) %>%
  #            left_join(nanodrop.meta, by=c('Sample.ID'='Sample_Name')) %>%
  #            filter(Sample.ID %in% selected.samples),
  #          aes(x=255,y=0.3, label = paste0('260/280=',round(A260A280,2),'\n260/230=',round(A260A230,2),'\n',round(Nucleic_AcidnguL,2),'ng/uL')))+
  geom_vline(xintercept=230,linetype='dashed',color='grey')+
  ggplot2::xlab('Wavelength (nm)')+
  ggplot2::ylab('Absorbance')+
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0),breaks = c(220,240,260,280,300,320,340))+
  ggplot2::theme_classic()+fdb_style(aspect.ratio=0.5)+
  scale_fill_manual(values = colorRampPalette(RColorBrewer::brewer.pal(8, "Set1"))(length(selected.samples)))+
  scale_color_manual(values = colorRampPalette(RColorBrewer::brewer.pal(8, "Set1"))(length(selected.samples)))+
  theme(legend.position = 'none')


p


s






df.nan <- nanodrop.csv %>%
  dplyr::left_join(nanodrop.meta, by=c('Sample_Name'='Sample_Name'))

df.nan <- nanodrop.csv %>%
  dplyr::left_join(nanodrop.meta, by=c('Sample_Name'='Sample_Name')) %>%
  dplyr::filter(Sample_Name %in% selected.samples)


res.aov2 <- aov(A260A230 ~ Protocol, data = df.nan)
summary(res.aov2)


res.aov2 <- aov(A260A230 ~ NaCl + PEG, data = df.nan)
summary(res.aov2)

#with interaction term
res.aov3 <- aov(A260A230 ~ NaCl + PEG + NaCl:PEG, data = df.nan)
summary(res.aov3)


#Check homogeneity of variance
car::leveneTest(A260A230 ~ NaCl*PEG, data = df.nan)

#>0.05 is considered homogenious variance

#Check normality
plot(res.aov3, 2)

# Extract the residuals
aov_residuals <- residuals(object = res.aov3)
# Run Shapiro-Wilk test
shapiro.test(x = aov_residuals )

#Shapiro-Wilk test on the ANOVA residuals (W = 0.98, p = 0.11) which finds no indication that normality is violated.
#http://www.sthda.com/english/wiki/two-way-anova-test-in-r
#Compute two-way ANOVA test in R for unbalanced designs
library(car)

car::leveneTest(A260A230 ~ NaCl * PEG * Taxon, data = df.nan)


res.aov3 <- aov(A260A230 ~ NaCl + PEG + NaCl:PEG, data = df.nan)

res.aov4 <- aov(A260A230 ~ NaCl + PEG + Taxon, data = df.nan)
aov_residuals <- residuals(object = res.aov4)
# Run Shapiro-Wilk test
shapiro.test(x = aov_residuals )


car::Anova(res.aov4, type = "III")
car::Anova(res.aov4, type = "I")
car::Anova(res.aov4, type = "II")


my_anova <- aov(A260A280 ~ NaCl + PEG + Taxon, data = df.nan)
car::Anova(my_anova, type = "III")

my_anova <- aov(A260A280 ~ NaCl + PEG + Taxon, data = df.nan)
car::Anova(my_anova, type = "III")





plot(res.aov3, 1)

#========================
#Summary Tables

#non-filtered
df.summary.non.filtered <- df.nan %>%
  dplyr::group_by(NaCl, PEG,Taxon) %>%
  dplyr::summarise(
    count = n(),
    mean_A260A230 = mean(A260A230, na.rm = TRUE),
    sd_A260A230 = sd(A260A230, na.rm = TRUE),
    mean_A260A280 = mean(A260A280, na.rm = TRUE),
    sd_A260A280 = sd(A260A280, na.rm = TRUE),
    mean_Nucleic_AcidnguL = mean(Nucleic_AcidnguL, na.rm = TRUE),
    sd_Nucleic_AcidnguL = sd(Nucleic_AcidnguL, na.rm = TRUE),
  )

df.summary.non.filtered

#filtered
df.summary.filtered <- df.nan %>%
  dplyr::filter(Sample_Name %in% selected.samples) %>%
  dplyr::group_by(NaCl, PEG,Taxon) %>%
  dplyr::summarise(
    count = n(),
    mean_A260A230 = mean(A260A230, na.rm = TRUE),
    sd_A260A230 = sd(A260A230, na.rm = TRUE),
    mean_A260A280 = mean(A260A280, na.rm = TRUE),
    sd_A260A280 = sd(A260A280, na.rm = TRUE),
    mean_Nucleic_AcidnguL = mean(Nucleic_AcidnguL, na.rm = TRUE),
    sd_Nucleic_AcidnguL = sd(Nucleic_AcidnguL, na.rm = TRUE),
  )

df.summary.filtered

#====================




df.nan %>% group_by(NaCl, PEG) %>%
  summarise(
    count = n(),
    mean = mean(A260A230, na.rm = TRUE),
    sd = sd(A260A230, na.rm = TRUE)
  )




#Bar Chart if needed
p.barchrt <- df.nan %>% group_by(NaCl, PEG,Taxon) %>%
  summarise(
    count = n(),
    mean = mean(A260A230, na.rm = TRUE),
    sd = sd(A260A230, na.rm = TRUE)
  ) %>% ggplot(aes(x=Taxon,
           y= mean,
           ymin=mean-sd,
           ymax=mean+sd,
           fill=PEG)) +
  ggplot2::geom_bar(position="dodge", stat = "identity") +
  ggplot2::geom_errorbar( position = position_dodge(0.9), colour="black", width = 0.3) +
  ggplot2::geom_point(position=position_dodge(.9), aes(y=mean), colour='black')+facet_wrap(~NaCl) +
  ggplot2::scale_fill_manual(values = c('grey60','grey30'))+
  ggplot2::scale_color_manual(values = c('grey60','grey30'))+
  ggplot2::scale_y_continuous(expand = c(0, 0))+
    ggplot2::theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  fdb_style(aspect.ratio=2)


p.barchrt
ggsave(plot = p.barchrt, file='/Users/sa01fd/DATA/HMW_DNA_Seaweeds/Method_development/barchart_230_non_filtered.pdf',width=6,height = 4)


p.barchrt <- df.nan %>%
  dplyr::filter(Sample_Name %in% selected.samples) %>%
  dplyr::group_by(NaCl, PEG,Taxon) %>%
  summarise(
    count = n(),
    mean = mean(A260A230, na.rm = TRUE),
    sd = sd(A260A230, na.rm = TRUE)
  ) %>% ggplot(aes(x=Taxon,
                   y= mean,
                   ymin=mean-sd,
                   ymax=mean+sd,
                   fill=PEG)) +
  ggplot2::geom_bar(position="dodge", stat = "identity") +
  ggplot2::geom_errorbar( position = position_dodge(0.9), colour="black", width = 0.3) +
  ggplot2::geom_point(position=position_dodge(.9), aes(y=mean), colour='black')+facet_wrap(~NaCl) +
  ggplot2::scale_fill_manual(values = c('grey60','grey30'))+
  ggplot2::scale_color_manual(values = c('grey60','grey30'))+
  ggplot2::scale_y_continuous(expand = c(0, 0))+
  fdb_style(aspect.ratio=2) +
  ggplot2::theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

ggsave(plot = p.barchrt, file='/Users/sa01fd/DATA/HMW_DNA_Seaweeds/Method_development/barchart_230_filtered.pdf',width=6,height = 4)



p.barchrt <- df.nan %>%
  dplyr::filter(Sample_Name %in% selected.samples) %>%
  dplyr::group_by(NaCl, PEG,Taxon) %>%
  summarise(
    count = n(),
    mean = mean(A260A280, na.rm = TRUE),
    sd = sd(A260A280, na.rm = TRUE)
  ) %>% ggplot(aes(x=Taxon,
                   y= mean,
                   ymin=mean-sd,
                   ymax=mean+sd,
                   fill=PEG)) +
  #ggplot2::geom_bar(position="dodge", stat = "identity") +
  ggplot2::geom_errorbar( position = position_dodge(0.9), colour="black", width = 0.3) +
  ggplot2::geom_point(position=position_dodge(.9), aes(y=mean), colour='black')+facet_wrap(~NaCl) +
  ggplot2::scale_fill_manual(values = c('grey60','grey30'))+
  ggplot2::scale_color_manual(values = c('grey60','grey30'))+
  ggplot2::scale_y_continuous(expand = c(0, 0))+
  fdb_style(aspect.ratio=2) +
  ggplot2::theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

ggsave(plot = p.barchrt, file='/Users/sa01fd/DATA/HMW_DNA_Seaweeds/Method_development/barchart_280_filtered.pdf',width=6,height = 4)









model.tables(res.aov3, type="means", se = TRUE)

#==========================================================
# BOXPLOTS FOR OVERALL DATASET PER TREATMENT

#A260/A230
p.boxplt.A260A230 <- nanodrop.csv %>%
  dplyr::left_join(nanodrop.meta, by=c('Sample_Name'='Sample_Name')) %>%
  ggplot2::ggplot(aes(x=Protocol,y=A260A230))+
  stat_boxplot( aes(x=Protocol,y=A260A230),
                geom='errorbar', linetype=1, width=0.5,color='grey30')+
  ggplot2::geom_boxplot(aes(fill=Protocol),outlier.shape=NA)+
  ggplot2::geom_jitter(aes(x=Protocol,y=A260A230),alpha=.2,width = 0.15) +
  ggplot2::xlab('')+
  ggplot2::ylab('A260/A230')+
  ggplot2::theme_classic()+fdb_style(aspect.ratio=2)+
  ggplot2::theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  ggplot2::scale_fill_manual(values=c('white','grey90','grey70','grey50'))+
  ggplot2::theme(legend.position = 'none')
p.boxplt.A260A230

ggsave(plot = p.boxplt.A260A230, file='/Users/sa01fd/DATA/HMW_DNA_Seaweeds/Method_development/boxplot_A260A230_non_filtered_per_treatment.pdf',width=6,height = 4)

#A260/A280
p.boxplt.A260A280 <- nanodrop.csv %>%
  dplyr::left_join(nanodrop.meta, by=c('Sample_Name'='Sample_Name')) %>%
  ggplot2::ggplot(aes(x=Protocol,y=A260A280))+
  stat_boxplot( aes(x=Protocol,y=A260A280),
                geom='errorbar', linetype=1, width=0.5,color='grey30')+
  ggplot2::geom_boxplot(aes(fill=Protocol),outlier.shape=NA)+
  ggplot2::geom_jitter(aes(x=Protocol,y=A260A280),alpha=.2,width = 0.15) +
  ggplot2::xlab('')+
  ggplot2::ylab('A260/A280')+
  ggplot2::theme_classic()+fdb_style(aspect.ratio=2)+
  ggplot2::theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  ggplot2::scale_fill_manual(values=c('white','grey90','grey70','grey50'))+
  ggplot2::theme(legend.position = 'none')
p.boxplt.A260A280

ggsave(plot = p.boxplt.A260A280, file='/Users/sa01fd/DATA/HMW_DNA_Seaweeds/Method_development/boxplot_A260A280_non_filtered_per_treatment.pdf',width=6,height = 4)

#Nucleic_AcidnguL
p.boxplt.Nucleic_AcidnguL <- nanodrop.csv %>%
  dplyr::left_join(nanodrop.meta, by=c('Sample_Name'='Sample_Name')) %>%
  ggplot2::ggplot(aes(x=Protocol,y=Nucleic_AcidnguL))+
  stat_boxplot( aes(x=Protocol,y=Nucleic_AcidnguL),
                geom='errorbar', linetype=1, width=0.5,color='grey30')+
  ggplot2::geom_boxplot(aes(fill=Protocol),outlier.shape=NA)+
  ggplot2::geom_jitter(aes(x=Protocol,y=Nucleic_AcidnguL),alpha=.2,width = 0.15) +
  ggplot2::xlab('')+
  ggplot2::ylab('Nucleic Acid (ng/uL)')+
  ggplot2::theme_classic()+fdb_style(aspect.ratio=2)+
  ggplot2::theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  ggplot2::scale_fill_manual(values=c('white','grey90','grey70','grey50'))+
  ggplot2::theme(legend.position = 'none')
p.boxplt.Nucleic_AcidnguL

ggsave(plot = p.boxplt.Nucleic_AcidnguL, file='/Users/sa01fd/DATA/HMW_DNA_Seaweeds/Method_development/boxplot_Nucleic_AcidnguL_non_filtered_per_treatment.pdf',width=6,height = 4)

#combined plot
p.box.comb <- gridExtra::grid.arrange(p.boxplt.A260A230, p.boxplt.A260A280, p.boxplt.Nucleic_AcidnguL,nrow=1)
ggsave(plot = p.box.comb, file='/Users/sa01fd/DATA/HMW_DNA_Seaweeds/Method_development/boxplot_collage_non_filtered_per_treatment.pdf',width=6,height = 4)


#FILTERED
#A260/A230
p.boxplt.A260A230 <- nanodrop.csv %>%
  dplyr::filter(Sample_Name %in% selected.samples) %>%
  dplyr::left_join(nanodrop.meta, by=c('Sample_Name'='Sample_Name')) %>%
  ggplot2::ggplot(aes(x=Protocol,y=A260A230))+
  stat_boxplot( aes(x=Protocol,y=A260A230),
                geom='errorbar', linetype=1, width=0.5,color='grey30')+
  ggplot2::geom_boxplot(aes(fill=Protocol),outlier.shape=NA)+
  ggplot2::geom_jitter(aes(x=Protocol,y=A260A230),alpha=.2,width = 0.15) +
  ggplot2::xlab('')+
  ggplot2::ylab('A260/A230')+
  ggplot2::theme_classic()+fdb_style(aspect.ratio=2)+
  ggplot2::theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  ggplot2::scale_fill_manual(values=c('white','grey90','grey70','grey50'))+
  ggplot2::theme(legend.position = 'none')
p.boxplt.A260A230

ggsave(plot = p.boxplt.A260A230, file='/Users/sa01fd/DATA/HMW_DNA_Seaweeds/Method_development/boxplot_A260A230_filtered_per_treatment.pdf',width=6,height = 4)


#A260/A280
p.boxplt.A260A280 <- nanodrop.csv %>%
  dplyr::filter(Sample_Name %in% selected.samples) %>%
  dplyr::left_join(nanodrop.meta, by=c('Sample_Name'='Sample_Name')) %>%
  ggplot2::ggplot(aes(x=Protocol,y=A260A280))+
  stat_boxplot( aes(x=Protocol,y=A260A280),
                geom='errorbar', linetype=1, width=0.5,color='grey30')+
  ggplot2::geom_boxplot(aes(fill=Protocol),outlier.shape=NA)+
  ggplot2::geom_jitter(aes(x=Protocol,y=A260A280),alpha=.2,width = 0.15) +
  ggplot2::xlab('')+
  ggplot2::ylab('A260/A280')+
  ggplot2::theme_classic()+fdb_style(aspect.ratio=2)+
  ggplot2::theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  ggplot2::scale_fill_manual(values=c('white','grey90','grey70','grey50'))+
  ggplot2::theme(legend.position = 'none')
p.boxplt.A260A280

ggsave(plot = p.boxplt.A260A280, file='/Users/sa01fd/DATA/HMW_DNA_Seaweeds/Method_development/boxplot_A260A280_filtered_per_treatment.pdf',width=6,height = 4)

#Nucleic_AcidnguL
p.boxplt.Nucleic_AcidnguL <- nanodrop.csv %>%
  dplyr::filter(Sample_Name %in% selected.samples) %>%
  dplyr::left_join(nanodrop.meta, by=c('Sample_Name'='Sample_Name')) %>%
  ggplot2::ggplot(aes(x=Protocol,y=Nucleic_AcidnguL))+
  stat_boxplot( aes(x=Protocol,y=Nucleic_AcidnguL),
                geom='errorbar', linetype=1, width=0.5,color='grey30')+
  ggplot2::geom_boxplot(aes(fill=Protocol),outlier.shape=NA)+
  ggplot2::geom_jitter(aes(x=Protocol,y=Nucleic_AcidnguL),alpha=.2,width = 0.15) +
  ggplot2::xlab('')+
  ggplot2::ylab('Nucleic Acid (ng/uL)')+
  ggplot2::theme_classic()+fdb_style(aspect.ratio=2)+
  ggplot2::theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  ggplot2::scale_fill_manual(values=c('white','grey90','grey70','grey50'))+
  ggplot2::theme(legend.position = 'none')
p.boxplt.Nucleic_AcidnguL

ggsave(plot = p.boxplt.Nucleic_AcidnguL, file='/Users/sa01fd/DATA/HMW_DNA_Seaweeds/Method_development/boxplot_Nucleic_AcidnguL_filtered_per_treatment.pdf',width=6,height = 4)


#combined plot
p.box.comb <- gridExtra::grid.arrange(p.boxplt.A260A230, p.boxplt.A260A280, p.boxplt.Nucleic_AcidnguL,nrow=1)
ggsave(plot = p.box.comb, file='/Users/sa01fd/DATA/HMW_DNA_Seaweeds/Method_development/boxplot_collage_filtered_per_treatment.pdf',width=6,height = 4)

#====







#without the pigmented ones?
p.boxplt <- nanodrop.csv %>%
  dplyr::left_join(nanodrop.meta, by=c('Sample_Name'='Sample_Name')) %>%
  dplyr::filter(!Strain %in% pigmented.sw) %>%
  ggplot2::ggplot(aes(x=Protocol,y=A260A230))+
  stat_boxplot( aes(x=Protocol,y=A260A230),
                geom='errorbar', linetype=1, width=0.5,color='grey30')+
  ggplot2::geom_boxplot(aes(fill=Protocol),outlier.shape=NA)+
  ggplot2::geom_jitter(aes(x=Protocol,y=A260A230),alpha=.2,width = 0.15) +
  ggplot2::xlab('Protocol')+
  ggplot2::ylab('A260/A230')+
  ggplot2::theme_classic()+fdb_style(aspect.ratio=1)+
  ggplot2::theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  ggplot2::scale_fill_manual(values=c('white','grey90','grey70','grey50'))
p.boxplt








#A260/A230
p.boxplt.A260A230 <- nanodrop.csv %>%
  dplyr::left_join(nanodrop.meta, by=c('Sample_Name'='Sample_Name')) %>%
  ggplot2::ggplot(aes(x=Protocol,y=A260A230))+
  stat_boxplot( aes(x=Protocol,y=A260A230),
                geom='errorbar', linetype=1, width=0.5,color='grey30')+
  ggplot2::geom_boxplot(aes(fill=Protocol),outlier.shape=NA)+
  ggplot2::geom_jitter(aes(x=Protocol,y=A260A230),alpha=.3,width = 0.1) +
  ggplot2::xlab('')+
  ggplot2::ylab('A260/A230')+
  ggplot2::theme_classic()+fdb_style(aspect.ratio=2)+
  ggplot2::facet_wrap(~Taxon)+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  scale_fill_manual(values=c('white','grey90','grey70','grey50'))+
  theme(legend.position='none')
p.boxplt.A260A230

ggsave(plot = p.boxplt.A260A230, file='/Users/sa01fd/DATA/HMW_DNA_Seaweeds/Method_development/boxplot_A260A230_non_filtered_per_Taxon.pdf',width=6,height = 4)

#A260/A280
p.boxplt.A260A280 <- nanodrop.csv %>%
  dplyr::left_join(nanodrop.meta, by=c('Sample_Name'='Sample_Name')) %>%
  ggplot2::ggplot(aes(x=Protocol,y=A260A280))+
  stat_boxplot( aes(x=Protocol,y=A260A280),
                geom='errorbar', linetype=1, width=0.5,color='grey30')+
  ggplot2::geom_boxplot(aes(fill=Protocol),outlier.shape=NA)+
  ggplot2::geom_jitter(aes(x=Protocol,y=A260A280),alpha=.3,width = 0.1) +
  ggplot2::xlab('')+
  ggplot2::ylab('A260/A280')+
  ggplot2::theme_classic()+fdb_style(aspect.ratio=2)+
  ggplot2::facet_wrap(~Taxon)+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  scale_fill_manual(values=c('white','grey90','grey70','grey50'))+
  theme(legend.position='none')
p.boxplt.A260A280

ggsave(plot = p.boxplt.A260A280, file='/Users/sa01fd/DATA/HMW_DNA_Seaweeds/Method_development/boxplot_A260A280_non_filtered_per_Taxon.pdf',width=6,height = 4)

#A260/A280
p.boxplt.Nucleic_AcidnguL <- nanodrop.csv %>%
  dplyr::left_join(nanodrop.meta, by=c('Sample_Name'='Sample_Name')) %>%
  ggplot2::ggplot(aes(x=Protocol,y=Nucleic_AcidnguL))+
  stat_boxplot( aes(x=Protocol,y=Nucleic_AcidnguL),
                geom='errorbar', linetype=1, width=0.5,color='grey30')+
  ggplot2::geom_boxplot(aes(fill=Protocol),outlier.shape=NA)+
  ggplot2::geom_jitter(aes(x=Protocol,y=Nucleic_AcidnguL),alpha=.3,width = 0.1) +
  ggplot2::xlab('')+
  ggplot2::ylab('Nucleic Acid (ng/uL)')+
  ggplot2::theme_classic()+fdb_style(aspect.ratio=2)+
  ggplot2::facet_wrap(~Taxon)+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  scale_fill_manual(values=c('white','grey90','grey70','grey50'))+
  theme(legend.position='none')
p.boxplt.Nucleic_AcidnguL

ggsave(plot = p.boxplt.Nucleic_AcidnguL, file='/Users/sa01fd/DATA/HMW_DNA_Seaweeds/Method_development/boxplot_Nucleic_AcidnguL_non_filtered_per_Taxon.pdf',width=6,height = 4)

p.box.comb <- gridExtra::grid.arrange(p.boxplt.A260A230, p.boxplt.A260A280, p.boxplt.Nucleic_AcidnguL,nrow=3)
ggsave(plot = p.box.comb, file='/Users/sa01fd/DATA/HMW_DNA_Seaweeds/Method_development/boxplot_collage_non_filtered_per_treatment.pdf',width=6,height = 4)

p.box.comb <- gridExtra::grid.arrange(p.boxplt.A260A230, p.boxplt.A260A280, p.boxplt.Nucleic_AcidnguL,nrow=1)
ggsave(plot = p.box.comb, file='/Users/sa01fd/DATA/HMW_DNA_Seaweeds/Method_development/boxplot_collage_h_non_filtered_per_treatment.pdf',width=12,height = 4)



#filtered
p.boxplt.A260A230 <- nanodrop.csv %>%
  dplyr::left_join(nanodrop.meta, by=c('Sample_Name'='Sample_Name')) %>%
  dplyr::filter(!Strain %in% selected.samples) %>%
  ggplot2::ggplot(aes(x=Protocol,y=A260A230))+
  stat_boxplot( aes(x=Protocol,y=A260A230),
                geom='errorbar', linetype=1, width=0.5,color='grey30')+
  ggplot2::geom_boxplot(aes(fill=Protocol),outlier.shape=NA)+
  ggplot2::geom_jitter(aes(x=Protocol,y=A260A230),alpha=.3,width = 0.1) +
  ggplot2::xlab('')+
  ggplot2::ylab('A260/A230')+
  ggplot2::theme_classic()+fdb_style(aspect.ratio=2)+
  ggplot2::facet_wrap(~Taxon)+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  scale_fill_manual(values=c('white','grey90','grey70','grey50'))+
  theme(legend.position='none')
p.boxplt.A260A230

ggsave(plot = p.boxplt.A260A230, file='/Users/sa01fd/DATA/HMW_DNA_Seaweeds/Method_development/boxplot_A260A230_filtered_per_Taxon.pdf',width=6,height = 4)

#A260/A280 #EXTRA FILTER AT <2.2
p.boxplt.A260A280 <- nanodrop.csv %>%
  dplyr::left_join(nanodrop.meta, by=c('Sample_Name'='Sample_Name')) %>%
  dplyr::filter(!Strain %in% selected.samples) %>%
  dplyr::filter(A260A280<2.2) %>%
  ggplot2::ggplot(aes(x=Protocol,y=A260A280))+
  stat_boxplot( aes(x=Protocol,y=A260A280),
                geom='errorbar', linetype=1, width=0.5,color='grey30')+
  ggplot2::geom_boxplot(aes(fill=Protocol),outlier.shape=NA)+
  ggplot2::geom_jitter(aes(x=Protocol,y=A260A280),alpha=.3,width = 0.1) +
  ggplot2::xlab('')+
  ggplot2::ylab('A260/A280')+
  ggplot2::theme_classic()+fdb_style(aspect.ratio=2)+
  ggplot2::facet_wrap(~Taxon)+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  scale_fill_manual(values=c('white','grey90','grey70','grey50'))+
  theme(legend.position='none')
p.boxplt.A260A280

ggsave(plot = p.boxplt.A260A280, file='/Users/sa01fd/DATA/HMW_DNA_Seaweeds/Method_development/boxplot_A260A280_filtered_per_Taxon.pdf',width=6,height = 4)

#A260/A280   #EXTRA FILTER AT <500
p.boxplt.Nucleic_AcidnguL <- nanodrop.csv %>%
  dplyr::left_join(nanodrop.meta, by=c('Sample_Name'='Sample_Name')) %>%
  dplyr::filter(!Strain %in% selected.samples) %>%
  dplyr::filter(Nucleic_AcidnguL<500) %>%

  ggplot2::ggplot(aes(x=Protocol,y=Nucleic_AcidnguL))+
  stat_boxplot( aes(x=Protocol,y=Nucleic_AcidnguL),
                geom='errorbar', linetype=1, width=0.5,color='grey30')+
  ggplot2::geom_boxplot(aes(fill=Protocol),outlier.shape=NA)+
  ggplot2::geom_jitter(aes(x=Protocol,y=Nucleic_AcidnguL),alpha=.3,width = 0.1) +
  ggplot2::xlab('')+
  ggplot2::ylab('Nucleic Acid (ng/uL)')+
  ggplot2::theme_classic()+fdb_style(aspect.ratio=2)+
  ggplot2::facet_wrap(~Taxon)+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  scale_fill_manual(values=c('white','grey90','grey70','grey50'))+
  theme(legend.position='none')
p.boxplt.Nucleic_AcidnguL

ggsave(plot = p.boxplt.Nucleic_AcidnguL, file='/Users/sa01fd/DATA/HMW_DNA_Seaweeds/Method_development/boxplot_Nucleic_AcidnguL_non_filtered_per_Taxon.pdf',width=6,height = 4)

p.box.comb <- gridExtra::grid.arrange(p.boxplt.A260A230, p.boxplt.A260A280, p.boxplt.Nucleic_AcidnguL,nrow=3)
ggsave(plot = p.box.comb, file='/Users/sa01fd/DATA/HMW_DNA_Seaweeds/Method_development/boxplot_collage_filtered_per_treatment.pdf',width=6,height = 10)

p.box.comb <- gridExtra::grid.arrange(p.boxplt.A260A230, p.boxplt.A260A280, p.boxplt.Nucleic_AcidnguL,nrow=1)
ggsave(plot = p.box.comb, file='/Users/sa01fd/DATA/HMW_DNA_Seaweeds/Method_development/boxplot_collage_h_filtered_per_treatment.pdf',width=12,height = 4)





nanodrop.csv %>% dplyr::left_join(nanodrop.meta, by=c('Sample_Name'='Sample_Name')) %>%
  dplyr::group_by(Strain) %>% arrange(-A260A230,.by_group=TRUE) %>%
  dplyr::select(Strain, A260A230, Protocol) %>%
  tidyr::pivot_wider(names_from=Protocol, values_from=A260A230) %>% collect() %>% data.frame


nanodrop.heatmap <- nanodrop.csv %>% dplyr::left_join(nanodrop.meta, by=c('Sample_Name'='Sample_Name')) %>%
  dplyr::group_by(Strain) %>% arrange(-A260A230,.by_group=TRUE) %>%
  dplyr::select(Strain, A260A230, Protocol) %>%
  tidyr::pivot_wider(names_from=Protocol, values_from=A260A230) %>% collect() %>% data.frame

rownames(nanodrop.heatmap) <- nanodrop.heatmap$Strain
nanodrop.heatmap <- nanodrop.heatmap[,2:ncol(nanodrop.heatmap)]

gplots::heatmap.2(as.matrix(nanodrop.heatmap),trace='none')

nanodrop.matrix

gplots::heatmap.2(as.matrix((nanodrop.matrix[,rev(colnames(nanodrop.matrix))])),trace='none', scale = c("row"),Colv=FALSE,dendrogram='row')


df.nan2 <- df.nan
rownames(df.nan2) <- df.nan2$Sample_Name
df.nan2[rownames(nanodrop.matrix),]
df.nan2[rownames(nanodrop.matrix),]$Protocol


#Pearson correlation
gplots::heatmap.2(as.matrix((nanodrop.matrix[selected.samples,])),
                  scale = c("row"),
                  Colv=FALSE,
                  dendrogram='row',
                  distfun = function(x) as.dist(1-cor(t(x))),
                  hclustfun = function(x) hclust(x, method="ward.D2"),
                  col=viridis::cividis(100),
                  RowSideColors=c('grey60','grey30','firebrick','deepskyblue3')[df.nan2[rownames(nanodrop.matrix[selected.samples,]),]$Protocol],
                  trace='none'
                  )







ggplot2::scale_color_manual(values =  )


nanodrop.csv %>% left_join(nanodrop.meta, by=c('Sample_Name'='Sample_Name')) %>%
  #filter(Sample_Name %in% selected.samples) %>%
  ggplot2::ggplot(aes(x=Protocol,y=A260A230, color=Taxon,group=Strain))+
  ggplot2::geom_line() +
  ggplot2::geom_point(aes(x=Protocol,y=A260A230, color=Taxon)) +
  ggplot2::xlab('Protocol')+
  ggplot2::ylab('A260/A230')+
  ggplot2::theme_classic()+fdb_style(aspect.ratio=1)+
  ggplot2::facet_wrap(~Taxon)+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))



p.lineplt <- nanodrop.csv %>% left_join(nanodrop.meta, by=c('Sample_Name'='Sample_Name')) %>%
  #filter(Sample_Name %in% selected.samples) %>%
  ggplot2::ggplot(aes(x=PEG,y=A260A230, color=Taxon,group=Strain))+
  ggplot2::geom_line() +
  ggplot2::geom_point(aes(x=PEG,y=A260A230, color=Taxon)) +
  ggplot2::xlab('Protocol')+
  ggplot2::ylab('A260/A230')+
  ggplot2::theme_classic()+fdb_style(aspect.ratio=1)+
  ggplot2::facet_wrap(~Taxon+NaCl,nrow=1)+
  ggplot2::theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
p.lineplt
ggsave(plot = p.lineplt, file='/Users/sa01fd/DATA/HMW_DNA_Seaweeds/Method_development/pre-post-PEG_per_Taxon.pdf',width=10,height = 4)



nanodrop.csv %>% left_join(nanodrop.meta, by=c('Sample_Name'='Sample_Name')) %>%
  filter(Nucleic_AcidnguL < 250) %>%
  ggplot2::ggplot(aes(x=PEG,y=Nucleic_AcidnguL, color=Taxon,group=Strain))+
  ggplot2::geom_line() +
  ggplot2::geom_point(aes(x=PEG,y=Nucleic_AcidnguL, color=Taxon)) +
  ggplot2::xlab('Protocol')+
  ggplot2::ylab('A260/A230')+
  ggplot2::theme_classic()+fdb_style(aspect.ratio=1)+
  ggplot2::facet_wrap(~Taxon+NaCl,nrow=1)+
  ggplot2::theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


nanodrop.csv %>% left_join(nanodrop.meta, by=c('Sample_Name'='Sample_Name')) %>%
  filter(A260A280 < 3) %>%
  ggplot2::ggplot(aes(x=PEG,y=A260A280, color=Taxon,group=Strain))+
  ggplot2::geom_line() +
  ggplot2::geom_point(aes(x=PEG,y=A260A280, color=Taxon)) +
  ggplot2::xlab('Protocol')+
  ggplot2::ylab('A260/A230')+
  ggplot2::theme_classic()+fdb_style(aspect.ratio=1)+
  ggplot2::facet_wrap(~Taxon+NaCl,nrow=1)+
  ggplot2::theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))





nanodrop.csv %>% left_join(nanodrop.meta, by=c('Sample_Name'='Sample_Name')) %>%
  filter(Sample_Name %in% selected.samples) %>%
  ggplot2::ggplot(aes(x=Protocol,y=A260A280, color=Protocol,group=Strain))+
  ggplot2::geom_line() +
  ggplot2::geom_point(aes(x=Protocol,y=A260A280, color=Protocol)) +
  ggplot2::xlab('Protocol')+
  ggplot2::ylab('A260/A280')+
  ggplot2::theme_classic()+fdb_style(aspect.ratio=1)



nanodrop.csv %>% left_join(nanodrop.meta, by=c('Sample_Name'='Sample_Name')) %>%
  #filter(Sample_Name %in% selected.samples) %>%
  filter(Nucleic_AcidnguL < 250) %>%
  ggplot2::ggplot()+
  ggplot2::geom_point(aes(x=Nucleic_AcidnguL,y=A260A230, color=Protocol)) +
  ggplot2::xlab('nucleic acids (ng/uL)')+
  ggplot2::ylab('A260/A230')+
  ggplot2::theme_classic()+fdb_style(aspect.ratio=1)


nanodrop.csv %>% left_join(nanodrop.meta, by=c('Sample_Name'='Sample_Name')) %>%
  filter(A260A280 < 2.5) %>%
  ggplot2::ggplot(aes(x=A260A280,y=A260A230, color=Protocol))+
  ggplot2::geom_point() +
  stat_ellipse()+
  ggplot2::xlab('A260/A280')+
  ggplot2::ylab('A260/A230')+
  ggplot2::theme_classic()+fdb_style(aspect.ratio=1)




#=================================================
# 15th of may
#==================================================
input <- "dsDNA 5_14_2021 3_24_36"

nanodrop <- read_nanodrop(file = paste0(dir, input, " AM_table.tsv"))
nanodrop.csv <- read_nanodrop_csv(file = paste0(dir, input, " AM.csv"))
nanodrop.meta <- read.delim(file = paste0(dir, input, " metadata.txt"))

nanodrop.long <- nanodrop_long(df = nanodrop)
nanodrop.summary <- nanodrop_ratios(df=nanodrop)


selected.samples <- nanodrop.long %>% left_join(nanodrop.meta, by=c('Sample.ID'='Sample_Name')) %>%
  filter(!Sample.ID %in% c('Sample 6','Sample 21','Sample 23','Sample 24','Sample 25', 'Sample 10','Sample 15','Sample 28','Sample 29','Sample 19','Sample 30')) %>%
  filter(LifeStage=='Sporophyte') %>%
  filter(homogenisation == 'LN2') %>% select(Sample.ID) %>% unique() %>% pull %>% as.character()

p <- nanodrop.long %>% left_join(nanodrop.meta, by=c('Sample.ID'='Sample_Name')) %>%
  filter(Sample.ID %in% selected.samples) %>%

  ggplot2::ggplot()+
  ggplot2::geom_line(aes(x=variable,y=value, color=Description)) +
  ggplot2::geom_area(aes(x=variable,y=value, fill=Description),alpha=0.5) +
  ggplot2::facet_wrap(Organism~Description,scales = 'free') +
  geom_text(data= nanodrop.csv  %>% mutate(Sample.ID = Sample_Name) %>%
              left_join(nanodrop.meta, by=c('Sample.ID'='Sample_Name')) %>%
              filter(Sample.ID %in% selected.samples),
            aes(x=255,y=0.3, label = paste0('260/280=',round(A260A280,2),'\n260/230=',round(A260A230,2),'\n',round(Nucleic_AcidnguL,2),'ng/uL')))+
  geom_vline(xintercept=230,linetype='dashed',color='grey')+
  ggplot2::xlab('Wavelength (nm)')+
  ggplot2::ylab('Absorbance')+
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0),breaks = c(220,240,260,280,300,320,340))+
  ggplot2::theme_classic()+fdb_style(aspect.ratio=0.5)+
  scale_fill_manual(values = colorRampPalette(RColorBrewer::brewer.pal(8, "Set1"))(length(selected.samples)))+
  scale_color_manual(values = colorRampPalette(RColorBrewer::brewer.pal(8, "Set1"))(length(selected.samples)))+
  theme(legend.position = 'none')

p
ggsave('~/DATA/CTAB_in_eppendorf_sporophyte_LN2_14may.pdf',width=9,height = 9)


selected.samples <- nanodrop.long %>% left_join(nanodrop.meta, by=c('Sample.ID'='Sample_Name')) %>%
  filter(!Sample.ID %in% c('Sample 6','Sample 21','Sample 23','Sample 24','Sample 25', 'Sample 26','Sample 15','Sample 28','Sample 29','Sample 19','Sample 30')) %>%
  filter(LifeStage=='Gametophyte') %>%
  filter(homogenisation == 'pestle') %>% select(Sample.ID) %>% unique() %>% pull %>% as.character()


p <- nanodrop.long %>% left_join(nanodrop.meta, by=c('Sample.ID'='Sample_Name')) %>%
  filter(Sample.ID %in% selected.samples) %>%
  ggplot2::ggplot()+
  ggplot2::geom_line(aes(x=variable,y=value, color=Description)) +
  ggplot2::geom_area(aes(x=variable,y=value, fill=Description),alpha=0.5) +
  ggplot2::facet_wrap(Organism~Description,scales = 'free') +
  geom_text(data= nanodrop.csv  %>% mutate(Sample.ID = Sample_Name) %>%
              left_join(nanodrop.meta, by=c('Sample.ID'='Sample_Name')) %>%
              filter(Sample.ID %in% selected.samples),
            aes(x=255,y=0.3, label = paste0('260/280=',round(A260A280,2),'\n260/230=',round(A260A230,2),'\n',round(Nucleic_AcidnguL,2),'ng/uL')))+
  geom_vline(xintercept=230,linetype='dashed',color='grey')+
  ggplot2::xlab('Wavelength (nm)')+
  ggplot2::ylab('Absorbance')+
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0),breaks = c(220,240,260,280,300,320,340))+
  ggplot2::theme_classic()+fdb_style(aspect.ratio=0.5)+
  scale_fill_manual(values = colorRampPalette(RColorBrewer::brewer.pal(8, "Set1"))(length(selected.samples)))+
  scale_color_manual(values = colorRampPalette(RColorBrewer::brewer.pal(8, "Set1"))(length(selected.samples)))+theme(legend.position = 'none')

ggsave('~/DATA/CTAB_in_eppendorf_gametophyte_pestle_14may.pdf',width=5,height = 5)


selected.samples <- nanodrop.long %>% left_join(nanodrop.meta, by=c('Sample.ID'='Sample_Name')) %>%
  filter(!Sample.ID %in% c('Sample 6','Sample 21','Sample 23','Sample 24','Sample 25', 'Sample 26','Sample 15','Sample 28','Sample 29','Sample 19','Sample 30')) %>%
  filter(LifeStage=='Gametophyte') %>%
  filter(homogenisation == 'LN2') %>% select(Sample.ID) %>% unique() %>% pull %>% as.character()


p <- nanodrop.long %>% left_join(nanodrop.meta, by=c('Sample.ID'='Sample_Name')) %>%
  filter(Sample.ID %in% selected.samples) %>%
  ggplot2::ggplot()+
  ggplot2::geom_line(aes(x=variable,y=value, color=Description)) +
  ggplot2::geom_area(aes(x=variable,y=value, fill=Description),alpha=0.5) +
  ggplot2::facet_wrap(Organism~Description,scales = 'free') +
  geom_text(data= nanodrop.csv  %>% mutate(Sample.ID = Sample_Name) %>%
              left_join(nanodrop.meta, by=c('Sample.ID'='Sample_Name')) %>%
              filter(Sample.ID %in% selected.samples),
            aes(x=255,y=0.3, label = paste0('260/280=',round(A260A280,2),'\n260/230=',round(A260A230,2),'\n',round(Nucleic_AcidnguL,2),'ng/uL')))+
  geom_vline(xintercept=230,linetype='dashed',color='grey')+
  ggplot2::xlab('Wavelength (nm)')+
  ggplot2::ylab('Absorbance')+
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0),breaks = c(220,240,260,280,300,320,340))+
  ggplot2::theme_classic()+fdb_style(aspect.ratio=0.5)+
  scale_fill_manual(values = colorRampPalette(RColorBrewer::brewer.pal(8, "Set1"))(length(selected.samples)))+
  scale_color_manual(values = colorRampPalette(RColorBrewer::brewer.pal(8, "Set1"))(length(selected.samples)))+
  theme(legend.position = 'none')

ggsave('~/DATA/CTAB_in_eppendorf_gametophyte_LN2_14may.pdf',width=5,height = 5)





nanodrop.csv  %>%
  left_join(nanodrop.meta, by=c('Sample_Name'='Sample_Name')) %>%
  filter(!Sample_Name %in% c('Sample 6','Sample 21','Sample 23','Sample 24','Sample 25', 'Sample 26','Sample 15','Sample 28','Sample 29','Sample 19','Sample 30')) %>%
  #filter(Sample.ID %in% selected.samples) %>%
  ggplot(aes(x=Nucleic_AcidnguL,A260A230)) + geom_point(aes(color=Description)) + fdb_style()


nanodrop.csv  %>%
  left_join(nanodrop.meta, by=c('Sample_Name'='Sample_Name')) %>%
  filter(!Sample_Name %in% c('Sample 6','Sample 21','Sample 23','Sample 24','Sample 25', 'Sample 26','Sample 15','Sample 28','Sample 29','Sample 19','Sample 30')) %>%
  filter(Description %in% c('Ect_CTAB','NaCl_CTAB')) %>%
  mutate(NaCL = ifelse(grepl('Ect',Description),'n','y')) %>% #filter(Sample.ID %in% selected.samples) %>%
  ggplot(aes(x=NaCL,A260A230)) +  geom_boxplot(aes(color=Description,fill=Description),alpha=0.4) +
  geom_jitter(aes(color=Description))  + fdb_style(aspect.ratio = 2)+
  scale_fill_manual(values = colorRampPalette(RColorBrewer::brewer.pal(8, "Set1"))(4))+
  scale_color_manual(values = colorRampPalette(RColorBrewer::brewer.pal(8, "Set1"))(4))




#==========================================================================
input <- "dsDNA 5_10_2021 8_12_03"

nanodrop <- read_nanodrop(file = paste0(dir, input, " AM_table.tsv"))
nanodrop.csv <- read_nanodrop_csv(file = paste0(dir, input, " AM.csv"))
nanodrop.meta <- read.delim(file = paste0(dir, input, " metadata.txt"))

nanodrop.long <- nanodrop_long(df = nanodrop)
nanodrop.summary <- nanodrop_ratios(df=nanodrop)

selected.samples <- nanodrop.long %>% left_join(nanodrop.meta, by=c('Sample.ID'='Sample_Name')) %>%
  filter(!Sample.ID %in% c('Sample 1','Sample 2')) %>%
  filter(homogenisation == 'pestle') %>% select(Sample.ID) %>% unique() %>% pull %>% as.character()



p <- nanodrop.long %>% left_join(nanodrop.meta, by=c('Sample.ID'='Sample_Name')) %>%
  filter(Sample.ID %in% selected.samples) %>%

  ggplot2::ggplot()+
  ggplot2::geom_line(aes(x=variable,y=value, color=Description)) +
  ggplot2::geom_area(aes(x=variable,y=value, fill=Description),alpha=0.5) +
  ggplot2::facet_wrap(Organism~Description,scales = 'free') +
  geom_text(data= nanodrop.csv  %>% mutate(Sample.ID = Sample_Name) %>%
              left_join(nanodrop.meta, by=c('Sample.ID'='Sample_Name')) %>%
              filter(Sample.ID %in% selected.samples),
            aes(x=255,y=0.3, label = paste0('260/280=',round(A260A280,2),'\n260/230=',round(A260A230,2),'\n',round(Nucleic_AcidnguL,2),'ng/uL')))+
  geom_vline(xintercept=230,linetype='dashed',color='grey')+
  ggplot2::xlab('Wavelength (nm)')+
  ggplot2::ylab('Absorbance')+
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0),breaks = c(220,240,260,280,300,320,340))+
  ggplot2::theme_classic()+fdb_style(aspect.ratio=0.5)+
  scale_fill_manual(values = colorRampPalette(RColorBrewer::brewer.pal(8, "Set1"))(9))+
  scale_color_manual(values = colorRampPalette(RColorBrewer::brewer.pal(8, "Set1"))(9))+
  theme(legend.position = 'none')

p
ggsave(paste0('~/DATA/',input,'CTAB_in_eppendorf_gametophyte_pestle.pdf'),width=6,height = 5)







selected.samples <- nanodrop.long %>% left_join(nanodrop.meta, by=c('Sample.ID'='Sample_Name')) %>%
  filter(Sample.ID %in% c('Sample 1','Sample 2')) %>%
  filter(homogenisation == 'pestle') %>% select(Sample.ID) %>% unique() %>% pull %>% as.character()

p <- nanodrop.long %>%
  filter(value>=0) %>%
  left_join(nanodrop.meta, by=c('Sample.ID'='Sample_Name')) %>%
  filter(Sample.ID %in% selected.samples) %>%

  ggplot2::ggplot()+
  ggplot2::geom_line(aes(x=variable,y=value, color=Biomass)) +
  ggplot2::geom_area(aes(x=variable,y=value, fill=Biomass),alpha=0.5) +
  ggplot2::facet_wrap(Organism~Biomass,scales = 'free') +
  geom_text(data= nanodrop.csv  %>% mutate(Sample.ID = Sample_Name) %>%
              left_join(nanodrop.meta, by=c('Sample.ID'='Sample_Name')) %>%
              filter(Sample.ID %in% selected.samples),
            aes(x=255,y=0.3, label = paste0('260/280=',round(A260A280,2),'\n260/230=',round(A260A230,2),'\n',round(Nucleic_AcidnguL,2),'ng/uL')))+
  geom_vline(xintercept=230,linetype='dashed',color='grey')+
  ggplot2::xlab('Wavelength (nm)')+
  ggplot2::ylab('Absorbance')+
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0),breaks = c(220,240,260,280,300,320,340))+
  ggplot2::theme_classic()+fdb_style(aspect.ratio=0.5)+
  scale_fill_manual(values = colorRampPalette(RColorBrewer::brewer.pal(8, "Set1"))(3))+
  scale_color_manual(values = colorRampPalette(RColorBrewer::brewer.pal(8, "Set1"))(3))+
  theme(legend.position = 'none')

p
ggsave(paste0('~/DATA/',input,'CTAB_high_low_biomass_in_eppendorf_Gametophyte_peslte.pdf'),width=6,height = 5)




#On a third day I scaled-up the CTAB protocol in 50mL falcon tubes by increasing all reagents 10-fold. One tube was washed x3 with Sortbitol buffer, while the other one was not washed. In this case I used Maryam’s best growing (‘the fat one’) Saccharine Latissima gametophyte culture.
selected.samples <- nanodrop.long %>% left_join(nanodrop.meta, by=c('Sample.ID'='Sample_Name')) %>%
  filter(homogenisation == 'glass_pestle') %>%
  select(Sample.ID) %>%
  unique() %>% pull %>% as.character()

p <- nanodrop.long %>%
  filter(value>=0) %>%
  left_join(nanodrop.meta, by=c('Sample.ID'='Sample_Name')) %>%
  filter(Sample.ID %in% selected.samples) %>%
  ggplot2::ggplot()+
  ggplot2::geom_line(aes(x=variable,y=value, color=Description)) +
  ggplot2::geom_area(aes(x=variable,y=value, fill=Description),alpha=0.5) +
  ggplot2::facet_wrap(Organism~Description,scales = 'free') +
  geom_text(data= nanodrop.csv  %>% mutate(Sample.ID = Sample_Name) %>%
              left_join(nanodrop.meta, by=c('Sample.ID'='Sample_Name')) %>%
              filter(Sample.ID %in% selected.samples),
            aes(x=255,y=0.3, label = paste0('260/280=',round(A260A280,2),'\n260/230=',round(A260A230,2),'\n',round(Nucleic_AcidnguL,2),'ng/uL')))+
  geom_vline(xintercept=230,linetype='dashed',color='grey')+
  ggplot2::xlab('Wavelength (nm)')+
  ggplot2::ylab('Absorbance')+
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0),breaks = c(220,240,260,280,300,320,340))+
  ggplot2::theme_classic()+fdb_style(aspect.ratio=0.5)+
  scale_fill_manual(values = colorRampPalette(RColorBrewer::brewer.pal(8, "Set1"))(3))+
  scale_color_manual(values = colorRampPalette(RColorBrewer::brewer.pal(8, "Set1"))(3))+
  theme(legend.position = 'none')

p
ggsave(paste0('~/DATA/',input,'CTAB_50mLFalcon_glass_pestle_Gametophyte_SL.pdf'),width=6,height = 5)



selected.samples <- nanodrop.long %>% left_join(nanodrop.meta, by=c('Sample.ID'='Sample_Name')) %>%
  filter(homogenisation == 'glass_pestle') %>% select(Sample.ID) %>% unique() %>% pull %>% as.character()

p <- nanodrop.long %>%
  filter(value>=0) %>%
  left_join(nanodrop.meta, by=c('Sample.ID'='Sample_Name')) %>%
  filter(Sample.ID %in% selected.samples) %>%

  ggplot2::ggplot()+
  ggplot2::geom_line(aes(x=variable,y=value, color=Description)) +
  ggplot2::geom_area(aes(x=variable,y=value, fill=Description),position = "identity",alpha=0.5) +
  #ggplot2::facet_wrap(Organism~Description,scales = 'free') +
  geom_text(data= nanodrop.csv  %>% mutate(Sample.ID = Sample_Name) %>%
              left_join(nanodrop.meta, by=c('Sample.ID'='Sample_Name')) %>%
              filter(Sample.ID %in% selected.samples),
            aes(x=255,y=0.3, label = paste0('260/280=',round(A260A280,2),'\n260/230=',round(A260A230,2),'\n',round(Nucleic_AcidnguL,2),'ng/uL')))+
  geom_vline(xintercept=230,linetype='dashed',color='grey')+
  ggplot2::xlab('Wavelength (nm)')+
  ggplot2::ylab('Absorbance')+
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0),breaks = c(220,240,260,280,300,320,340))+
  ggplot2::theme_classic()+fdb_style(aspect.ratio=0.5)+
  scale_fill_manual(values = colorRampPalette(RColorBrewer::brewer.pal(8, "Set1"))(10))+
  scale_color_manual(values = colorRampPalette(RColorBrewer::brewer.pal(8, "Set1"))(10))+
  theme(legend.position = 'none')

p
ggsave(paste0('~/DATA/',input,'CTAB_50mLFalcon_glass_pestle_Gametophyte_SL_overlayed.pdf'),width=4,height = 4)


#==============================================================================
# 11th of may
#==============================================================================

input <- "dsDNA 5_11_2021 10_44_03"
nanodrop <- read_nanodrop(file = paste0(dir, input, " AM_table.tsv"))
nanodrop.csv <- read_nanodrop_csv(file = paste0(dir, input, " AM.csv"))
nanodrop.meta <- read.delim(file = paste0(dir, input, " metadata.txt"))

input2 <- "dsDNA 5_10_2021 8_12_03"
nanodrop2 <- read_nanodrop(file = paste0(dir, input2, " AM_table.tsv"))
nanodrop.csv2 <- read_nanodrop_csv(file = paste0(dir, input2, " AM.csv"))
nanodrop.meta2 <- read.delim(file = paste0(dir, input2, " metadata.txt"))

nanodrop <- rbind(nanodrop2[c(1,2),],nanodrop[c(1,3),])
nanodrop.csv <- rbind(nanodrop.csv2[c(1,2),],nanodrop.csv[c(1,3),])
nanodrop.meta <- rbind(nanodrop.meta2[c(1,2),],nanodrop.meta[c(1,3),])

nanodrop$Sample.ID = c('Sample 1','Sample 2', 'Sample 3', 'Sample 4')
nanodrop.csv$Sample_Name = c('Sample 1','Sample 2', 'Sample 3', 'Sample 4')
nanodrop.meta$Sample_Name = c('Sample 1','Sample 2', 'Sample 3', 'Sample 4')

nanodrop.long <- nanodrop_long(df = nanodrop)
nanodrop.summary <- nanodrop_ratios(df=nanodrop)

selected.samples <- nanodrop.long %>% left_join(nanodrop.meta, by=c('Sample.ID'='Sample_Name')) %>%
  filter(homogenisation == 'pestle') %>% select(Sample.ID) %>% unique() %>% pull %>% as.character()

p <- nanodrop.long %>% left_join(nanodrop.meta, by=c('Sample.ID'='Sample_Name')) %>%
  ggplot2::ggplot()+
  ggplot2::geom_line(aes(x=variable,y=value, color=Description)) +
  ggplot2::geom_area(aes(x=variable,y=value, fill=Description),alpha=0.5) +
  ggplot2::facet_wrap(Biomass~ Description,scales = 'free') +
  geom_text(data= nanodrop.csv  %>% mutate(Sample.ID = Sample_Name) %>%
              left_join(nanodrop.meta, by=c('Sample.ID'='Sample_Name')) %>%
              filter(Sample.ID %in% selected.samples),
            aes(x=255,y=0.3, label = paste0('260/280=',round(A260A280,2),'\n260/230=',round(A260A230,2),'\n',round(Nucleic_AcidnguL,2),'ng/uL')))+
  geom_vline(xintercept=230,linetype='dashed',color='grey')+
  ggplot2::xlab('Wavelength (nm)')+
  ggplot2::ylab('Absorbance')+
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0),breaks = c(220,240,260,280,300,320,340))+
  ggplot2::theme_classic()+fdb_style(aspect.ratio=0.5)+
  scale_fill_manual(values = colorRampPalette(RColorBrewer::brewer.pal(8, "Set1"))(9))+
  scale_color_manual(values = colorRampPalette(RColorBrewer::brewer.pal(8, "Set1"))(9))+
  theme(legend.position = 'none')

p
ggsave(paste0('~/DATA/',input,'CTAB_post_cleanup.pdf'),width=5,height = 5)


#==============================================================================


input <- "dsDNA 5_12_2021 7_49_09"
nanodrop <- read_nanodrop(file = paste0(dir, input, " AM_table.tsv"))
nanodrop.csv <- read_nanodrop_csv(file = paste0(dir, input, " AM.csv"))
nanodrop.meta <- read.delim(file = paste0(dir, input, " metadata.txt"))




nanodrop.long <- nanodrop_long(df = nanodrop)
nanodrop.summary <- nanodrop_ratios(df=nanodrop)

selected.samples <- nanodrop.long %>%
  filter(!(Sample.ID %in% c('Sample 7','Sample 8','Sample 9', 'Sample 10','Sample 11'))) %>%
  left_join(nanodrop.meta, by=c('Sample.ID'='Sample_Name')) %>%
  filter(extractionVolume=='Eppendorf') %>%
  filter(homogenisation == 'LN2_mortar') %>% select(Sample.ID) %>% unique() %>% pull %>% as.character()

p <- nanodrop.long %>% left_join(nanodrop.meta, by=c('Sample.ID'='Sample_Name')) %>%
  filter(Sample.ID %in% selected.samples) %>%
  ggplot2::ggplot()+
  ggplot2::geom_line(aes(x=variable,y=value, color=Description)) +
  ggplot2::geom_area(aes(x=variable,y=value, fill=Description),alpha=0.5) +
  ggplot2::facet_wrap(Biomass~ Description+Sample.ID,scales = 'free') +
  geom_text(data= nanodrop.csv  %>% mutate(Sample.ID = Sample_Name) %>%
              left_join(nanodrop.meta, by=c('Sample.ID'='Sample_Name')) %>%
              filter(Sample.ID %in% selected.samples),
            aes(x=255,y=0.3, label = paste0('260/280=',round(A260A280,2),'\n260/230=',round(A260A230,2),'\n',round(Nucleic_AcidnguL,2),'ng/uL')))+
  geom_vline(xintercept=230,linetype='dashed',color='grey')+
  ggplot2::xlab('Wavelength (nm)')+
  ggplot2::ylab('Absorbance')+
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0),breaks = c(220,240,260,280,300,320,340))+
  ggplot2::theme_classic()+fdb_style(aspect.ratio=0.5)+
  scale_fill_manual(values = colorRampPalette(RColorBrewer::brewer.pal(8, "Set1"))(9))+
  scale_color_manual(values = colorRampPalette(RColorBrewer::brewer.pal(8, "Set1"))(9))+
  theme(legend.position = 'none')

p
ggsave(paste0('~/DATA/',input,'CTAB_test_EthWash.pdf'),width=9,height = 9)




selected.samples <- nanodrop.long %>%
  filter(Sample.ID %in% c('Sample 7','Sample 11')) %>%
  left_join(nanodrop.meta, by=c('Sample.ID'='Sample_Name')) %>%
  filter(extractionVolume=='50mL falcon') %>%
  filter(homogenisation == 'LN2_mortar') %>% select(Sample.ID) %>% unique() %>% pull %>% as.character()

p <- nanodrop.long %>% left_join(nanodrop.meta, by=c('Sample.ID'='Sample_Name')) %>%
  filter(Sample.ID %in% selected.samples) %>%
  ggplot2::ggplot()+
  ggplot2::geom_line(aes(x=variable,y=value, color=Description)) +
  ggplot2::geom_area(aes(x=variable,y=value, fill=Description),alpha=0.5) +
  ggplot2::facet_wrap(Biomass~ Description,scales = 'free') +
  geom_text(data= nanodrop.csv  %>% mutate(Sample.ID = Sample_Name) %>%
              left_join(nanodrop.meta, by=c('Sample.ID'='Sample_Name')) %>%
              filter(Sample.ID %in% selected.samples),
            aes(x=255,y=0.3, label = paste0('260/280=',round(A260A280,2),'\n260/230=',round(A260A230,2),'\n',round(Nucleic_AcidnguL,2),'ng/uL')))+
  geom_vline(xintercept=230,linetype='dashed',color='grey')+
  ggplot2::xlab('Wavelength (nm)')+
  ggplot2::ylab('Absorbance')+
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0),breaks = c(220,240,260,280,300,320,340))+
  ggplot2::theme_classic()+fdb_style(aspect.ratio=0.5)+
  scale_fill_manual(values = colorRampPalette(RColorBrewer::brewer.pal(8, "Set1"))(9))+
  scale_color_manual(values = colorRampPalette(RColorBrewer::brewer.pal(8, "Set1"))(9))+
  theme(legend.position = 'none')

p
ggsave(paste0('~/DATA/',input,'CTAB_NaCl_falcon_AL_Sporophyte.pdf'),width=5,height = 5)




selected.samples <- nanodrop.long %>%
  filter(Sample.ID %in% c('Sample 9','Sample 10')) %>%
  left_join(nanodrop.meta, by=c('Sample.ID'='Sample_Name')) %>%
  filter(extractionVolume=='Eppendorf') %>%
  filter(homogenisation == 'pestle') %>% select(Sample.ID) %>% unique() %>% pull %>% as.character()


p <-nanodrop.long %>% left_join(nanodrop.meta, by=c('Sample.ID'='Sample_Name')) %>%
  filter(Sample.ID %in% selected.samples) %>%
  ggplot2::ggplot()+
  ggplot2::geom_line(aes(x=variable,y=value, color=Description)) +
  ggplot2::geom_area(aes(x=variable,y=value, fill=Description),alpha=0.5) +
  ggplot2::facet_wrap(Biomass~ Description,scales = 'free') +
  geom_text(data= nanodrop.csv  %>% mutate(Sample.ID = Sample_Name) %>%
              left_join(nanodrop.meta, by=c('Sample.ID'='Sample_Name')) %>%
              filter(Sample.ID %in% selected.samples),
            aes(x=255,y=0.3, label = paste0('260/280=',round(A260A280,2),'\n260/230=',round(A260A230,2),'\n',round(Nucleic_AcidnguL,2),'ng/uL')))+
  geom_vline(xintercept=230,linetype='dashed',color='grey')+
  ggplot2::xlab('Wavelength (nm)')+
  ggplot2::ylab('Absorbance')+
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0),breaks = c(220,240,260,280,300,320,340))+
  ggplot2::theme_classic()+fdb_style(aspect.ratio=0.5)+
  scale_fill_manual(values = colorRampPalette(RColorBrewer::brewer.pal(8, "Set1"))(9))+
  scale_color_manual(values = colorRampPalette(RColorBrewer::brewer.pal(8, "Set1"))(9))+
  theme(legend.position = 'none')

p





#============================================

input <- 'dsDNA 5_19_2021 6_09_13'

nanodrop <- read_nanodrop(file = paste0(dir, input, " AM_table.tsv"))
nanodrop.csv <- read_nanodrop_csv(file = paste0(dir, input, " AM.csv"))
nanodrop.meta <- read.delim(file = paste0(dir, input, " metadata.txt"))

nanodrop.long <- nanodrop_long(df = nanodrop)
nanodrop.summary <- nanodrop_ratios(df=nanodrop)

selected.samples <- nanodrop.long %>%
  filter(!(Sample.ID %in% c('Sample 1','Sample 13','Sample 4'))) %>%
  left_join(nanodrop.meta, by=c('Sample.ID'='Sample_Name')) %>%
  filter(Organism !='Saccharina latissima') %>%
  filter(extractionVolume=='Eppendorf') %>%
  filter(homogenisation == 'pestle') %>% select(Sample.ID) %>% unique() %>% pull %>% as.character()

p <- nanodrop.long %>% left_join(nanodrop.meta, by=c('Sample.ID'='Sample_Name')) %>%
  filter(Sample.ID %in% selected.samples) %>%
  ggplot2::ggplot()+
  ggplot2::geom_line(aes(x=variable,y=value, color=Description)) +
  ggplot2::geom_area(aes(x=variable,y=value, fill=Description),alpha=0.5) +
  ggplot2::facet_wrap(Biomass~ Organism+Strain+Sample.ID,scales = 'free') +
  geom_text(data= nanodrop.csv  %>% mutate(Sample.ID = Sample_Name) %>%
              left_join(nanodrop.meta, by=c('Sample.ID'='Sample_Name')) %>%
              filter(Sample.ID %in% selected.samples),
            aes(x=290,y=2, label = paste0('260/280=',round(A260A280,2),'\n260/230=',round(A260A230,2),'\n',round(Nucleic_AcidnguL,2),'ng/uL')))+
  geom_vline(xintercept=230,linetype='dashed',color='grey')+
  ggplot2::xlab('Wavelength (nm)')+
  ggplot2::ylab('Absorbance')+
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0),breaks = c(220,240,260,280,300,320,340))+
  ggplot2::theme_classic()+fdb_style(aspect.ratio=0.5)+
  scale_fill_manual(values = c("#47A265","#E41A1C","#7B7281"))+
  scale_color_manual(values = c("#47A265","#E41A1C","#7B7281"))+
  theme(legend.position = 'none')

p
ggsave(paste0('~/DATA/',input,'CTAB_NaCl_reds_browns_and_green.pdf'),width=9,height = 9)





selected.samples <- nanodrop.long %>%
  filter((Sample.ID %in% c('Sample 7','Sample 8','Sample 9','Sample 10'))) %>%
  left_join(nanodrop.meta, by=c('Sample.ID'='Sample_Name')) %>%
  filter(Organism =='Saccharina latissima') %>%
  filter(extractionVolume=='Eppendorf') %>%
  filter(homogenisation == 'pestle') %>% select(Sample.ID) %>% unique() %>% pull %>% as.character()

p <- nanodrop.long %>% left_join(nanodrop.meta, by=c('Sample.ID'='Sample_Name')) %>%
  filter(Sample.ID %in% selected.samples) %>%
  ggplot2::ggplot()+
  ggplot2::geom_line(aes(x=variable,y=value, color=Description)) +
  ggplot2::geom_area(aes(x=variable,y=value, fill=Description),alpha=0.5) +
  ggplot2::facet_wrap(Biomass~ Organism+Strain+Sample.ID,scales = 'free') +
  geom_text(data= nanodrop.csv  %>% mutate(Sample.ID = Sample_Name) %>%
              left_join(nanodrop.meta, by=c('Sample.ID'='Sample_Name')) %>%
              filter(Sample.ID %in% selected.samples),
            aes(x=255,y=0.3, label = paste0('260/280=',round(A260A280,2),'\n260/230=',round(A260A230,2),'\n',round(Nucleic_AcidnguL,2),'ng/uL')))+
  geom_vline(xintercept=230,linetype='dashed',color='grey')+
  ggplot2::xlab('Wavelength (nm)')+
  ggplot2::ylab('Absorbance')+
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0),breaks = c(220,240,260,280,300,320,340))+
  ggplot2::theme_classic()+fdb_style(aspect.ratio=0.5)+
  scale_fill_manual(values = c("#47A265","#E41A1C","#7B7281"))+
  scale_color_manual(values = c("#47A265","#E41A1C","#7B7281"))+
  theme(legend.position = 'none')

p

ggsave(paste0('~/DATA/',input,'CTAB_NaCl_SL_metagenome_samples_dave.pdf'),width=5,height = 5)






#==============

input <- 'dsDNA 5_20_2021 5_54_51'

nanodrop <- read_nanodrop(file = paste0(dir, input, " AM_table.tsv"))
nanodrop.csv <- read_nanodrop_csv(file = paste0(dir, input, " AM.csv"))
nanodrop.meta <- read.delim(file = paste0(dir, input, " metadata.txt"))

input2 <- "dsDNA 5_19_2021 6_09_13"
nanodrop2 <- read_nanodrop(file = paste0(dir, input2, " AM_table.tsv"))
nanodrop.csv2 <- read_nanodrop_csv(file = paste0(dir, input2, " AM.csv"))
nanodrop.meta2 <- read.delim(file = paste0(dir, input2, " metadata.txt"))

sel.row.1 <- c(1,2,3,4)
sel.row.2 <- c(1,4,9,10)
#sel.row.2 <- c(11,4,9,10)

nanodrop <- rbind(nanodrop2[sel.row.2,],nanodrop[sel.row.1,])
nanodrop.csv <- rbind(nanodrop.csv2[sel.row.2,],nanodrop.csv[sel.row.1,])
nanodrop.meta <- rbind(nanodrop.meta2[sel.row.2,],nanodrop.meta[sel.row.1,])

nanodrop$Sample.ID = paste('Sample',1:nrow(nanodrop))
nanodrop.csv$Sample_Name = paste('Sample',1:nrow(nanodrop))
nanodrop.meta$Sample_Name = paste('Sample',1:nrow(nanodrop))



nanodrop.long <- nanodrop_long(df = nanodrop)
nanodrop.summary <- nanodrop_ratios(df=nanodrop)

selected.samples <- nanodrop.long %>%
  filter((Sample.ID %in% c('Sample 1','Sample 2','Sample 5','Sample 6'))) %>%
  left_join(nanodrop.meta, by=c('Sample.ID'='Sample_Name')) %>%
  select(Sample.ID) %>% unique() %>% pull %>% as.character()

p <- nanodrop.long %>% left_join(nanodrop.meta, by=c('Sample.ID'='Sample_Name')) %>%
  filter(Sample.ID %in% selected.samples) %>%
  ggplot2::ggplot()+
  ggplot2::geom_line(aes(x=variable,y=value, color=Description)) +
  ggplot2::geom_area(aes(x=variable,y=value, fill=Description),alpha=0.5) +
  ggplot2::facet_wrap(Biomass~ Organism+Strain+Sample.ID,scales = 'free') +
  geom_text(data= nanodrop.csv  %>% mutate(Sample.ID = Sample_Name) %>%
              left_join(nanodrop.meta, by=c('Sample.ID'='Sample_Name')) %>%
              filter(Sample.ID %in% selected.samples),
            aes(x=290,y=0.3, label = paste0('260/280=',round(A260A280,2),'\n260/230=',round(A260A230,2),'\n',round(Nucleic_AcidnguL,2),'ng/uL')))+
  geom_vline(xintercept=230,linetype='dashed',color='grey')+
  ggplot2::xlab('Wavelength (nm)')+
  ggplot2::ylab('Absorbance')+
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0),breaks = c(220,240,260,280,300,320,340))+
  ggplot2::theme_classic()+fdb_style(aspect.ratio=0.5)+
  scale_fill_manual(values = c("#47A265","#E41A1C","#7B7281"))+
  scale_color_manual(values = c("#47A265","#E41A1C","#7B7281"))+
  theme(legend.position = 'none')

p
ggsave(paste0('~/DATA/',input,'CTAB_CCAP_seaweeds_NaCl_wash2.pdf'),width=7,height = 7)



selected.samples <- nanodrop.long %>%
  filter((Sample.ID %in% c('Sample 1','Sample 5'))) %>%
  left_join(nanodrop.meta, by=c('Sample.ID'='Sample_Name')) %>%
  select(Sample.ID) %>% unique() %>% pull %>% as.character()

p <- nanodrop.long %>% left_join(nanodrop.meta, by=c('Sample.ID'='Sample_Name')) %>%
  filter(Sample.ID %in% selected.samples) %>%
  ggplot2::ggplot()+
  ggplot2::geom_line(aes(x=variable,y=value, color=Sample.ID)) +
  ggplot2::geom_area(aes(x=variable,y=value, fill=Sample.ID),alpha=0.5,position = "identity") +
  #ggplot2::facet_wrap(Biomass~ Organism+Strain+Sample.ID,scales = 'free') +
  geom_text(data= nanodrop.csv  %>% mutate(Sample.ID = Sample_Name) %>%
              left_join(nanodrop.meta, by=c('Sample.ID'='Sample_Name')) %>%
              filter(Sample.ID %in% selected.samples),
            aes(x=290,y=0.3, label = paste0('260/280=',round(A260A280,2),'\n260/230=',round(A260A230,2),'\n',round(Nucleic_AcidnguL,2),'ng/uL')))+
  geom_vline(xintercept=230,linetype='dashed',color='grey')+
  ggplot2::xlab('Wavelength (nm)')+
  ggplot2::ylab('Absorbance')+
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0),breaks = c(220,240,260,280,300,320,340))+
  ggplot2::theme_classic()+fdb_style(aspect.ratio=0.5)+
  scale_fill_manual(values = c("#47A265","#E41A1C","#7B7281"))+
  scale_color_manual(values = c("#47A265","#E41A1C","#7B7281"))+
  theme(legend.position = 'none')

p
ggsave(paste0('~/DATA/',input,'CTAB_CCAP_sporh_overlay_NaCl_wash.pdf'),width=5,height = 5)



selected.samples <- nanodrop.long %>%
  filter((Sample.ID %in% c('Sample 3','Sample 4','Sample 7','Sample 8'))) %>%
  left_join(nanodrop.meta, by=c('Sample.ID'='Sample_Name')) %>%
  select(Sample.ID) %>% unique() %>% pull %>% as.character()

p <- nanodrop.long %>% left_join(nanodrop.meta, by=c('Sample.ID'='Sample_Name')) %>%
  filter(Sample.ID %in% selected.samples) %>%
  ggplot2::ggplot()+
  ggplot2::geom_line(aes(x=variable,y=value, color=Description)) +
  ggplot2::geom_area(aes(x=variable,y=value, fill=Description),alpha=0.5) +
  ggplot2::facet_wrap(Biomass~ Organism+Strain+Sample.ID,scales = 'free') +
  geom_text(data= nanodrop.csv  %>% mutate(Sample.ID = Sample_Name) %>%
              left_join(nanodrop.meta, by=c('Sample.ID'='Sample_Name')) %>%
              filter(Sample.ID %in% selected.samples),
            aes(x=290,y=0.3, label = paste0('260/280=',round(A260A280,2),'\n260/230=',round(A260A230,2),'\n',round(Nucleic_AcidnguL,2),'ng/uL')))+
  geom_vline(xintercept=230,linetype='dashed',color='grey')+
  ggplot2::xlab('Wavelength (nm)')+
  ggplot2::ylab('Absorbance')+
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0),breaks = c(220,240,260,280,300,320,340))+
  ggplot2::theme_classic()+fdb_style(aspect.ratio=0.5)+
  scale_fill_manual(values = c("#47A265","#E41A1C","#7B7281"))+
  scale_color_manual(values = c("#47A265","#E41A1C","#7B7281"))+
  theme(legend.position = 'none')

p
ggsave(paste0('~/DATA/',input,'CTAB_CCAP_metagenome_nacl_wash.pdf'),width=7,height = 7)


#==============







p <- nanodrop.long %>% #filter(Sample.ID %in% selected.samples) %>%
  ggplot2::ggplot()+
  ggplot2::geom_line(aes(x=variable,y=value, color=Biomass)) +
  ggplot2::geom_area(aes(x=variable,y=value, fill=Biomass),alpha=0.5) +
  ggplot2::facet_wrap(~Biomass,scales = 'free') +
  geom_text(data= nanodrop.summary  %>%
              filter(Sample.ID %in% selected.samples),
            aes(x=255,y=0.3, label = paste0('260/280=',round(r260_280,2),'; 260/230=',round(r260_230,2))))+
  geom_vline(xintercept=230,linetype='dashed',color='grey')+
  ggplot2::xlab('Wavelength (nm)')+
  ggplot2::ylab('Absorbance')+
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0),breaks = c(220,240,260,280,300,320,340)) +
  ggplot2::theme_classic()+
  fdb_style(aspect.ratio=0.5)+
  scale_fill_manual(values = colorRampPalette(RColorBrewer::brewer.pal(8, "Set1"))(length(selected.samples)))+
  scale_color_manual(values = colorRampPalette(RColorBrewer::brewer.pal(8, "Set1"))(length(selected.samples)))


p


nanodrop.long %>%# filter(Sample.ID %in% c('Sample 1','Sample 2')) %>%
  ggplot2::ggplot(aes(x=variable,y=value, color=Sample.ID))+
  ggplot2::geom_line()+
  geom_vline(xintercept=230,linetype='dashed',color='grey')+
  ggplot2::xlab('Wavelength (nm)')+
  ggplot2::ylab('Absorbance')+
  #scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0))+
  ggplot2::theme_classic()+theme(aspect.ratio=0.5)

nanodrop.long %>% filter(Sample.ID %in% c('Sample 1','Sample 2')) %>%
  ggplot2::ggplot(aes(x=variable,y=value, color=Sample.ID))+
  ggplot2::geom_line()+
  geom_vline(xintercept=230,linetype='dashed',color='grey')+
  ggplot2::xlab('Wavelength (nm)')+
  ggplot2::ylab('Absorbance')+
  #scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0))+
  ggplot2::theme_classic()+theme(aspect.ratio=0.5)+facet_wrap(~Sample.ID,scales='free')


selected.samples <- c('Sample 1','Sample 2')
selected.samples <- c('Sample 7','Sample 8')
selected.samples <- c('Sample 5','Sample 6')

selected.samples <- c('Sample 1','Sample 2','Sample 3','Sample 4')
selected.samples <- c('Sample 5','Sample 6','Sample 7','Sample 8')
selected.samples <- c('Sample 8','Sample 9','Sample 10','Sample 11')


selected.samples <- nanodrop.long %>% select(Sample.ID) %>% unique() %>% pull %>% as.character()



selected.samples <- c('Sample 10','Sample 23','Sample 24','Sample 25','Sample 26')
selected.samples <- c('Sample 16','Sample 28','Sample 29')
selected.samples <- c('Sample 6','Sample 21','Sample 22')
selected.samples <- c('Sample 15','Sample 27')
selected.samples <- c('Sample 19','Sample 30','Sample 31')


p <- nanodrop.long %>% filter(Sample.ID %in% selected.samples) %>%
  ggplot2::ggplot()+
  ggplot2::geom_line(aes(x=variable,y=value, color=Sample.ID)) +
  ggplot2::geom_area(aes(x=variable,y=value, fill=Sample.ID),alpha=0.5) +
  ggplot2::facet_wrap(~Sample.ID,scales = 'free') +
  geom_text(data= nanodrop.summary  %>%
                    filter(Sample.ID %in% selected.samples),
            aes(x=255,y=0.3, label = paste0('260/280=',round(r260_280,2),'; 260/230=',round(r260_230,2))))+
  geom_vline(xintercept=230,linetype='dashed',color='grey')+
  ggplot2::xlab('Wavelength (nm)')+
  ggplot2::ylab('Absorbance')+
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0),breaks = c(220,240,260,280,300,320,340)) +
  ggplot2::theme_classic()+fdb_style(aspect.ratio=0.5)+
  scale_fill_manual(values = colorRampPalette(RColorBrewer::brewer.pal(8, "Set1"))(length(selected.samples)))+
  scale_color_manual(values = colorRampPalette(RColorBrewer::brewer.pal(8, "Set1"))(length(selected.samples)))


p

ggsave('~/DATA/CTAB_in_eppendorf_preliminary_exported_12may...pdf',width=9,height = 9)





nanodrop.long %>% filter(Sample.ID %in% c('Sample 1','Sample 3')) %>%
  ggplot2::ggplot(aes(x=variable,y=value, color=Date.and.Time))+
  ggplot2::geom_line()+
  ggplot2::theme_classic()+facet_wrap(~Sample.ID)


p <- nanodrop.long %>% filter(Sample.ID %in% c('Sample 1','Sample 2')) %>%
  ggplot2::ggplot(aes(x=variable,y=value, color=Sample.ID))+
  ggplot2::geom_line()+
  geom_vline(xintercept=230,linetype='dashed',color='grey')+
  ggplot2::xlab('Wavelength (nm)')+
  ggplot2::ylab('Absorbance')+
  #scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0))+
  ggplot2::theme_classic()+theme(aspect.ratio=0.5)

ggsave('~/DATA/CTAB_in_eppendorf_preliminary.pdf',width=6,height = 5)

p <- nanodrop.long %>%
  ggplot2::ggplot()+
  ggplot2::geom_line(aes(x=variable,y=value, color=Sample.ID)) +
  ggplot2::facet_wrap(~Sample.ID,scales = 'free') +
  geom_text(data=nanodrop.summary, aes(x=320,y=1, label = paste0('260/280=',round(r260_280,2),'\n 260/230=',round(r260_230,2))))+
  #geom_vline(xintercept=280,linetype='dashed',color='grey')+
  #geom_vline(xintercept=260,linetype='dashed',color='grey')+
  geom_vline(xintercept=230,linetype='dashed',color='grey')+
  ggplot2::xlab('Wavelength (nm)')+
  ggplot2::ylab('Absorbance')+
  #scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  ggplot2::theme_classic()+theme(aspect.ratio=0.5)



  ggsave('~/DATA/all_preliminary.pdf',width=9,height = 9)


  p <- nanodrop.long %>% filter(Sample.ID %in% c('Sample 1','Sample 3'))  %>%
    ggplot2::ggplot()+
    ggplot2::geom_line(aes(x=variable,y=value, color=Sample.ID)) +
    ggplot2::facet_wrap(~Sample.ID,scales = 'free') +
    geom_text(data=nanodrop.summary%>% filter(Sample.ID %in% c('Sample 1','Sample 3')), aes(x=320,y=0.3, label = paste0('260/280=',round(r260_280,2),'\n 260/230=',round(r260_230,2))))+
    #geom_vline(xintercept=280,linetype='dashed',color='grey')+
    #geom_vline(xintercept=260,linetype='dashed',color='grey')+
    geom_vline(xintercept=230,linetype='dashed',color='grey')+
    ggplot2::xlab('Wavelength (nm)')+
    ggplot2::ylab('Absorbance')+
    #scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    ggplot2::theme_classic()+theme(aspect.ratio=0.5)

ggsave('~/DATA/clean_facet_preliminary.pdf',width=9,height = 9)


p<-nanodrop.long %>% filter(Sample.ID %in% c('Sample 7','Sample 8'))%>%
  ggplot2::ggplot(aes(x=variable,y=value, color=Sample.ID))+
  ggplot2::geom_line()+
  ggplot2::theme_classic()+geom_vline(xintercept=230,linetype='dashed',color='grey')+
  ggplot2::xlab('Wavelength (nm)')+
  ggplot2::ylab('Absorbance')+
  #scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  ggplot2::theme_classic()+theme(aspect.ratio=0.5)

ggsave('~/DATA/sorbitol_preliminary.pdf',width=9,height = 9)



nanodrop.long %>% filter(Sample.ID %in% c('Sample 1','Sample 2','Sample 3'))%>%
  ggplot2::ggplot(aes(x=variable,y=value, color=Sample.ID))+
  ggplot2::geom_line()+
  ggplot2::theme_classic()+geom_vline(xintercept=230,linetype='dashed',color='grey')+
  ggplot2::xlab('Wavelength (nm)')+
  ggplot2::ylab('Absorbance')+
  #scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  ggplot2::theme_classic()

p <- nanodrop.long %>% filter(Sample.ID %in% c('Sample 3','Sample 4','Sample 5','Sample 6'))%>%
  ggplot2::ggplot(aes(x=variable,y=value, color=Sample.ID))+
  ggplot2::geom_line()+
  ggplot2::theme_classic()+geom_vline(xintercept=230,linetype='dashed',color='grey')+
  ggplot2::xlab('Wavelength (nm)')+
  ggplot2::ylab('Absorbance')+
  #scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  ggplot2::theme_classic()+theme(aspect.ratio=0.5)


ggsave('~/DATA/microwave_sorbitol_preliminary.pdf',width=6,height = 5)


nanodrop.long %>% filter(Sample.ID %in% c('Sample 1','Sample 3'))%>%
  ggplot2::ggplot(aes(x=variable,y=value, color=Sample.ID))+
  ggplot2::geom_line()+
  ggplot2::theme_classic()+geom_vline(xintercept=230,linetype='dashed',color='grey')+
  ggplot2::xlab('Wavelength (nm)')+
  ggplot2::ylab('Absorbance')+
  #scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  ggplot2::theme_classic()


