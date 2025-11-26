library(phyloseq)
library(picante) #for pd()
library(vegan)
library(MuMIn)
library(DHARMa)
library(ggsignif)
library(scales)
library(ggplot2)
library(tidyr)
library(dplyr)
library(purrr)
library(forcats)

capture_dat <- read.delim("~/Dropbox/PanamaBats/Bat2017_02_27_original_Panama_Metadata.csv",
                          header = TRUE,
                          stringsAsFactors = FALSE)
capture_dat<-subset(capture_dat, species=="C.perspicillata" | species=="A.jamaicensis")
capture_dat<-subset(capture_dat, astro!="NA")
capture_dat$location<-as.factor(capture_dat$location)
levels(capture_dat$location) <- c("AS4", #waterfall is close to AS4
                                           "AS1", "AS2", "AS3", "AS4", "AS5", 
                                           "NA",
                                           "BCI-NE", "BCI-NE", "BCI-NE", "BCI-NE", "BCI-NE", "BCI-NE", "BCI-NE", #BCI #BCI_Brbour_lathrop
                                           "BCI-NE", "BCI-NE", #BCI-SW BCIA
                                           "BCI-NE", 
                                           "Bohio", "Bohio", 
                                           "Chicha", 
                                           "Frijoles",
                                           "Gigante", "Gigante", 
                                           "Guanabano", "Guava", "NA", "Mona Grita", "Mona Grita", "Mona Grita", "Orchid", "Pato Horqueta", 
                                           "Pena Blanca", "Pato Horqueta", "Tres Almendras", "Tres Almendras", "Tres Almendras", "NA")
capture_dat<-subset(capture_dat, location!="NA")
capture_dat$landscape<-as.factor(capture_dat$landscape)
levels(capture_dat$landscape) <- c("A", "A", 
                                  "C", "C", "C", "C", "C", 
                                  "I", "I", "I", "I", "C", 
                                  "I", "I", "C", "I", "I")

                    
library(dplyr)
library(tidyr)

# Step 0: Clean astro column
capture_dat <- capture_dat %>%
  mutate(astro = toupper(astro))  # fix NEg → NEG

# Step 1: Count POS and NEG per species × location × landscape
location_counts <- capture_dat %>%
  group_by(species, landscape, location, astro) %>%
  summarise(n = n(), .groups = "drop") %>%
  pivot_wider(names_from = astro, values_from = n, values_fill = 0)

# Step 2: Calculate prevalence per location
location_counts <- location_counts %>%
  mutate(prevalence = POS / (POS + NEG))

# Step 3: Aggregate per landscape: mean ± SE of location-level prevalence
landscape_summary <- location_counts %>%
  group_by(species, landscape) %>%
  summarise(
    mean_prev = mean(prevalence, na.rm = TRUE),
    se_prev   = sd(prevalence, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  )

landscape_summary


# Add location counts per species × landscape for sample size annotation
landscape_summary <- location_counts %>%
  group_by(species, landscape) %>%
  summarise(
    mean_prev    = mean(prevalence, na.rm = TRUE),
    se_prev      = sd(prevalence, na.rm = TRUE) / sqrt(n()),
    n_locations  = n(),                                 # number of locations
    n_bats       = sum(POS + NEG, na.rm = TRUE),       # total number of bats
    .groups      = "drop"
  )




# Prepare data: location-level prevalence
location_prev <- location_counts %>%
  select(species, landscape, location, prevalence)

location_prev_AJ<-subset(location_prev, species=="A.jamaicensis")
location_prev<-subset(location_prev, landscape!="I")
location_prev<-subset(location_prev, species=="C.perspicillata")

# One-sided t-test, raw prevalence
t.test(prevalence ~ landscape,
       data = location_prev,
       alternative = "greater") 

anova(lm(prevalence ~ landscape, data=location_prev_AJ))

##### some basis for just looking at carollia but more generally the prevalence is just higher in agricultural sites ####

levels(landscape_summary$landscape)<-c("Forest Fragments", "Continuous Forest", "Forested Islands")

# Reorder factor levels BEFORE plotting
landscape_summary$landscape <- factor(
  landscape_summary$landscape,
  levels = c("Continuous Forest", "Forest Fragments", "Forested Islands")
)

# Add annotation to plot
astro_prev_plot <- ggplot(subset(landscape_summary, landscape != "Forested Islands"), 
                          aes(x = (landscape), y = mean_prev, fill=landscape)) +
  geom_col(position = position_dodge(0.8), width = 0.7, colour = "black") +
  geom_errorbar(aes(ymin = mean_prev - se_prev, ymax = mean_prev + se_prev),
                position = position_dodge(0.8), width = 0.2) +
  geom_text(aes(label = paste0("nL=", n_locations, "\n", "nB=", n_bats)),
            position = position_dodge(0.8), vjust = 3, size = 4, color = "black") +
  labs(x = "Landscape", y = "Mean AstV prevalence", fill = "Landscape; p-value = 0.002") +
  scale_fill_manual(values = c("#0E7860", "#F3D112")) +
  scale_x_discrete(labels = c("Forest Fragments" = "Forest Fragments", "Continuous Forest" = "Continuous Forest")) +
  ylim(0,.3)+
  theme_classic() +
  theme(
    axis.title = element_text(face = "bold", size = 14),
    legend.text = element_text(face = "italic", size = 12),
    legend.title = element_text(face = "bold", size = 12),
    legend.position = c(.3, .9)
  )

astro_prev_plot



#################################################################
########              microbiome data                    ########
#################################################################

library(phyloseq)
ps_unrarefied <- readRDS("~/Dropbox/PanamaBats/workingdata/ps_131025.rds")
ps_unrarefied <- subset_samples(ps_unrarefied, species=="C.perspicillata")
ps_unrarefied <- subset_samples(ps_unrarefied, astro!="NA" & astro!="")
ps_unrarefied <- subset_samples(ps_unrarefied, landscape!="I")
ps_unrarefied<-subset_samples(ps_unrarefied, age!="NA")
ps_unrarefied<-subset_samples(ps_unrarefied, sex!="NA")

microbiome::summarize_phyloseq(ps_unrarefied)

###### adding season ####### 
library(lubridate)

sample_data(ps_unrarefied)$date<-as.Date(sample_data(ps_unrarefied)$date, format = "%d.%m.%y")

# make sure it's Date
sample_data(ps_unrarefied)$date <- as.Date(sample_data(ps_unrarefied)$date)

# extract day-of-year (1–365/366)
doy <- yday(sample_data(ps_unrarefied)$date)

# wet season: April 15 (~day 106) – Nov 15 (~day 320), else dry
sample_data(ps_unrarefied)$season <- ifelse(doy >= 106 & doy <= 320, "wet", "dry")

######
meta<-as(sample_data(ps_unrarefied), "data.frame")
ps_rarefied <- rarefy_even_depth(ps_unrarefied, sample.size = 5000, rngseed = 123, verbose = TRUE)
meta_rarefied<-as(sample_data(ps_rarefied), "data.frame")

sample_sizes_cp <- plyr::ddply(meta_rarefied, c("landscape", "astro"), summarise, n = length(landscape))
sample_sizes_unrarefied_cp <- plyr::ddply(meta, c("landscape", "astro"), summarise, n = length(landscape))

library(ggh4x)

#### sample size #####
# set landscape order: C first, then A
sample_sizes_cp$landscape <- factor(
  sample_sizes_cp$landscape,
  levels = c("C", "A")
)

sample_size_plot_cp <- ggplot(sample_sizes_cp, aes(landscape, n, fill = astro)) +
  geom_bar(stat = "identity", colour = "black") +
  scale_fill_manual(values = c("#3177B5", "#D14538")) +
  theme_classic() +
  theme(
    axis.title = element_text(face = "bold", size = 14),
    legend.title = element_text(face = "bold", size = 12),
    legend.text = element_text(size = 10),
    legend.position = c(0.3, .9)
  ) +
  labs(
    x = "Landscape",
    y = "Sample Size",
    fill = "Astrovirus"
  ) +
  scale_x_discrete(
    labels = c(
      "A" = "Forest Fragments",
      "C" = "Continuous Forest"
    )
  )


ggpubr::ggarrange(
  astro_prev_plot,
  sample_size_plot_cp,
  labels = c("B", "C"),
  widths = c(2, 1)
)


#################################################################
########      compare rarefied vs unrarefied reads       ########
#################################################################

# --- unrarefied alpha diversity ---
alpha_unrarefied <- estimate_richness(ps_unrarefied) 
rownames(alpha_unrarefied) <- sub("^X", "", rownames(alpha_unrarefied))
rownames(alpha_unrarefied) <- gsub("\\.", "-", rownames(alpha_unrarefied))
df.pd_unrare <- picante::pd(t(as.data.frame(otu_table(ps_unrarefied))),
                            phy_tree(ps_unrarefied), include.root=TRUE)
meta$Shannon_unrare <- alpha_unrarefied$Shannon 
meta$Observed_unrare <- alpha_unrarefied$Observed
meta$Faiths_unrare <- df.pd_unrare$PD

# --- rarefied alpha diversity ---
alpha_rarefied <- estimate_richness(ps_rarefied) 
rownames(alpha_rarefied) <- sub("^X", "", rownames(alpha_rarefied))
rownames(alpha_rarefied) <- gsub("\\.", "-", rownames(alpha_rarefied))

df.pd_rare <- picante::pd(t(as.data.frame(otu_table(ps_rarefied))),
                          phy_tree(ps_rarefied), include.root=TRUE)
rownames(df.pd_rare) <- sub("^X", "", rownames(df.pd_rare))
rownames(df.pd_rare) <- gsub("\\.", "-", rownames(df.pd_rare))

# --- build meta_rarefied ---
meta_rarefied <- as(sample_data(ps_rarefied), "data.frame")

# Make sure everything has a SampleID column
alpha_rarefied$SampleID <- rownames(alpha_rarefied)
alpha_unrarefied$SampleID <- rownames(alpha_unrarefied)
df.pd_rare$SampleID <- rownames(df.pd_rare)
df.pd_unrare$SampleID <- rownames(df.pd_unrare)
meta_rarefied$SampleID <- rownames(meta_rarefied)

# Merge rarefied first
meta_rarefied <- meta_rarefied %>%
  left_join(alpha_rarefied %>% 
              select(SampleID,
                     Shannon_rare = Shannon,
                     Observed_rare = Observed),
            by = "SampleID") %>%
  left_join(df.pd_rare %>% 
              select(SampleID,
                     Faiths_rare = PD),
            by = "SampleID")

# Add unrarefied metrics for the same samples
meta_rarefied <- meta_rarefied %>%
  left_join(alpha_unrarefied %>% 
              select(SampleID,
                     Shannon_unrare = Shannon,
                     Observed_unrare = Observed),
            by = "SampleID") %>%
  left_join(df.pd_unrare %>% 
              select(SampleID,
                     Faiths_unrare = PD),
            by = "SampleID")



# --- 3. Correlations between rarefied & unrarefied ---
make_corrplot <- function(df, x, y, xlab, ylab) {
  test <- cor.test(df[[x]], df[[y]])
  label_text <- sprintf("R² = %.3f, p = %.2g", test$estimate, test$p.value)
  
  ggplot(df, aes_string(x, y)) +
    geom_point(alpha = 0.7) +
    geom_smooth(method = "lm", se = FALSE, color = "black") +
    theme_bw(base_size = 12) +
    labs(x = xlab, y = ylab) +
    annotate("text", x = -Inf, y = Inf, hjust = -0.1, vjust = 1.1,
             label = label_text, size = 4)
}

Obs_unrareVSrare <- make_corrplot(meta_rarefied, "Observed_rare", "Observed_unrare",
                                  "Rarefied Observed ASVs", "Unrarefied Observed ASVs")

Sha_unrareVSrare <- make_corrplot(meta_rarefied, "Shannon_rare", "Shannon_unrare",
                                  "Rarefied Shannon Index", "Unrarefied Shannon Index")

Pd_unrareVSrare <- make_corrplot(meta_rarefied, "Faiths_rare", "Faiths_unrare",
                                 "Rarefied Faith's PD", "Unrarefied Faith's PD")

# --- 4. Correlations with sequencing depth ---
make_seqplot <- function(df, metric_rare, metric_unrare, ylab) {
  plot_data <- df %>%
    select(seq_depth, !!metric_rare, !!metric_unrare) %>%
    pivot_longer(cols = c(!!metric_rare, !!metric_unrare),
                 names_to = "type", values_to = "value")
  
  # Run correlations separately
  test_rare   <- cor.test(df[[metric_rare]],   df$seq_depth)
  test_unrare <- cor.test(df[[metric_unrare]], df$seq_depth)
  label_rare   <- sprintf("R² = %.3f, p = %.2g", test_rare$estimate, test_rare$p.value)
  label_unrare <- sprintf("R² = %.3f, p = %.2g", test_unrare$estimate, test_unrare$p.value)
  
  ggplot(plot_data, aes(x = seq_depth, y = value, color = type)) +
    geom_point(alpha = 0.6) +
    geom_smooth(method = "lm", se = FALSE) +
    theme_bw(base_size = 12) +
    labs(x = "Sequencing depth", y = ylab) +
    scale_color_manual(values = c("Shannon_unrare" = "darkgrey", "Shannon_rare" = "lightgrey",
                                  "Observed_unrare" = "darkgrey", "Observed_rare" = "lightgrey",
                                  "Faiths_unrare" = "darkgrey", "Faiths_rare" = "lightgrey"),
                       labels = c("Rarefied", "Unrarefied")) +
    annotate("text", x = Inf, y = -Inf, hjust = 1.1, vjust = -0.5,
             label = label_rare, color = "lightgrey", size = 4) +
    annotate("text", x = Inf, y = -Inf, hjust = 1.1, vjust = -2,
             label = label_unrare, color = "darkgrey", size = 4)
}

ObsVSseq <- make_seqplot(meta_rarefied, "Observed_rare", "Observed_unrare", "Observed ASVs")
ShaVSseq <- make_seqplot(meta_rarefied, "Shannon_rare", "Shannon_unrare", "Shannon Index")
PDVSseq <- make_seqplot(meta_rarefied, "Faiths_rare", "Faiths_unrare", "Faith's PD Index")

# --- 5. Arrange plots in a publication-ready panel ---
ggpubr::ggarrange(Obs_unrareVSrare, Sha_unrareVSrare, Pd_unrareVSrare, ObsVSseq, ShaVSseq, PDVSseq,
                  labels = c("A", "B", "C", "D", "E", "F"), common.legend = TRUE,
                  legend = "bottom", align = "v")





#################################################################
########      compositional plots                        ########
#################################################################


library(phyloseq)
library(microbiome)
library(ggplot2)
library(dplyr)
library(forcats)

#----------------------------
# Subset to Carollia perspicillata only
#----------------------------
ps_CP <- subset_samples(ps_rarefied, species == "C.perspicillata")

#----------------------------
# Collapse to Genus level
#----------------------------
ps_genus <- tax_glom(ps_CP, taxrank = "Genus")

# Convert to relative abundance (compositional)
ps_genus <- microbiome::transform(ps_genus, "compositional")

#----------------------------
# Melt for ggplot
#----------------------------
genus_per_sample <- psmelt(ps_genus)

#----------------------------
# Identify dominant genera (as in your original list)
#----------------------------
top_genus <- c(
  "Enterobacteriaceae Family", "Streptococcus", "Escherichia-Shigella",
  "Ureaplasma", "Mannheimia", "Staphylococcus", 
  "Mycoplasma", "Pseudomonas", "Clostridium_sensu_stricto_1",
  "Pantoea", "Acinetobacter", "Haemophilus", "Stenotrophomonas",
  "Sphingomonas", "Methylobacterium-Methylorubrum"
)

genus_colours <- c(
  "Enterobacteriaceae Family" = "#D55E00",
  "Streptococcus" = "#E69F00",
  "Escherichia-Shigella" = "#56B4E9",
  "Ureaplasma" = "#01AE80",
  "Mannheimia" = "#CC79A7",
  "Staphylococcus" = "#F5E834",
  "Mycoplasma" = "#0072B2",
  "Pseudomonas" = "#999999",
  "Clostridium_sensu_stricto_1" = "#E41A1C",
  "Pantoea" = "#7CA9F1",
  "Acinetobacter" = "#48AB44",
  "Haemophilus" = "#984EA3",
  "Stenotrophomonas" = "#F5EE9E",
  "Sphingomonas" = "#F781BF",
  "Methylobacterium-Methylorubrum" = "#8781BF",
  "Other" = "grey80"
)

#----------------------------
# Prepare dataframe for plotting
#----------------------------
genus_per_sample <- genus_per_sample %>%
  mutate(
    Genus_plot = ifelse(Genus %in% top_genus, Genus, "Other"),
    Genus_plot = factor(Genus_plot, levels = c(top_genus, "Other")),
    landscape = fct_recode(
      landscape,
      "Continuous forest" = "C",
      "Forest fragments" = "A"
    ),
    # enforce order: Continuous first
    landscape = factor(landscape, levels = c("Continuous forest", "Forest fragments"))
  )


#----------------------------
# Final plot (Carollia only)
#----------------------------
Carollia_com_plot <- ggplot(genus_per_sample, 
                            aes(x = astro, y = Abundance, fill = Genus_plot)) +
  geom_col(position = "fill") +
  facet_wrap(~ landscape) +
  scale_fill_manual(values = genus_colours, drop = FALSE) +
  labs(
    y = "Relative Abundance",
    x = "Astrovirus Infection",
    fill = "Bacterial Genus"
  ) +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(face = "italic", size = 14, hjust = 0.5),
    legend.position = "right",
    legend.title = element_text(face = "bold"),
    axis.title.y = element_text(face = "bold"),
    axis.title.x = element_text(face = "bold"),
    strip.text = element_text(face = "bold", size = 11)
  )

Carollia_com_plot

##### family stats #####

# Collapse to Family level
ps_family <- tax_glom(ps_CP, taxrank = "Family")

# Convert to compositional abundance
ps_family <- microbiome::transform(ps_family, "compositional")

# Melt for inspection
family_per_sample <- psmelt(ps_family)

# Print ranked family abundances
print(
  microbiomeutilities::get_group_abundances(
    ps_family, 
    level = "Family", 
    group = "Family", 
    transform = "compositional"
  ) %>% arrange(desc(mean_abundance)),
  n = 35
)


#######################################################
##########        alpha diversity          ############
#######################################################

meta_CP<-subset(meta_rarefied, species=="C.perspicillata")

# --- Packages ---
library(lmerTest)
library(MuMIn)
library(DHARMa)
library(car)

######## shannon #########
# --- Fit full model ---
lm_full <- lmer(Shannon_rare ~ astro * landscape + 
                  sex + astro * age + season + seq_depth + 
                  (1 | location),
                data = meta_CP, na.action = na.fail)

# --- Check multicollinearity ---
vif_vals <- car::vif(lm_full)
print(vif_vals)

# --- Backward stepwise selection (manual or using step()) ---
# step() is not directly compatible with lmerTest, but drop1() gives the same logic
drop1(lm_full, test = "Chisq")

# Suppose astro:age is nonsignificant
lm_step1 <- update(lm_full, . ~ . - astro:age)
drop1(lm_step1, test = "Chisq")

# Next, if astro:landscape is nonsignificant
#lm_step2 <- update(lm_step1, . ~ . - astro:landscape)
#drop1(lm_step2, test = "Chisq")

# You can iteratively remove nonsignificant terms:
# e.g. lm_reduced <- update(lm_full, . ~ . - astro:landscape:species)

lm_final <- lm_step1  # last model after removing interactions
drop1(lm_final, test = "Chisq")
summary(lm_final)
car::vif(lm_final)

# --- Check model assumptions ---
res_full <- simulateResiduals(fittedModel = lm_final, n = 1000)
plot(res_full)


######## observed #####
# --- Fit full model ---
lm_full <- lmer(log(Observed_rare) ~ astro * landscape + 
                  sex + astro * age + season + seq_depth + 
                  (1 | location),
                data = meta_CP, na.action = na.fail)

# --- Check multicollinearity ---
vif_vals <- car::vif(lm_full)
print(vif_vals)

# --- Backward stepwise selection (manual or using step()) ---
# step() is not directly compatible with lmerTest, but drop1() gives the same logic
drop1(lm_full, test = "Chisq")

# Suppose astro:age is nonsignificant
lm_step1 <- update(lm_full, . ~ . - astro:age)
drop1(lm_step1, test = "Chisq")

# Next, if astro:landscape is nonsignificant
lm_step2 <- update(lm_step1, . ~ . - astro:landscape)
drop1(lm_step2, test = "Chisq")

# You can iteratively remove nonsignificant terms:
# e.g. lm_reduced <- update(lm_full, . ~ . - astro:landscape:species)

lm_final <- lm_step2  # last model after removing interactions
drop1(lm_final, test = "Chisq")
summary(lm_final)
car::vif(lm_final)

# --- Check model assumptions ---
res_full <- simulateResiduals(fittedModel = lm_final, n = 1000)
plot(res_full)



######## Faiths #####
# --- Fit full model ---
lm_full <- lmer(log(Faiths_rare) ~ astro * landscape + 
                  sex + astro * age + season + seq_depth + 
                  (1 | location),
                data = meta_CP, na.action = na.fail)

# --- Check multicollinearity ---
vif_vals <- car::vif(lm_full)
print(vif_vals)

# --- Backward stepwise selection (manual or using step()) ---
# step() is not directly compatible with lmerTest, but drop1() gives the same logic
drop1(lm_full, test = "Chisq")

# Suppose astro:age is nonsignificant
lm_step1 <- update(lm_full, . ~ . - astro:age)
drop1(lm_step1, test = "Chisq")

# Next, if astro:landscape is nonsignificant
lm_step2 <- update(lm_step1, . ~ . - astro:landscape)
drop1(lm_step2, test = "Chisq")

# You can iteratively remove nonsignificant terms:
# e.g. lm_reduced <- update(lm_full, . ~ . - astro:landscape:species)

lm_final <- lm_step2  # last model after removing interactions
drop1(lm_final, test = "Chisq")
summary(lm_final)
car::vif(lm_final)

# --- Check model assumptions ---
res_full <- simulateResiduals(fittedModel = lm_final, n = 1000)
plot(res_full)



####### visualisation #########

# Create a grouping variable combining landscape and astro
meta_CP$group <- paste(meta_CP$landscape, meta_CP$astro, sep = "_")

# Reorder so landscape C groups appear first, landscape A groups second
meta_CP$group <- factor(
  meta_CP$group,
  levels = c(
    "C_NEG", "C_POS",
    "A_NEG", "A_POS"
  )
)

library(beanplot)
library(ggplotify)   # for as.ggplot()
library(ggpubr)

# 1. Create the beanplot inside a function (required)
bean_fun <- function() {
  beanplot(
    Shannon_rare ~ group,
    data = meta_CP,
    col = list("#3177B5", "#D14538", "#3177B5", "#D14538"),
    what = c(1,1,1,1),
    border = "black",
    side = "both",
    overallline = "median",
    xlab = "Landscape",
    ylab = "Shannon Diversity",
    cex.lab = 1.4,
    font.lab = 3,
    xaxt = "n"
  )
  axis(
    1,
    at = 1:2,
    labels = c("Continuous Forest", "Forest Fragments"),
    cex.axis = 1.2
  )
}

# 2. Convert base plot → ggplot-like object
bean_grob <- as.ggplot(bean_fun)

# 3. Combine with ggplot objects
ggarrange(Carollia_com_plot,
  bean_grob,
  labels=c("A","B"),
  ncol = 2
)


#####################################################

###### beta diversity ######
# ============================================================
# 1. Setup
# ============================================================
library(phyloseq)
library(vegan)
library(compositions)
library(tidyverse)

# Focus species
ps_carollia <- subset_samples(ps_rarefied, species == "C.perspicillata")

# Metadata
meta <- as(sample_data(ps_carollia), "data.frame")

# ============================================================
# 2. Prepare taxonomic levels
# ============================================================
ps_list <- list(
  ASV     = ps_carollia,
  Genus   = ps_carollia %>%
    tax_glom("Genus") %>% subset_taxa(!is.na(Genus) & !Genus %in% c("uncultured"))
  )

# ============================================================
# 3. Distance function
# ============================================================
get_distances <- function(ps_obj) {
  otu <- as.data.frame(t(otu_table(ps_obj)))
  
  list(
    bray      = phyloseq::distance(ps_obj, method = "bray"),
    jaccard   = phyloseq::distance(ps_obj, method = "jaccard"),
    wunifrac  = phyloseq::distance(ps_obj, method = "wunifrac"),
    unifrac   = phyloseq::distance(ps_obj, method = "unifrac"),
    aitchison = {
      clr_mat <- as.matrix(clr(otu + 1))
      dist(clr_mat, method = "euclidean")
    }
  )
}

# ============================================================
# 4. Model runner for adonis2
# ============================================================
run_adonis_models <- function(dist_obj, meta) {
  # Ensure metadata alignment
  sample_ids <- rownames(as.matrix(dist_obj))
  meta_sub <- meta[meta$SampleID %in% sample_ids, ]
  meta_sub <- meta_sub[match(sample_ids, meta_sub$SampleID), ]
  if (nrow(meta_sub) != length(sample_ids)) stop("Sample mismatch.")
  
  # Full model (with interactions)
  mod_full <- adonis2(
    dist_obj ~ astro + landscape*astro + sex+astro + age*astro + season + seq_depth,
    data = meta_sub, permutations = 999, strata = meta_sub$location, by = "margin"
  )
  
  mod_reduced <- adonis2(
    dist_obj ~ landscape*astro + sex + age+astro + season + seq_depth,
    data = meta_sub, permutations = 999, strata = meta_sub$location, by = "margin"
  )

  # Main-effects-only model
  mod_simple <- adonis2(
    dist_obj ~ astro + landscape + sex + age + season + seq_depth,
    data = meta_sub, permutations = 999, strata = meta_sub$location, by = "margin"
  )
  
  list(full = mod_full, reduced = mod_reduced, main = mod_simple)
}

# ============================================================
# 5. Run all analyses and collect summary table
# ============================================================
summary_table <- list()

for (lvl in names(ps_list)) {
  cat("Processing:", lvl, "\n")
  ps_sub <- ps_list[[lvl]]
  dists <- get_distances(ps_sub)
  
  for (dname in names(dists)) {
    dist_obj <- dists[[dname]]
    res <- run_adonis_models(dist_obj, meta)
    
    # Extract clean results
    df_full <- as.data.frame(res$full)[, c("Df","SumOfSqs","R2","F","Pr(>F)")]
    df_full$model <- "Full"
    df_full$Distance <- dname
    df_full$Level <- lvl
    
    df_reduced <- as.data.frame(res$reduced)[, c("Df","SumOfSqs","R2","F","Pr(>F)")]
    df_reduced$model <- "Full"
    df_reduced$Distance <- dname
    df_reduced$Level <- lvl
    
    df_main <- as.data.frame(res$main)[, c("Df","SumOfSqs","R2","F","Pr(>F)")]
    df_main$model <- "Main"
    df_main$Distance <- dname
    df_main$Level <- lvl
    
    summary_table[[paste(lvl, dname, sep = "_")]] <- bind_rows(df_full, df_reduced, df_main)
  }
}

summary_df <- bind_rows(summary_table)

# ============================================================
# 6. Export results
# ============================================================
write.csv(summary_df, "~/Dropbox/PanamaBats/Stefans_manuscript/Carollia_adonis_summary.csv", row.names = FALSE)
print(summary_df)





library(phyloseq)
library(dplyr)
library(ggplot2)
library(patchwork)
library(vegan)
library(compositions)

# ---------------------------
# Aitchison distance function
# ---------------------------
get_aitchison_dist <- function(ps_obj) {
  otu <- as.data.frame(t(otu_table(ps_obj)))
  clr_mat <- as.matrix(compositions::clr(otu + 1))
  dist(clr_mat, method = "euclidean")
}

# ---------------------------
# PCoA plotting function
# ---------------------------
plot_pcoa_aitchison <- function(ps_obj, level_label, species_name, banner_color) {
  # Aggregate at desired taxonomic level
  ps_level <- if (tolower(level_label) == "genus") {
    tax_glom(ps_obj, taxrank = "Genus")
  } else if (tolower(level_label) == "asv") {
    ps_obj
  } else stop("level_label must be either 'ASV' or 'Genus'")
  
  # Compute distance and ordination
  aitch <- get_aitchison_dist(ps_level)
  ord <- ordinate(ps_level, method = "PCoA", distance = aitch)
  
  # Extract ordination scores
  df <- data.frame(ord$vectors[, 1:2])
  df$SampleID <- rownames(df)
  meta <- as(sample_data(ps_level), "data.frame")
  df <- left_join(df, meta, by = "SampleID")
  
  # Variance explained
  eigvals <- ord$values$Eigenvalues
  var_expl <- round(100 * eigvals / sum(eigvals), 1)
  
  # Compute centroids by astro
  centroids <- df %>%
    group_by(astro) %>%
    summarise(
      Axis.1 = mean(Axis.1, na.rm = TRUE),
      Axis.2 = mean(Axis.2, na.rm = TRUE)
    )
  
  # Main PCoA plot
  p_main <- ggplot(df, aes(x = Axis.1, y = Axis.2, fill = astro)) +
    geom_point(size = 2, alpha = 0.8, shape = 21, colour = "black") +
    stat_ellipse(aes(fill = astro), geom = "polygon", alpha = 0.2, colour = NA) +
    geom_point(
      data = centroids,
      aes(x = Axis.1, y = Axis.2, fill = astro),
      shape = 23, size = 6,
      colour = "white", stroke = 1.5
    ) +
    scale_fill_manual(values = c("NEG" = "#3177B5", "POS" = "#D14538")) +
    labs(
      x = paste0("PCoA1 (", var_expl[1], "%)"),
      y = paste0("PCoA2 (", var_expl[2], "%)")
    ) +
    theme_bw(base_size = 14) +
    theme(
      panel.grid = element_blank(),
      legend.position = "none",
      legend.title = element_blank(), 
      axis.title = element_text(face = "bold", size = 14)
    )
  
}

# ---------------------------
# Apply PCoA for Carollia only
# ---------------------------
p_CP_genus <- plot_pcoa_aitchison(
  ps_obj = ps_CP,
  level_label = "genus",
  species_name = "Carollia perspicillata",
  banner_color = "#F3721C"
)

p_CP_asv <- plot_pcoa_aitchison(
  ps_obj = ps_CP,
  level_label = "ASV",
  species_name = "Carollia perspicillata",
  banner_color = "#F3721C"
)


# ---------------------------
# Beta-dispersion analysis
# ---------------------------
# Compute Aitchison distance
dist_aitch <- get_aitchison_dist(ps_CP)
meta_CP <- as(sample_data(ps_CP), "data.frame")

# Beta-dispersion for landscape
bd_landscape <- betadisper(dist_aitch, meta_CP$landscape)
anova(bd_landscape)          # test if dispersion differs between landscapes
permutest(bd_landscape)      # permutation test
plot(bd_landscape)            # visualize spread by landscape

# Extract distances and metadata
bd_landscape_df <- data.frame(
  distance = bd_landscape$distances,
  landscape = meta_CP$landscape
)

# Reorder factor levels so Continuous Forest comes first
bd_landscape_df$landscape <- factor(
  bd_landscape_df$landscape,
  levels = c("C", "A")  # C = Continuous Forest, A = Forest Fragments
)

# Boxplot
p_bd_landscape <- ggplot(bd_landscape_df, aes(x = landscape, y = distance, fill = landscape)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.1) +
  scale_fill_manual(
    values = c("C" = "#0E7860", "A" = "#F3D112"),
    labels = c("Continuous Forest", "Forest Fragments")
  ) +
  scale_x_discrete(labels = c("A" = "Forest Fragments", "C" = "Continuous Forest"))+
  labs(x = "Landscape", y = "Distance to Centroid") +
  theme_classic(base_size = 14) +
  theme(
    axis.title = element_text(face = "bold", size = 14),
    legend.position = "none"
  )

# Beta-dispersion for astro
bd_astro <- betadisper(dist_aitch, meta_CP$astro)
anova(bd_astro)
permutest(bd_astro)
plot(bd_astro)

# Extract distances and metadata
bd_astro_df <- data.frame(
  distance = bd_astro$distances,
  astro = meta_CP$astro
)

# Boxplot
p_bd_astro <- ggplot(bd_astro_df, aes(x = astro, y = distance, fill = astro)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.1) +
  scale_fill_manual(values = c("NEG" = "#3177B5", "POS" = "#D14538")) +
  labs(x = "Astrovirus infection", y = "Distance to Centroid") +
  theme_classic(base_size = 14) +
  theme(
    axis.title = element_text(face = "bold", size = 14),
    legend.position = "none"
  )

# ---- Combine with patchwork ----
ggpubr::ggarrange(p_CP_genus, p_bd_astro, p_bd_landscape, labels=c("A", "B", "C"), ncol=3,widths=c(2,1,1))


#################
#### DA #########
#################

##### because the interaction with landscape is rarely if ever significant we only would look at AstV

library(phyloseq)
library(ANCOMBC)
library(dplyr)
library(tibble)

ps_CP <- subset_samples(ps_rarefied, species == "C.perspicillata")
ps_CP_10 <- subset_taxa(ps_CP, rowSums(otu_table(ps_CP) > 0) > (0.1 * nsamples(ps_CP)))

ancom_CP_asv <- ancombc2(
  data = ps_CP_10,
  assay_name = "counts",           # default OTU table slot name
  tax_level = NULL,                # NULL = ASV level
  fix_formula = "astro + landscape",
  rand_formula = NULL,             # no random effects here
  p_adj_method = "holm",
  struc_zero = TRUE,               # structural zeros (recommended)
  neg_lb = TRUE,                   # more conservative
  global = TRUE,                   # global tests for each factor
  alpha = 0.05,
  lib_cut = 0,                     # keep all samples
  group = "astro",                 # main grouping factor (for global tests)
  verbose = TRUE
)

res_asv <- ancom_CP_asv$res
sig_asv <- res_asv %>%
  filter(q_astroPOS < 0.05 | q_landscapeC < 0.05 ) %>%
  arrange(q_astroPOS)

head(sig_asv)

# Extract taxonomy table from your phyloseq object
tax_df <- as.data.frame(tax_table(ps_CP))  # or ps_ASV, whichever matches your sig_asv
tax_df$taxon <- rownames(tax_df)

# Check what ranks are available
head(tax_df)
# e.g., Kingdom, Phylum, Class, Order, Family, Genus, Species, taxon

# Join taxonomy info with sig_asv
sig_asv_annotated <- sig_asv %>%
  left_join(tax_df %>% select(taxon, Genus, Species), by = "taxon")

# Inspect result
head(sig_asv_annotated)



#### at genus level 
ps_genus <- tax_glom(ps_CP, taxrank = "Genus")
ps_genus_10 <- subset_taxa(ps_genus, rowSums(otu_table(ps_genus) > 0) > (0.1 * nsamples(ps_genus)))

ancom_CP_genus <- ancombc2(
  data = ps_genus_10,
  tax_level = "Genus",
  fix_formula = "astro + landscape",
  rand_formula = NULL,
  p_adj_method = "holm",
  struc_zero = TRUE,
  neg_lb = TRUE,
  global = TRUE,
  alpha = 0.05,
  lib_cut = 0,
  group = "astro",
  verbose = TRUE
)

res_genus <- ancom_CP_genus$res
sig_genus <- res_genus %>%
  filter(q_astroPOS < 0.05 | q_landscapeC < 0.05 ) %>%
  arrange(q_astroPOS)

head(sig_genus)

sig_plot <- res_genus %>%
  filter( q_astroPOS < 0.05) %>%
  mutate(
    direction = ifelse(lfc_astroPOS > 0, "Infected ↑", "Infected ↓")
  )

# --- Save ASV-level results ---
write.csv(sig_asv_annotated, "~/Dropbox/PanamaBats/Stefans_manuscript/ANCOMBC2_results_ASV_all.csv", row.names = FALSE)

# --- Save Genus-level results ---
write.csv(sig_genus, "~/Dropbox/PanamaBats/Stefans_manuscript/ANCOMBC2_results_Genus_all.csv", row.names = FALSE)


ggplot(sig_plot, aes(x = reorder(taxon, lfc_astroPOS), y = lfc_astroPOS, fill = direction)) +
  geom_col() +
  coord_flip() +
  theme_minimal(base_size = 12) +
  scale_fill_manual(values = c("Infected ↑" = "#D14538", "Infected ↓" = "#3177B5")) +
  labs(
    x = NULL, y = "Log Fold Change (Astro+ vs Astro–)",
    title = expression(italic("Carollia perspicillata")),
    subtitle = "Differentially abundant genera (ANCOM-BC2)"
  )






library(phyloseq)
library(tidyverse)

# Collapse to genus level
ps_genus <- tax_glom(ps_CP, taxrank = "Genus")

# Filter: keep genera present in >10% of samples
ps_genus_10 <- subset_taxa(ps_genus, rowSums(otu_table(ps_genus) > 0) > (0.1 * nsamples(ps_genus)))

# Melt to long format
genus_long <- psmelt(ps_genus_10)

# Combine infection and landscape into one factor
genus_long <- genus_long %>%
  mutate(group = paste(astro, landscape, sep = "_")) %>%
  group_by(Genus, group) %>%
  summarise(mean_abundance = mean(Abundance), .groups = "drop") %>%
  pivot_wider(names_from = group, values_from = mean_abundance, values_fill = 0)

# Check what columns we actually have
print(colnames(genus_long))

# Rename consistently to "Fragmented_AstroPOS", etc.
genus_long <- genus_long %>%
  rename(
    Fragmented_AstroPOS = POS_A,
    Fragmented_AstroNEG = NEG_A,
    Continuous_AstroPOS = POS_C,
    Continuous_AstroNEG = NEG_C
  )

# Normalize rows to proportions
genus_long <- genus_long %>%
  mutate(total = rowSums(across(where(is.numeric)))) %>%
  mutate(across(where(is.numeric), ~ .x / total)) %>%
  select(-total)


sig_genus <- res_genus %>%
  filter(lfc_astroPOS > 1.5 | lfc_astroPOS  < -1.5) %>%
  arrange(lfc_astroPOS)

sig_genus_frag <- res_genus %>%
  filter(lfc_landscapeC > 1.5 | lfc_landscapeC  < -1.5) %>%
  arrange(lfc_landscapeC)

sig_genus<-rbind(sig_genus_frag, sig_genus)

# Ensure sig_genus has a column "taxon" matching genus_long$Genus
sig_genera <- sig_genus$taxon

# Compute ternary coordinates and change magnitude
ternary_data <- genus_long %>%
  rowwise() %>%
  mutate(
    reference = Continuous_AstroNEG,
    astro_change = Continuous_AstroPOS + Fragmented_AstroPOS,
    fragment_change = Fragmented_AstroNEG + Fragmented_AstroPOS,
    total = reference + astro_change + fragment_change,
    # normalize for ternary plot
    reference = reference / total,
    astro_change = astro_change / total,
    fragment_change = fragment_change / total,
    change_magnitude = sqrt((astro_change)^2 + (fragment_change)^2) # Euclidean-like measure
  ) %>%
  ungroup() %>%
  mutate(
    label = ifelse(Genus %in% sig_genera, Genus, "")
  )

# Plot ternary
tenaryCP<-plot_ly(
  ternary_data,
  type = 'scatterternary',
  mode = 'markers+text',
  a = ~reference,
  b = ~astro_change,
  c = ~fragment_change,
  text = ~label,
  textposition = 'top center',
  marker = list(
    size = ~change_magnitude * 20,  # scale factor for visibility
    color = "#F3721C",
    line = list(width = 1, color = 'black')
  )
) %>%
  layout(
    title = "Ternary plot: Change relative to Continuous Astro–",
    ternary = list(
      sum = 1,
      aaxis = list(title = "Continuous Astro– (reference)"),
      baxis = list(title = "Astrovirus positive"),
      caxis = list(title = "Fragmented Forest")
    )
  )


library(phyloseq)
library(tidyverse)

# Collapse to genus level
ps_genus <- tax_glom(ps_CP, taxrank = "Genus")

# Filter: keep genera present in >10% of samples
ps_genus_10 <- subset_taxa(ps_genus, rowSums(otu_table(ps_genus) > 0) > (0.1 * nsamples(ps_genus)))

# Mean relative abundance per genus across all samples
genus_abund <- psmelt(ps_genus_10) %>%
  group_by(Genus) %>%
  summarise(mean_rel_abundance = mean(Abundance), .groups = "drop")

# Merge DA results with abundances
res_plot <- res_genus %>%
  rename(Genus = taxon) %>%
  left_join(genus_abund, by = "Genus")

# Label + color for significant taxa
res_plot <- res_plot %>%
  mutate(
    sig = case_when(
      q_astroPOS < 0.05 & lfc_astroPOS > 0 ~ "up",
      q_astroPOS < 0.05 & lfc_astroPOS < 0 ~ "down",
      TRUE ~ "ns"
    ),
    colour = case_when(
      sig == "up" ~ "#D14538",
      sig == "down" ~ "#3177B5",
      TRUE ~ "grey70"
    ),
    label = ifelse(sig != "ns", Genus, "")
  )


library(ggrepel)

p_da <- ggplot(res_plot,
               aes(x = mean_rel_abundance,
                   y = lfc_astroPOS,
                   label = label,
                   fill = colour)) +
  geom_point(size = 3, alpha = 0.8, colour="black", shape=21) +
  geom_text_repel(size = 3,
                  max.overlaps = Inf,
                  box.padding = 0.4,
                  point.padding = 0.2,
                  min.segment.length = 0) +
  scale_x_log10() +
  scale_fill_identity() +
  labs(
    x = "Mean Relative Abundance (log10)",
    y = "Log Fold-Change (in AstV+)"
  ) +
  theme_classic(base_size = 14) +
  theme(
    axis.title = element_text(face = "bold", size = 14),
    plot.title = element_text(face = "bold")
  )

p_da








