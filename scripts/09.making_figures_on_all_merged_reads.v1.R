#` # Making Figures
#' Author: Laura Los
#' Date: 2025-08-07

#' # set file path
path.7.figures <- file.path(here("data", "7.figures_all_merged_reads"))
if(!dir.exists(path.7.figures)) dir.create(path.7.figures)

#' # Import data
load(here("data", "6.phyloseq_on_all_merged_reads", "species.abundance.RData"))

# create figure of genera counts ####
taxa.long <- species.abundance %>% 
  pivot_longer(cols = -c(ASV, Kingdom, Phylum, Class, Order, Family, Genus, Species),
               names_to = "sample",
               values_to = "count") %>% 
  select(c(Genus, sample, count)) %>% 
  replace_na(replace = list(Genus =  "other"))

print(taxa.long, n = 50)

genus.counts <- taxa.long %>% 
  group_by(sample, Genus) %>% 
  summarise(across(where(is.numeric), ~ sum(.x, na.rm = TRUE))) %>% 
  ungroup()

print(genus.counts, n = 50)
str(genus.counts)
length(unique(genus.counts$sample))
length(unique(genus.counts$Genus))
length(unique(genus.counts$sample))*length(unique(genus.counts$Genus))

summarise

# ➋ Plot
## 1. Get the distinct taxon names (order doesn’t matter)
taxa_levels <- sort(unique(genus.counts$Genus))

## 2. Create a default palette for *all but* "other"
library(scales)        
n_coloured  <- length(taxa_levels) - 1          # how many need automatic hues?
auto_cols   <- hue_pal()(n_coloured)            # ggplot’s default discrete palette
## 3. Turn that into a named vector
names(auto_cols) <- taxa_levels[taxa_levels != "other"]

## 4. Add the grey override for "other"
my_palette <- c(auto_cols, other = "grey")

## 5. plot
genus.plot <- ggplot(genus.counts, aes(x = sample, y = count, fill = Genus)) +
  geom_bar(stat = "identity") +
  labs(x = "Sample", y = "Read count", fill = "Taxon") +
  scale_fill_manual(values = my_palette)+
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
genus.plot
filename = file.path(path.7.figures, "genus_counts.png")
ggsave(plot = genus.plot, filename, width = 8, height = 4)

# create figure of oomycete species counts ####
species.long <- species.abundance %>% 
  pivot_longer(cols = -c(ASV, Kingdom, Phylum, Class, Order, Family, Genus, Species),
               names_to = "sample",
               values_to = "count") %>% 
  select(c(Genus, Species, sample, count)) %>%   
  mutate(taxa = coalesce(Species, Genus)) %>% # if Species is NA, then it takes data from Genus column 
  replace_na(replace = list(taxa = "other"))

print(species.long, n = 100)

species.counts <- species.long %>% 
  group_by(sample, taxa) %>% 
  summarise(across(where(is.numeric), ~ sum(.x, na.rm = TRUE))) %>% 
  ungroup()

print(species.counts, n = 50)
str(species.counts)
length(unique(species.counts$sample))
length(unique(species.counts$taxa))
length(unique(species.counts$sample))*length(unique(species.counts$taxa))
str(species.counts)
# rows match the mutliplication

# ➋ Plot
## 1. Get the distinct taxon names (order doesn’t matter)
taxa_levels <- sort(unique(species.counts$taxa))

## 2. Create a default palette for *all but* "other"
library(scales)        
n_coloured  <- length(taxa_levels) - 1          # how many need automatic hues?
auto_cols   <- hue_pal()(n_coloured)            # ggplot’s default discrete palette
## 3. Turn that into a named vector
names(auto_cols) <- taxa_levels[taxa_levels != "other"]

## 4. Add the grey override for "other"
my_palette <- c(auto_cols, other = "grey")

## 5. plot
taxa.plot <- ggplot(species.counts, aes(x = sample, y = count, fill = taxa)) +
  geom_bar(stat = "identity") +
  labs(x = "Sample", y = "Read count", fill = "Taxon") +
  scale_fill_manual(values = my_palette)+
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
taxa.plot

## 6. save
filename = file.path(path.7.figures, "species_counts.png")
ggsave(plot = taxa.plot, filename, width = 8, height = 4)

# create figure of Kingdoms ####
kingdom.long <- species.abundance %>% 
  pivot_longer(cols = -c(ASV, Kingdom, Phylum, Class, Order, Family, Genus, Species),
               names_to = "sample",
               values_to = "count") %>% 
  select(c(Kingdom, sample, count)) %>%   
  mutate(taxa = Kingdom) %>% 
  replace_na(replace = list(taxa = "other"))

print(kingdom.long, n = 100)

kingdom.counts <- kingdom.long %>% 
  group_by(sample, taxa) %>% 
  summarise(across(where(is.numeric), ~ sum(.x, na.rm = TRUE))) %>% 
  ungroup()

print(kingdom.counts, n = 50)
str(kingdom.counts)
length(unique(kingdom.counts$sample))
length(unique(kingdom.counts$taxa))
length(unique(kingdom.counts$sample))*length(unique(kingdom.counts$taxa))
str(kingdom.counts)
# rows match the mutliplication

# ➋ Plot
## 1. Get the distinct taxon names (order doesn’t matter)
taxa_levels <- sort(unique(kingdom.counts$taxa))

## 2. Create a default palette for *all but* "other"
library(scales)        
n_coloured  <- length(taxa_levels) - 1          # how many need automatic hues?
auto_cols   <- hue_pal()(n_coloured)            # ggplot’s default discrete palette
## 3. Turn that into a named vector
names(auto_cols) <- taxa_levels[taxa_levels != "other"]

## 4. Add the grey override for "other"
my_palette <- c(auto_cols, other = "grey")

## 5. plot
taxa.plot <- ggplot(kingdom.counts, aes(x = sample, y = count, fill = taxa)) +
  geom_bar(stat = "identity") +
  labs(x = "Sample", y = "Read count", fill = "Taxon") +
  scale_fill_manual(values = my_palette)+
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
taxa.plot

## 6. save
filename = file.path(path.7.figures, "Kingdom_counts.png")
ggsave(plot = taxa.plot, filename, width = 8, height = 4)
