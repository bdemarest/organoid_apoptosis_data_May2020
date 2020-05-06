# Apoptosis data from timelapse green fluorescent imaging of brain organoids,
# ~19 hours (5 minute time steps) after transfer to matrigel.

# Experiment design:
# 3 genotypes - WT (WTC11), KMT2D null (KMT2D5.6), KMT2D null (KMT2D10.1)
# null lines were separately derived, but using the same CRISPR guide rna.


library(data.table)
library(here)
library(ggplot2)

# Added to github using these commands.
# usethis::use_git()
# usethis::use_github()


tab = fread("Mean_Fluorescence_all_samples_20April2020_utf8.txt")

tab[, replicate_id:=paste("rep", replicate, sep="")]
tab[, id2:=paste(genotype, replicate_id, sep="_")]
set(tab, j="replicate", value=NULL)

# Add mean-centered data column.
tab[, mf_centered:=mean_fluorescence - mean(mean_fluorescence), by=id2]

# Add 'standardized' data column.
tab[, mf_standardized:=(mean_fluorescence - 
                        mean(mean_fluorescence)) / 
                        sd(mean_fluorescence),
    by=id2]





by_id2 = tab[, list(sample_mean=mean(mean_fluorescence), 
                    sample_stdev=sd(mean_fluorescence),
                    sample_var=var(mean_fluorescence)), 
             by=list(id2, genotype)]


# Add difference column.
setorder(tab, sample_id)

tab[, diff_log10_mean_fluor:=c(NA_real_, diff(log10(mean_fluorescence))),
    by=id2]


dtab = tab[, list(mean_diff=mean(diff_log10_mean_fluor, na.rm=TRUE)), 
           by=list(id2, genotype)]


#---------------------------
# Mean vs. variance scatterplot.

s1 = ggplot(by_id2, aes(x=log10(sample_mean), 
                        y=log10(sample_var), 
                        fill=genotype)) +
     theme_bw() +
     geom_point(size=4, colour="grey20", shape=21) +
     labs(title="Mean vs variance relationship for n = 9 brain organoids.") +
     labs(subtitle="Each point is a summary of 229 values per sample over a 19 hour time series.")

ggsave("mean_v_variance_n9_organoids_n3genotypes_20200505.pdf", plot=s1,
       width=6, height=5)

#---------------------------
# Dotplots of average growth per 5 minute time step.

d1 = ggplot(dtab, aes(x=genotype, y=mean_diff, fill=genotype)) +
     theme_bw() +
     geom_point(size=4, colour="grey20", shape=21)
  

#---------------------------
# Histograms.

h1 = ggplot(tab, aes(x=log10(mean_fluorescence), fill=genotype)) +
     theme_bw() +
     geom_histogram(size=0.3, colour="grey20", bins=100) +
     facet_grid(id2 ~ .) +
     xlab("Log10(Mean Fluorescence)") +
     labs(title="Distribution of mean fluorescence values. n = 9 brain organoids.") +
     labs(subtitle="n = 229 mean fluorescence value per organoid (1 for each time point).")

ggsave("mean_fluorescence_histograms_n9_organoids_20200505.pdf", plot=h1,
       width=8, height=11)



h2 = ggplot(tab, aes(x=mf_centered, fill=genotype)) +
  theme_bw() +
  geom_histogram(size=0.3, colour="grey20", bins=50) +
  facet_wrap(~ id2) +
  xlab("Mean Fluorescence, Centered") +
  labs(title="Distribution of centered mean fluorescence values. n = 9 brain organoids.") +
  labs(subtitle="n = 229 mean fluorescence value per organoid (1 for each time point).")

  

ggsave("mean_fluorescence_histograms_centered_n9_organoids_20200505.pdf", 
       plot=h2, width=10, height=8)


h3 = ggplot(tab, aes(x=mf_standardized, fill=genotype)) +
  theme_bw() +
  geom_histogram(size=0.3, colour="grey20", bins=50) +
  facet_wrap(~ id2) +
  xlab("Mean Fluorescence, Standardized") +
  labs(title="Distribution of standardized mean fluorescence values. n = 9 brain organoids.") +
  labs(subtitle="n = 229 mean fluorescence value per organoid (1 for each time point).")

ggsave("mean_fluorescence_histograms_standardized_n9_organoids_20200505.pdf", 
       plot=h3, width=10, height=8)

#---------------------------

h4 = ggplot(tab, aes(x=diff_log10_mean_fluor, fill=genotype)) +
  theme_bw() +
  geom_histogram(size=0.3, colour="grey20", bins=100) +
  facet_grid(id2 ~ .) +
  xlab("Difference Log10(Mean Fluorescence)") +
  labs(title="Distribution of difference in mean fluorescence values. n = 9 brain organoids.") 

ggsave("diff_log10_mean_fluorescence_histograms_n9_organoids_20200505.pdf", plot=h4,
       width=8, height=11)





#---------------------------
# Time series plots.

p1 = ggplot(tab, aes(x=time, 
                     y=log10(mean_fluorescence), 
                     colour=genotype, 
                     group=id2)) +
     theme_bw() +
     geom_line(size=0.2) +
     geom_point(size=0.4)

ggsave("mean_fluorescence_log10_v_time_n9_organoids_20200505.pdf", 
       plot=p1, width=10, height=8)

p12 = ggplot(tab, aes(x=time, 
                     y=mf_standardized, 
                     colour=genotype, 
                     group=id2)) +
      theme_bw() +
      geom_line(size=0.2) +
      geom_point(size=0.4)

ggsave("mean_fluorescence_standardized_v_time_n9_organoids_20200505.pdf", 
       plot=p12, width=10, height=8)


#---------------------------

p2 = ggplot(tab, aes(x=time, 
                     y=diff_log10_mean_fluor, 
                     colour=genotype, 
                     group=id2)) +
     theme_bw() +
     geom_line(size=0.2) +
     geom_point(size=0.4)

ggsave("diff_mean_fluorescence_log10_v_time_n9_organoids_20200505.pdf", 
       plot=p2, width=10, height=8)

p3 = ggplot(tab, aes(x=time, 
                     y=diff_log10_mean_fluor, 
                     colour=genotype, 
                     group=id2)) +
     theme_bw() +
     geom_line(size=0.2) +
     geom_point(size=0.4) +
     geom_smooth(size=0.8) +
     facet_grid(id2 ~ .)

ggsave("diff_mean_fluorescence_log10_v_time_by_sampleid_n9_organoids_20200505.pdf", 
       plot=p3, width=10, height=8)

