library(targets)
library(cooccur)
library(reshape2)
library(ggplot2)
library(Hmsc)
library(vegan)

tar_load(comunidad_original_spcor)
tar_load(comunidad_muestreada_spcor_0cce3dc9)
tar_load(comunidad_muestreada_spcor_cf3c4d18)
tar_load(comunidad_muestreada_spcor_9e1f9bf7)
tar_load(comunidad_muestreada_spcor_affac5fc)
tar_load(comunidad_muestreada_spcor_728a5c6b)

sp_names <- paste0("sp", 1:10)

# coocurrence_pattern <- function(comunidad, sp_names) {
#   cooccur.original <- cooccur(mat = t(comunidad[,sp_names]), type = "spp_site", thresh = TRUE, spp_names = TRUE)
#   plot(cooccur.original)
# }
# 
# 
# 
# coocurrence_pattern(comunidad_original_spcor, sp_names)
# coocurrence_pattern(comunidad_muestreada_spcor_728a5c6b, sp_names)
# coocurrence_pattern(comunidad_muestreada_spcor_affac5fc, sp_names)
# coocurrence_pattern(comunidad_muestreada_spcor_9e1f9bf7, sp_names)
# coocurrence_pattern(comunidad_muestreada_spcor_cf3c4d18, sp_names)
# coocurrence_pattern(comunidad_muestreada_spcor_0cce3dc9, sp_names)


original_prob <- cooccur(mat = t(comunidad_original_spcor[,sp_names]), type = "spp_site", thresh = TRUE,
          spp_names = TRUE, only_effects = TRUE, eff_standard = TRUE, eff_matrix = TRUE)

prob_10 <- cooccur(mat = t(comunidad_muestreada_spcor_0cce3dc9[,sp_names]), type = "spp_site", thresh = TRUE,
        spp_names = TRUE, only_effects = TRUE, eff_standard = TRUE, eff_matrix = TRUE)

prob_25 <-  cooccur(mat = t(comunidad_muestreada_spcor_cf3c4d18[,sp_names]), type = "spp_site", thresh = TRUE,
                    spp_names = TRUE, only_effects = TRUE, eff_standard = TRUE, eff_matrix = TRUE)

prob_50 <- cooccur(mat = t(comunidad_muestreada_spcor_9e1f9bf7[,sp_names]), type = "spp_site", thresh = TRUE,
                   spp_names = TRUE, only_effects = TRUE, eff_standard = TRUE, eff_matrix = TRUE)

prob_75 <- cooccur(mat = t(comunidad_muestreada_spcor_affac5fc[,sp_names]), type = "spp_site", thresh = TRUE,
                   spp_names = TRUE, only_effects = TRUE, eff_standard = TRUE, eff_matrix = TRUE)

prob_90 <- cooccur(mat = t(comunidad_muestreada_spcor_728a5c6b[,sp_names]), type = "spp_site", thresh = TRUE,
                   spp_names = TRUE, only_effects = TRUE, eff_standard = TRUE, eff_matrix = TRUE)



plot_format <- function(prob, sample) {
  prob_matrix <- as.matrix(prob)
  prob_matrix[upper.tri(prob_matrix, diag = TRUE)] <- NA
  prob_matrix <- na.omit(melt(prob_matrix))
  prob_matrix <- cbind(prob_matrix, paste0(prob_matrix[,1], "-", prob_matrix[,2]))
  prob_matrix[,5]<-  sample
  colnames(prob_matrix) <- c("var1", "var2", "prob", "sp-sp", "sample")
  prob_matrix
}


# original_ggplot <- plot_format(original_prob, "original")
original_ggplot <- plot_format(original_prob, 1)
ggplot_10 <- plot_format(prob_10, 0.10)
ggplot_25 <- plot_format(prob_25, 0.25)
ggplot_50 <- plot_format(prob_50, 0.50)
ggplot_75 <- plot_format(prob_75, 0.75)
ggplot_90 <- plot_format(prob_90, 0.90)

coo_plot <- rbind(original_ggplot, ggplot_10, ggplot_25, ggplot_50, ggplot_75, ggplot_90)

ggplot(coo_plot) + 
  geom_line(aes(x= sample, y=prob, group = `sp-sp`, colour = prob), size = 1) +
  scale_colour_distiller(palette = "RdYlGn", values = c(0, 0.17, 1)) +
  theme_bw() + scale_x_reverse()

coo_plot$sample <- factor(coo_plot$sample)
coo_plot$sample <-factor(coo_plot$sample, levels = sort(levels(coo_plot$sample), decreasing= TRUE))

coo_plot$`sp-sp` <- factor(coo_plot$`sp-sp`)
coo_plot$`sp-sp` <-factor(coo_plot$`sp-sp`, levels = original_ggplot$`sp-sp` [order(original_ggplot$prob)])


ggplot(coo_plot) + 
  geom_tile(aes(x= sample, y=`sp-sp`, fill = prob)) +
  scale_fill_distiller(palette = "RdYlGn", values = c(0, 0.17, 1)) +
  theme_bw()

###########################

#original
tar_load(datos_artificiales_spcor)
original_corr = datos_artificiales_spcor$param$spCor
plot_original <- plot_format(original_corr, "original")
  
#samples
tar_load(modelo_hmsc_spcor_8cd04534)
tar_load(modelo_hmsc_spcor_dad2cdf7)
tar_load(modelo_hmsc_spcor_f525c10a)
tar_load(modelo_hmsc_spcor_05aaf6b0)
tar_load(modelo_hmsc_spcor_8f69786a)
tar_load(full_modelo_hmsc_spcor)


OmegaCor_0.1 <- computeAssociations(modelo_hmsc_spcor_8cd04534)
OmegaCor_0.25 <- computeAssociations(modelo_hmsc_spcor_dad2cdf7)
OmegaCor_0.50 <- computeAssociations(modelo_hmsc_spcor_f525c10a)
OmegaCor_0.75 <- computeAssociations(modelo_hmsc_spcor_05aaf6b0)
OmegaCor_0.90 <- computeAssociations(modelo_hmsc_spcor_8f69786a)
OmegaCor_full <- computeAssociations(full_modelo_hmsc_spcor)

corr_0.1 <- OmegaCor_0.1[[1]]$mean
corr_0.25 <- OmegaCor_0.25[[1]]$mean
corr_0.50 <- OmegaCor_0.50[[1]]$mean
corr_0.75 <- OmegaCor_0.75[[1]]$mean
corr_0.90 <- OmegaCor_0.90[[1]]$mean
corr_full <- OmegaCor_full[[1]]$mean


plot_0.1 <- plot_format(corr_0.1, 0.1)
plot_0.25 <- plot_format(corr_0.25, 0.25)
plot_0.50 <- plot_format(corr_0.50, 0.50)
plot_0.75 <- plot_format(corr_0.75, 0.75)
plot_0.90 <- plot_format(corr_0.90, 0.90)
plot_full <- plot_format(corr_full, 1)

corr_plot <- rbind(plot_0.1, plot_0.25, plot_0.50, plot_0.75, plot_0.90, plot_full, plot_original)

corr_plot$sample <- factor(corr_plot$sample)
corr_plot$sample <-factor(corr_plot$sample, levels = sort(levels(corr_plot$sample), decreasing= TRUE))
levels(corr_plot$sample)   

corr_plot$`sp-sp` <- factor(corr_plot$`sp-sp`)
corr_plot$`sp-sp` <-factor(corr_plot$`sp-sp`, levels = plot_full$`sp-sp` [order(plot_full$prob)])

ggplot(corr_plot) + 
  geom_raster(aes(x= sample, y=`sp-sp`, fill = prob)) +
  scale_fill_distiller(palette = "RdYlGn") +
  theme_bw() 
  

ggplot(corr_plot) + 
  geom_raster(aes(x= sample, y=`sp-sp`, fill = sign(prob))) +
  scale_fill_distiller(palette = "RdYlGn") +
  theme_bw() 

ggplot(corr_plot) + 
  geom_raster(aes(x= sample, y=`sp-sp`, fill = sign(prob))) +
  scale_fill_distiller(palette = "RdYlGn") +
  geom_text( x= 2, y = "sp8-sp5", label = "1.918333") +
  geom_text( x= 3, y = "sp8-sp5", label = "1.918333") +
  geom_text( x= 4, y = "sp8-sp5", label = "1.918333") +
  geom_text( x= 5, y = "sp8-sp5", label = "2.284362") +
  geom_text( x= 6, y = "sp8-sp5", label = "2.720173") +
  geom_text( x= 7, y = "sp8-sp5", label = "2.330424") +
  theme_bw() 

######PROCRUSTES##########

procrus_1 <- procrustes(sign(original_corr), sign(corr_full))
procrus_0.90 <- procrustes(sign(original_corr), sign(corr_0.90))
procrus_0.75 <- procrustes(sign(original_corr), sign(corr_0.75))
procrus_0.50 <- procrustes(sign(original_corr), sign(corr_0.50))
procrus_0.25 <- procrustes(sign(original_corr), sign(corr_0.25))
procrus_0.10 <- procrustes(sign(original_corr), sign(corr_0.1))


muestreos <- c(0.1, 0.25, 0.50, 0.75, 0.90, 1)
rmse <- c(summary(procrus_0.10)$rmse, summary(procrus_0.25)$rmse, summary(procrus_0.50)$rmse,
                summary(procrus_0.75)$rmse, summary(procrus_0.90)$rmse, summary(procrus_1)$rmse)

procrustes <- data.frame(cbind(muestreos, rmse))

procrustes$muestreos <- factor(procrustes$muestreos)
procrustes$muestreos <-factor(procrustes$muestreos, levels = sort(levels(procrustes$muestreos), decreasing= TRUE))

ggplot(procrustes) +
  geom_point(aes(x= muestreos, y= rmse)) + 
  theme_bw()




