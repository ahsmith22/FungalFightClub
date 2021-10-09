final_data_ca = read_csv("~/FungalFightClub/finaldata1.csv")

# format and clean data
final_data_ca = final_data_ca[,c(1:4,7)]
final_data_ca$cont_ID = ifelse(substr(final_data_ca$Plate_ID,2,2) == "v", "v",
                                ifelse(substr(final_data_ca$Plate_ID,2,2) == "s", "s", 0)) # v for svs control, s for single control, 0 for pair
final_data_cap = final_data_ca[final_data_ca$cont_ID==0,]
final_data_cac = final_data_ca[final_data_ca$cont_ID=="v"|final_data_ca$cont_ID=="s",]

# format pairs
final_data_cap$Fun_ID = toupper(final_data_cap$Plug_ID) # Fungal ID (A, C, H, L, or P)
final_data_cap$Plug_ID = paste(final_data_cap$Plug_ID, final_data_cap$Plate_ID, sep = "") # For pair plates, ID of plug
final_data_cap$PH_ID = substr(final_data_cap$Plate_ID, 3, 3) # pH of medium
final_data_cap$opp_ID = ifelse(final_data_cap$Fun_ID == substr(final_data_cap$Plate_ID, 1, 1), substr(final_data_cap$Plate_ID, 2, 2), substr(final_data_cap$Plate_ID, 1, 1)) # For pair plates, ID of opponent plug

# format controls
final_data_cac$Fun_ID = substr(final_data_cac$Plate_ID, 1, 1)
final_data_cac$cont_ID = substr(final_data_cac$Plate_ID, 2, 2) # either s or v
final_data_cac$PH_ID = substr(final_data_cac$Plate_ID, 3, 3)
final_data_cac$Plug_ID = paste(substr(final_data_cac$Plate_ID, 1, 1), final_data_cac$Plate_ID, sep = "")
final_data_cac$opp_ID = ifelse(final_data_cac$cont_ID=="v", substr(final_data_cac$Plate_ID, 1, 1), "N") # v for svs and N for none

# combine
final_data_ca = rbind(final_data_cap, final_data_cac)
final_data_ca = final_data_ca[,c(1,7,6,8,2,9,3,5)]

for (u in 1:length(value)) { # iterate for r and K
  for (x in 1:length(fungi)) { # iterate for each fungi
    for (w in 1:length(pH)) { # iterate for each pH
      funT = final_data_ca[final_data_ca$Fun_ID==fungi[x]&final_data_ca$PH_ID==pH[w],] # this Fun and pH
      sinT = funT[funT$cont_ID=="s",] # all single control for this Fun and pH
      svsT = funT[funT$cont_ID=="v",] # all svs for this Fun and pH
      parT = funT[funT$cont_ID==0,] # all pairs for this Fun and pH
      denomT = as.data.frame(matrix(nrow=nrow(parT),ncol=nrow(svsT))) # matrix rows = svs, cols = 1 opp pair
      for (y in 1:nrow(svsT)) {
        for (z in 1:nrow(parT)) {
          ainter = parT[z,value[u]]
          aintra = svsT[y,value[u]]
          denomT[z,y] = sqrt(ainter*aintra)  # compute denominator
        }
      }
      denomT$mean = rowMeans(denomT, na.rm=TRUE) # calculate denominator mean for each pair
      parT$denom = denomT$mean # add means back to pair df
      compabilT = as.data.frame(matrix(nrow=nrow(parT),ncol=nrow(sinT)))
      for (a in 1:nrow(sinT)) {
        for (b in 1:nrow(parT)) {
          compabilT[b,a] = abs((parT$denom[b])-((parT$denom[b])/sinT[a,value[u]]))  #(((sinT[a,value[u]])^2)-sinT[a,value[u]])/(parT$denom[b])
        }
      }
      colnam = paste("mean",value[u],sep="_")
      compabilT[[colnam]] = rowMeans(compabilT, na.rm=TRUE)
      compabilT$Plug_ID = parT$Plug_ID
      for (c in 1:nrow(compabilT)) {
        final_data_ca[final_data_ca$Plug_ID==compabilT$Plug_ID[c],colnam] = compabilT[c,ncol(compabilT)-1] # store values in parent df
      }
    }
  }
}

# caplots = vector(mode="list", length=5) # vector for storing individual ca plots
# # create plots of ca by r
# for(i in 1:length(fungi)){
#   subdat <- final_data_ca[final_data_ca$Fun_ID==fungi[i],]
#   caplots[[i]] = ggplot(subdat, aes(opp_ID,mean_r)) +
#     theme_cowplot() +
#     geom_boxplot(outlier.alpha = 0, aes(fill=PH_ID), color="black") +
#     geom_point(position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.9), aes(opp_ID,mean_r,fill = PH_ID,alpha=0.7)) +
#     scale_y_continuous(breaks=c(-1,-0.75,-0.5,-0.25,0,0.25,0.5,0.75,1)) +
#     geom_hline(yintercept=0, size=0.35) +
#     geom_vline(xintercept = c(1.5,2.5,3.5,4.5), linetype=2, size=0.25) +
#     theme(legend.position = "none", axis.line.x = element_blank(), axis.ticks.x = element_blank(), plot.margin = unit(c(1,1,0.5,1), "cm"), axis.title.x = element_blank(), axis.text.x = element_blank()) +
#     ylab("ca (r)")
# }
# caplot_r=plot_grid(caplots[[1]],caplots[[2]],caplots[[3]],caplots[[4]],caplots[[5]],nrow=5,ncol=1,labels=labs$main,label_size=14,label_x=labs$x,label_y=labs$y,rel_heights= c(1,1,1,1,1))
# 
# caplot_r

caplot_r = ggplot(final_data_ca[final_data_ca$cont_ID==0,], aes(Fun_ID, mean_r)) + 
  theme_cowplot() +
  geom_boxplot(outlier.alpha = 0, aes(fill=PH_ID), color="black", width=0.75) +
  geom_point(position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.9), aes(Fun_ID,mean_r,fill = PH_ID,alpha=0.7)) +
  ylab(bquote('growth rate ('*cm^2*'/day)')) +
  scale_x_discrete(labels = labs$nams) +
  geom_vline(xintercept = c(1.5,2.5,3.5,4.5), linetype=2, size=0.25) +
  theme(legend.position = "none", axis.line.x = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_blank(), plot.margin = unit(c(1,0,1,1), "cm"), axis.text.x = element_text(face="italic", angle=330, size=10), axis.title.y = element_text(size = 11))
caplot_r
# ggsave("caplot_r.png")

TukeyHSD(aov(lm(mean_r~Fun_ID*PH_ID, data=final_data_ca)))

caplot_K = ggplot(final_data_ca[final_data_ca$cont_ID==0,], aes(Fun_ID, mean_K)) + 
  theme_cowplot() +
  geom_boxplot(outlier.alpha = 0, aes(fill=PH_ID), color="black", width=0.75) +
  geom_point(position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.9), aes(Fun_ID,mean_K,fill = PH_ID,alpha=0.7)) +
  ylab(bquote('colony size ('*cm^2*')')) +
  scale_x_discrete(labels = labs$nams) +
  geom_vline(xintercept = c(1.5,2.5,3.5,4.5), linetype=2, size=0.25) +
  theme(legend.position = "none", axis.line.x = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_blank(), plot.margin = unit(c(1,0,1,1), "cm"), axis.text.x = element_text(face="italic", angle=330, size=10), axis.title.y = element_text(size = 11))
caplot_K
# ggsave("caplot_K.png")

TukeyHSD(aov(lm(mean_K~Fun_ID*PH_ID, data=final_data_ca)))


sinplot_r = ggplot(final_data_ca[final_data_ca$cont_ID=='s',], aes(Fun_ID, r)) +
  theme_cowplot() +
  geom_boxplot(outlier.alpha = 0, aes(fill=PH_ID), color="black", width=0.75) +
  geom_point(position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.9), aes(Fun_ID,r,fill = PH_ID,alpha=0.7)) +
  ylab(bquote('growth rate ('*cm^2*'/day)')) +
  scale_x_discrete(labels = labs$nams) +
  geom_vline(xintercept = c(1.5,2.5,3.5,4.5), linetype=2, size=0.25) +
  theme(plot.background = element_rect(fill = "white"), legend.position = "none", axis.line.x = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_blank(), plot.margin = unit(c(1,0,1,1), "cm"), axis.text.x = element_text(face="italic", angle=330, size=10), axis.title.y = element_text(size = 11))
sinplot_r
# ggsave("sinplot_r.png")

TukeyHSD(aov(lm(r~Fun_ID*PH_ID, data=final_data_ca[final_data_ca$cont_ID=="s",])))

sinplot_K = ggplot(final_data_ca[final_data_ca$cont_ID=='s',], aes(Fun_ID, K)) +
  theme_cowplot() +
  geom_boxplot(outlier.alpha = 0, aes(fill=PH_ID), color="black", width=0.75) +
  geom_point(position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.9), aes(Fun_ID,K,fill = PH_ID,alpha=0.7)) +
  ylab(bquote('colony size ('*cm^2*')')) +
  scale_x_discrete(labels = labs$nams) +
  geom_vline(xintercept = c(1.5,2.5,3.5,4.5), linetype=2, size=0.25) +
  theme(plot.background = element_rect(fill = "white"), legend.position = "none", axis.line.x = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_blank(), plot.margin = unit(c(1,0,1,1), "cm"), axis.text.x = element_text(face="italic", angle=330, size=10), axis.title.y = element_text(size = 11))
sinplot_K
# ggsave("sinplot_K.png")

TukeyHSD(aov(lm(K~Fun_ID*PH_ID, data=final_data_ca[final_data_ca$cont_ID=="s",])))


sin_r = summarySE(final_data_ca[final_data_ca$cont_ID=="s",], "r", c("Fun_ID", "PH_ID"), na.rm = TRUE)
sin_K = summarySE(final_data_ca[final_data_ca$cont_ID=="s",], "K", c("Fun_ID", "PH_ID"), na.rm = TRUE)
final_data_cap = final_data_ca[final_data_ca$cont_ID==0,]
final_data_cap = left_join(final_data_cap, sin_r[,c(1,2,4)], by = c("Fun_ID", "PH_ID"), suffix = c("","sin"))
final_data_cap = left_join(final_data_cap, sin_K[,c(1,2,4)], by = c("Fun_ID", "PH_ID"), suffix = c("","sin"))
final_data_cap$pair_ID = substr(final_data_cap$Plate_ID,1,2)
final_data_cap = left_join(final_data_cap, dist, by = "pair_ID")


mod_CAr_r = lmer(mean_r ~ rsin + PH_ID + (1|Fun_ID) + (1|Plate_ID), data = final_data_cap)
summary(mod_CAr_r)
class(mod_CAr_r) = "lmerMod"
r.squaredGLMM(mod_CAr_r)
CAr_r = as.data.frame(r.squaredGLMM(mod_CAr_r))
CAr_r_R2 = round(CAr_r$R2c,4)

mod_CAr_dist = lmer(mean_r ~ dist + PH_ID + (1|Fun_ID) + (1|Plate_ID), data = final_data_cap)
summary(mod_CAr_dist)
class(mod_CAr_dist) = "lmerMod"
r.squaredGLMM(mod_CAr_dist)
CAr_dist = as.data.frame(r.squaredGLMM(mod_CAr_dist))
CAr_dist_R2 = round(CAr_dist$R2c,4)

stargazer(mod_CAr_r, mod_CAr_dist,
          type = "text",
          #out = "mod_CAr.htm",
          dep.var.labels = c("CAr"),
          covariate.labels = c("Growth Rate", "Phylogenetic Distance", "pH 7", "constant"),
          column.labels = c("Trait-based", "Phylogeny-based"),
          # ci = TRUE,
          add.lines = list(c("Conditional pseudo-R2", CAr_r_R2, CAr_dist_R2)),
          summary = TRUE,
          digits = 3
)

mod_CAK_K = lmer(mean_K ~ Ksin + PH_ID + (1|Fun_ID) + (1|Plate_ID), data = final_data_cap)
summary(mod_CAK_K)
class(mod_CAK_K) = "lmerMod"
r.squaredGLMM(mod_CAK_K)
CAK_K = as.data.frame(r.squaredGLMM(mod_CAK_K))
CAK_K_R2 = round(CAK_K$R2c,4)

mod_CAK_dist = lmer(mean_K ~ dist + PH_ID + (1|Fun_ID) + (1|Plate_ID), data = final_data_cap)
summary(mod_CAK_dist)
class(mod_CAK_dist) = "lmerMod"
r.squaredGLMM(mod_CAK_dist)
CAK_dist = as.data.frame(r.squaredGLMM(mod_CAK_dist))
CAK_dist_R2 = round(CAK_dist$R2c,4)

stargazer(mod_CAK_K, mod_CAK_dist,
          type = "text",
          #out = "mod_CAr.htm",
          dep.var.labels = c("CAK"),
          covariate.labels = c("Colony Size", "Phylogenetic Distance", "pH 7", "constant"),
          column.labels = c("Trait-based", "Phylogeny-based"),
          # ci = TRUE,
          add.lines = list(c("Conditional pseudo-R2", CAK_K_R2, CAK_dist_R2)),
          summary = TRUE,
          digits = 3
)


pred_CAr_r = ggpredict(mod_CAr_r, terms = c("rsin", "PH_ID"))
plot_CAr_r = ggplot(pred_CAr_r) + 
  geom_line(aes(x = x, y = predicted, color = group)) +          # slope
  geom_ribbon(aes(x = x, ymin = predicted - std.error, ymax = predicted + std.error, fill = group), 
              alpha = 0.2) +  # error band
  geom_point(data = final_data_cap,                      # adding the raw data 
             aes(x = rsin, y = mean_r, color = PH_ID)) + 
  labs(x = "control growth rate", y = "CA_r") + 
  theme_cowplot() +
  theme(plot.background = element_rect(fill = "white"))
plot_CAr_r
# ggsave("plot_CAr_r.png")

pred_CAr_dist = ggpredict(mod_CAr_dist, terms = c("dist", "PH_ID"))
plot_CAr_dist = ggplot(pred_CAr_dist) + 
  geom_line(aes(x = x, y = predicted, color = group)) +          # slope
  geom_ribbon(aes(x = x, ymin = predicted - std.error, ymax = predicted + std.error, fill = group), 
              alpha = 0.2) +  # error band
  geom_point(data = final_data_cap,                      # adding the raw data 
             aes(x = dist, y = mean_r, color = PH_ID)) + 
  labs(x = "phylo dist", y = "CA_r") + 
  theme_cowplot() +
  theme(plot.background = element_rect(fill = "white"))
plot_CAr_dist
# ggsave("plot_CAr_dist.png")

pred_CAK_K = ggpredict(mod_CAK_K, terms = c("Ksin", "PH_ID"))
plot_CAK_K = ggplot(pred_CAK_K) + 
  geom_line(aes(x = x, y = predicted, color = group)) +          # slope
  geom_ribbon(aes(x = x, ymin = predicted - std.error, ymax = predicted + std.error, fill = group), 
              alpha = 0.2) +  # error band
  geom_point(data = final_data_cap,                      # adding the raw data 
             aes(x = Ksin, y = mean_K, color = PH_ID)) + 
  labs(x = "control colony size", y = "CA_K") + 
  theme_cowplot() +
  theme(plot.background = element_rect(fill = "white"))
plot_CAK_K
# ggsave("plot_CAK_K.png")

pred_CAK_dist = ggpredict(mod_CAK_dist, terms = c("dist", "PH_ID"))
plot_CAK_dist = ggplot(pred_CAK_dist) + 
  geom_line(aes(x = x, y = predicted, color = group)) +          # slope
  geom_ribbon(aes(x = x, ymin = predicted - std.error, ymax = predicted + std.error, fill = group), 
              alpha = 0.2) +  # error band
  geom_point(data = final_data_cap,                      # adding the raw data 
             aes(x = dist, y = mean_K, color = PH_ID)) + 
  labs(x = "phylo dist", y = "CA_K") + 
  theme_cowplot() +
  theme(plot.background = element_rect(fill = "white"))
plot_CAK_dist
# ggsave("plot_CAK_dist.png")

