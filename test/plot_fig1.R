devtools::install_github("johannesbjork/LaCroixColoR")
library(LaCroixColoR)
lacroix_palette("PeachPear", type = "discrete")

progs <- c("Helmsman", "SomaticSignatures", "deconstructSigs", "signeR")

performance_df <- data.frame(
  program=rep(progs, times=2), 
  metric=rep(c("processing time (s)", "memory usage (MB)"), each=4), 
  vals=c(8,227,2376,1740,125,18000,7500,10200))

performance_df$program <- factor(progs, levels=progs)
# levels(performance_df$progs) <- progs

ggplot(performance_df, aes(x=program, y=vals, fill=program))+
  geom_col()+
  facet_wrap(~metric, scales="free", strip.position="left")+
  scale_y_continuous(expand=c(0,0))+
  # coord_flip()+
  # scale_y_log10(expand=c(0,0), breaks=10^seq(0,4))+
  # scale_fill_manual(values=lacroix_palettes$Pamplemousse[1:4])+
  scale_fill_brewer(palette="Spectral", direction=-1)+
  ylab(NULL)+
  theme_classic()+
  theme(legend.position="bottom",
        strip.background = element_blank(),
        strip.placement = "outside",
        axis.text.y=element_text(size=14),
        strip.text=element_text(size=16),
        legend.title = element_blank(),
        legend.text = element_text(size=14),
        axis.ticks.x=element_blank(),
        axis.text.x=element_blank(),
        axis.title.x=element_blank())
ggsave("helmsman/assets/fig1_v3.png", width=12, height=6)
