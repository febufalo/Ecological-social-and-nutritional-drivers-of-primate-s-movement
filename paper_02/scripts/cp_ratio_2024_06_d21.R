library(ggplot2)
df <- read.csv("cp_ratio_2024_06_d21.csv");head(df)
df$cp_ratio <- df$CPs/df$Scans;head(df)
# Remove rows where CPs equals 0
df <- subset(df, CPs != 0)

# Order the Behaviors factor levels
df$Behaviors <- factor(df$Behaviors, levels = c("Frugivory", "Locomotion", 
                                                "Resting", "Faunivory", 
                                                "Gummivory", "Vigilant", 
                                                "Other"))

# Reorder the levels of "Area" variable
df$Area <- factor(df$Area, levels = c("Continuous forest", "Medium fragment", 
                                      "Small fragment", "Riparian forest"))

# Create the bar plot with uniform bar widths and no space between x-axis and y-axis=0
#tiff("cp_ratio_2024_09_d10.tiff", width = 11000, height = 6500, res = 1000)
ggplot(df, aes(x = Behaviors, y = cp_ratio, fill = Area)) +
  scale_color_manual(values = c("Continuous forest" = "gray25", 
                                "Medium fragment" = "gray50", 
                                "Small fragment" = "gray80", 
                                "Riparian forest" = "gray100")) +
  scale_fill_manual(values = c("Continuous forest" = "gray25", 
                               "Medium fragment" = "gray50", 
                               "Small fragment" = "gray80", 
                               "Riparian forest" = "gray100")) +
  geom_col(position = position_dodge2(width = 0.9, preserve = "single"), 
           width = 0.7, color = "black", linewidth = 0.5) +
  geom_text(aes(label = CPs), 
            position = position_dodge2(width = 0.7, preserve = "single"), 
            vjust = -1.9, lineheight = 0.6,  # Adjust vjust for CPs
            size = 4.2, color = "black") +  # Add text labels for CPs on top of bars
  geom_text(aes(label = "\n___\n"), 
            position = position_dodge2(width = 0.7, preserve = "single"), 
            vjust = -0.3, lineheight = 0.6,  # Adjust vjust for "\n___\n"
            size = 4, color = "black") +  # Add text labels for "\n___\n" on top of bars
  geom_text(aes(label = Scans), 
            position = position_dodge2(width = 0.7, preserve = "single"), 
            vjust = -.2, lineheight = 0.6,  # Adjust vjust for Scans
            size = 4.2, color = "black") +  # Add text labels for Scans on top of bars
  labs(x = "",
       y = "CP Ratios") +
  scale_y_continuous(limits = c(0, 0.12), expand = c(0, 0)) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text.y = element_text(size = 18, color = "black"),
        axis.text.x = element_text(size = 18, color = "black"),
        axis.title.y = element_text(size = 20),
        legend.position = "bottom",
        legend.title = element_blank(),  # Remove the legend title
        legend.text = element_text(size = 20),
        legend.margin = margin(t = -10, r = 0, b = 0, l = 0))
#dev.off()
