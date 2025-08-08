library(hexSticker)
library(showtext)

font_add_google("Merriweather", "merri")
showtext_auto()

sticker("BAM-Logo.png", package = "BAMexploreR",
        p_color = "#1D5258",
        p_size = 16,               # package name font size
        s_x = 1, s_y = 0.75,      # sticker image x/y position
        s_width = 0.6,            # sticker image size
        h_fill = "#F0EDE5",       # hex background color
        h_color = "#1D5258",      # hex border color
        p_family = "merri",
         filename = "bam_hex.png")

