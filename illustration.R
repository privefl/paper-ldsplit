corr <- readRDS(url("https://www.dropbox.com/s/65u96jf7y32j2mj/spMat.rds?raw=1"))
corrT <- as(corr, "dgTMatrix")
THR_R2 <- 0.05
ind <- 81:103
library(magrittr)
df2 <- corrT[ind, ind] %>%
  { dplyr::filter(data.frame(i = .@i + ind[1], j = .@j + ind[1], r2 = .@x^2),
                  r2 >= THR_R2) }


library(ggplot2)
ggplot(df2, aes(i, j)) +
  geom_raster(aes(fill = r2), alpha = 0.6) +
  scale_fill_viridis_b(direction = -1, breaks = c(0.02, 0.1, 0.3, 0.8)) +
  geom_text(aes(label = label),
            data = data.frame(i = 80, j = c(81, 91), label = c("i", "j"))) +
  geom_text(aes(label = label),
            data = data.frame(i = c(88, 99, 100), j = c(83, 83, 94),
                              label = c("C(i, k)", "E(i, j)", "C(j+1, k-1)"))) +
  coord_equal() +
  scale_x_continuous(breaks = 1:1000 - 0.5, minor_breaks = NULL) +
  scale_y_reverse(breaks = 1:1000 - 0.5, minor_breaks = NULL) +
  theme_minimal() +
  theme(axis.text.y = element_blank(), axis.text.x = element_blank()) +
  labs(x = NULL, y = NULL, fill = expression(r^2)) +
  geom_hline(yintercept = c(80.5, 91.5), linetype = 2) +
  geom_vline(xintercept = c(80.5, 91.5), linetype = 2)

