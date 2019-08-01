filter_group <- function (df, group_number) {
    df %>%
    filter(!is.na(group_number))  %>%
    group_by(group_number) %>%
    summarize(n = n())
}
#Count number of unique species in each group
grp_01 <- spp %>%
  filter(!is.na(groups01))

nspp_grp <- nrow(grp_01)

grp_01 <-  grp_01 %>%
  group_by(groups01) %>%
  summarize(n = n())

ngrp_grp <-nrow(grp_01)
```
