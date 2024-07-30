## July 4th 1b!!

res = July4$res
df.res = res %>% dplyr::select(taxon, contains("DiabetesType")) %>%
  dplyr::filter(diff_DiabetesTypelate == 1) 

res = July4$res
df.res = res %>% dplyr::select(taxon, contains("Diabetes")) %>%
  dplyr::filter(diff_Diabetesno == 1) 


## need to remove mouse inoculum!!

## JUNE 14TH PRIMARY ANALYSIS

# a -> Comparison 1a, otu level
a.res = a$res
df.a.res = a.res  %>% dplyr::select(taxon, contains("DiabetesType")) %>%
  dplyr::filter(diff_DiabetesTypelate == 1) 

# b -> Comparison 1b, otu level
b.res = b$res
df.b.res = b.res  %>% dplyr::select(taxon, contains("DiabetesType")) %>%
  dplyr::filter(diff_DiabetesTypelate == 1)

# c -> Comparison 2, otu level
c.res = c$res
df.c.res = c.res  %>% dplyr::select(taxon, contains("Diabetesno")) %>%
  
# d -> Comparison 1b, genus level
d.res = d$res
df.d.res = d.res %>% dplyr::select(taxon, contains("DiabetesType")) %>%
  dplyr::filter(diff_DiabetesTypelate == 1) 

# e -> Comparison 1b, family level
e.res = e$res
df.e.res = e.res %>% dplyr::select(taxon, contains("DiabetesType")) %>%
  dplyr::filter(diff_DiabetesTypelate == 1) 

# f -> Comparison 1b, order level
f.res = f$res
df.f.res = f.res %>% dplyr::select(taxon, contains("DiabetesType")) %>%
  dplyr::filter(diff_DiabetesTypelate == 1) 


# Comparison 2 genus level no sig

# Comparison 2, family level no sig

# Comparison 2, order level no sig



%>% filter(diff_DiabetesTypelate == 1)

res_prim = ancom_output$res # ancom_output with ps2.test.norm, distinct
res_prim_2 = ancom_output_2$res # ancom_output with ps2.norm, dinstinct
res_prim_all = ancom_output_all$res # ancom_output with ps2.norm, all not distinct
res_prim_all_family = ancom_output_all_family$res # tax level at family
res_prim_all_genus = ancom_output_all_genus$res # tax level at family
res_prim_all_order = ancom_output_all_order$res # tax level at family
res_prim_all_class = ancom_output_all_class$res # tax level at family
res_prim_all_phylum = ancom_output_all_phylum$res # tax level at family

res2_otu = ancom_output2$res
res2_phy = ancom_output2_phy$res
res2_gen = ancom_output2_gen$res
res2_ord = ancom_output2_ord$res
res2_fam = ancom_output2_fam$res

res_prim_J17 = ancom_output_J17
res2_J17 = ancom_output_J17$res
df_DT_J17 = res2_J17 %>% dplyr::select(taxon, contains("DiabetesType")) %>%
  dplyr::filter(diff_DiabetesTypelate == 1) 

res_prim_J18 = ancom_output_J18
res2_J18 = ancom_output_J18$res
df_DT_J18 = res2_J18 %>% dplyr::select(taxon, contains("DiabetesType")) %>%
  dplyr::filter(diff_DiabetesTypelate == 1) 

## Example from online tutorial

# results for diabetes type
df_DT2_phy = res2_phy %>%
  dplyr::select(taxon, contains("Diabetes")) %>%
  dplyr::filter(diff_Diabetesno == 1) 

df_fig_DT1 = df_DT_all %>%
  dplyr::filter(diff_DiabetesTypelate == 1) %>%
  dplyr::mutate(lfc1 = ifelse(diff_DiabetesTypelate == 1, 
                              round(lfc_DiabetesTypelate, 2), 0),
                lfc2 = ifelse(diff_bmilean == 1, 
                              round(lfc_bmilean, 2), 0)) %>%
  tidyr::pivot_longer(cols = lfc1:lfc2, 
                      names_to = "group", values_to = "value") %>%
  dplyr::arrange(taxon)

df_fig_bmi2 = df_bmi %>%
  dplyr::filter(diff_bmilean == 1 | 
                  diff_bmioverweight == 1) %>%
  dplyr::mutate(lfc1 = ifelse(passed_ss_bmioverweight == 1 & diff_bmioverweight == 1, 
                              "aquamarine3", "black"),
                lfc2 = ifelse(passed_ss_bmilean == 1 & diff_bmilean == 1, 
                              "aquamarine3", "black")) %>%
  tidyr::pivot_longer(cols = lfc1:lfc2, 
                      names_to = "group", values_to = "color") %>%
  dplyr::arrange(taxon)

df_fig_bmi = df_fig_bmi1 %>%
  dplyr::left_join(df_fig_bmi2, by = c("taxon", "group"))

df_fig_bmi$group = recode(df_fig_bmi$group, 
                          `lfc1` = "Overweight - Obese",
                          `lfc2` = "Lean - Obese")
df_fig_bmi$group = factor(df_fig_bmi$group, 
                          levels = c("Overweight - Obese",
                                     "Lean - Obese"))

lo = floor(min(df_fig_bmi$value))
up = ceiling(max(df_fig_bmi$value))
mid = (lo + up)/2
fig_bmi = df_fig_bmi %>%
  ggplot(aes(x = group, y = taxon, fill = value)) + 
  geom_tile(color = "black") +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       na.value = "white", midpoint = mid, limit = c(lo, up),
                       name = NULL) +
  geom_text(aes(group, taxon, label = value, color = color), size = 4) +
  scale_color_identity(guide = FALSE) +
  labs(x = NULL, y = NULL, title = "Log fold changes as compared to obese subjects") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))
fig_bmi

