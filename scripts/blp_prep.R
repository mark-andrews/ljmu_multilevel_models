library(osfr)
library(tidyverse)

# download BLP
# blp_osf_project <- 'b5sdk'
# blp_osf_node <- osf_retrieve_node(blp_osf_project)
# osf_ls_files(blp_osf_node) |> osf_download(path = 'data/blp')

# load trials
load('data/blp/blp-trials.Rdata')

blp_trials_df <- blp.trials %>% 
  select(subject = participant, lexicality, spelling, accuracy, rt = rt.raw) %>% 
  as_tibble() %>% 
  mutate(lex = case_when(
    lexicality == 'N' ~ 'nonword',
    lexicality == 'W' ~ 'word'),
    accuracy = accuracy == 1
  ) %>% select(subject, item = spelling, lex, accuracy, rt) %>% 
  filter(rt > 200) %>% 
  drop_na()


load('data/blp/blp-stimuli.Rdata')

blp_trials <- left_join(
  blp_trials_df, 
  select(blp.stimuli, item = spelling, freq = bnc.frequency.million, length = nletters),
  by = 'item')

blp_trials <- blp_trials %>% 
  filter(freq > 0) %>% 
  sample_frac(0.5) %>% 
  add_count(item) %>% 
  filter(n > 20)

save(file = 'data/blp_trials.Rda', blp_trials)
