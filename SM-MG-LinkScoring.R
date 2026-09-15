setwd("")

library(readxl)
library(readr)
library(dplyr)
library(stringr)
library(purrr)
library(tidyr)
library(circlize)

#-----------------------------
# 1. Read the xlsx file
#-----------------------------
features <- read_xlsx(
  "HighlyCorrelatedFeatures_3M_Sirius_formula_edit.xlsx"
) %>%
  mutate(across(starts_with("ClassyFire"), as.character))

classyfire_levels <- c(
  "ClassyFire_most_specific_class",
  "ClassyFire_level-5",
  "ClassyFire_subclass",
  "ClassyFire_class",
  "ClassyFire_superclass"
)

#-----------------------------
# 2. Read and parse the TXT file
#-----------------------------
txt_file <- "chamois_bgc_annotations.txt"
txt_lines <- read_lines(txt_file)

# Identify BGC headers (C#R#)
bgc_idx <- which(str_detect(txt_lines, "^C\\d+R\\d+$"))
bgc_idx <- c(bgc_idx, length(txt_lines) + 1)

parse_block <- function(start, end, lines) {
  bgc_code <- lines[start]
  
  # Collapse the whole block into one string
  block_text <- paste(lines[(start + 1):(end - 1)], collapse = " ")
  
  # Normalize whitespace and remove box-drawing clutter
  block_text <- block_text %>%
    str_replace_all("[│├└─]", " ") %>%
    str_replace_all("\\s+", " ")
  
  # Split on CHEMONTID entries
  entries <- str_split(block_text, "CHEMONTID:", simplify = FALSE)[[1]]
  
  # First element is garbage before the first CHEMONTID
  entries <- entries[-1]
  
  map_dfr(entries, function(entry) {
    
    # Extract class name (everything inside the first parentheses)
    class_match <- str_match(entry, "\\(([^)]+)\\)")
    if (is.na(class_match[1])) return(NULL)
    
    class_name <- class_match[2]
    
    # Extract the FIRST probability after the class name
    # (this handles probabilities on the same or later "line")
    prob_match <- str_match(entry, "\\):\\s*([0-9.]+)|\\s([0-9]\\.[0-9]+)")
    if (all(is.na(prob_match))) return(NULL)
    
    prob <- as.numeric(na.omit(prob_match[2:3])[1])
    
    tibble(
      BGC = bgc_code,
      classyfire_term = class_name,
      CHAMOIS_prob = prob
    )
  })
}

classyfire_long <- map_dfr(
  seq_len(length(bgc_idx) - 1),
  ~ parse_block(bgc_idx[.x], bgc_idx[.x + 1], txt_lines)
) %>%
  distinct()

#-----------------------------
# 3. Match xlsx vs TXT content
#-----------------------------
match_feature_to_bgc <- function(feature_row, classyfire_db, levels) {
  
  for (lvl in levels) {
    
    class_value <- feature_row[[lvl]]
    
    if (is.na(class_value) || class_value == "") next
    
    hits <- classyfire_db %>%
      filter(classyfire_term == class_value)
    
    if (nrow(hits) > 0) {
      return(
        hits %>%
          mutate(
            FeatureID = feature_row$FeatureID,
            matched_level = lvl,
            matched_class = class_value
          )
      )
    }
  }
  
  NULL
}

matches <- features %>%
  split(.$FeatureID) %>%
  map_dfr(
    ~ match_feature_to_bgc(
      feature_row = .x[1, ],
      classyfire_db = classyfire_long,
      levels = classyfire_levels
    )
  )

final_matches <- matches %>%
  left_join(
    features %>% select(FeatureID, Medium),
    by = "FeatureID"
  ) %>%
  select(
    FeatureID,
    Medium,
    matched_class,
    matched_level,
    BGC,
    CHAMOIS_prob
  ) %>%
  arrange(FeatureID, BGC)

final_matches %>%
  count(matched_level) %>%
  arrange(desc(n))

setdiff(features$FeatureID, final_matches$FeatureID)

# Write output TSV with all matches
write_tsv(final_matches, "FeatureID_BGC_all_class_iterative_matches.tsv")

#-----------------------------
# 4. Filter top matches CHAMOIS_prob > 0.7
#-----------------------------
#Filter matches by probability and level
links_filt <- final_matches %>%
  filter(
    CHAMOIS_prob > 0.7,
    matched_level != "ClassyFire_superclass"
  )

nrow(links_filt)

# Write output TSV with filtered matches
write_tsv(links_filt, "FeatureID_BGC_filtered_class_iterative_matches_.tsv")

#-----------------------------
# 5. Compute a support score
#-----------------------------
#Load supporting classes
support_raw <- read_xlsx(
  "HighlyCorrelatedFeatures_3M_Sirius_formula_allclassesprob.xlsx"
) %>%
  mutate(across(starts_with("ClassyFire"), as.character))

#Add an explicit row ID before parsing
support_raw <- support_raw %>%
  mutate(SupportRowID = row_number())

#Define hierarchy order
cf_levels <- tibble(
  level = c(
    "ClassyFire_level-6",
    "ClassyFire_level-5",
    "ClassyFire_subclass",
    "ClassyFire_class",
    "ClassyFire_superclass"
  ),
  level_rank = c(5, 4, 3, 2, 1),
  prob_col = paste0(level, "_probability")
)

#Pivot to long format
support_long <- map_dfr(seq_len(nrow(support_raw)), function(i) {
  
  row <- support_raw[i, ]
  
  map_dfr(seq_len(nrow(cf_levels)), function(j) {
    
    lvl  <- cf_levels$level[j]
    prob <- cf_levels$prob_col[j]
    
    class_val <- row[[lvl]]
    prob_val  <- row[[prob]]
    
    if (is.na(class_val) || class_val == "") return(NULL)
    
    tibble(
      FeatureID = row$FeatureID,
      Medium = row$Medium,
      SupportRowID = row$SupportRowID,
      classyfire_term = str_trim(class_val),
      P_feat = as.numeric(prob_val),
      level = lvl,
      level_rank = cf_levels$level_rank[j]
    )
  })
})

#Keep only the most specific class per row
support_most_specific <- support_long %>%
  group_by(FeatureID, Medium, SupportRowID) %>%
  slice_max(level_rank, n = 1, with_ties = FALSE) %>%
  ungroup()

###  Now to the proper computation, assuming we already have support_most_specific, classyfire_long and links_filt

###Fixing the ClassyFire_most_specific_class ambiguity
real_levels <- c(
  "ClassyFire_level-7",
  "ClassyFire_level-6",
  "ClassyFire_level-5",
  "ClassyFire_subclass",
  "ClassyFire_class",
  "ClassyFire_superclass"
)

feature_class_levels <- features %>%
  select(FeatureID, all_of(real_levels)) %>%
  pivot_longer(
    cols = all_of(real_levels),
    names_to = "true_level",
    values_to = "class_value"
  ) %>%
  filter(!is.na(class_value))

links_main <- links_filt %>%
  left_join(
    feature_class_levels,
    by = c(
      "FeatureID" = "FeatureID",
      "matched_class" = "class_value"
    )
  ) %>%
  mutate(
    matched_level_fixed = if_else(
      matched_level == "ClassyFire_most_specific_class",
      true_level,
      matched_level
    )
  )

#Assign numeric hierarchy levels consistently
hierarchy_map <- tibble(
  level = c(
    "ClassyFire_superclass",
    "ClassyFire_class",
    "ClassyFire_subclass",
    "ClassyFire_level-5",
    "ClassyFire_level-6",
    "ClassyFire_level-7"
  ),
  level_rank = c(1, 2, 3, 4, 5, 6)
)

#Attach this to the main match table
links_main <- links_main %>%
  left_join(
    hierarchy_map,
    by = c("matched_level_fixed" = "level")
  ) %>%
  rename(main_level_rank = level_rank) %>%
  rename(P_bgc_main = CHAMOIS_prob)

# Match supporting classes to BGC classes
support_hits <- support_most_specific %>%
  inner_join(
    classyfire_long,
    by = "classyfire_term"
  ) %>%
  rename(P_bgc = CHAMOIS_prob)

#Keep only support evidence for existing Feature–BGC links
support_hits <- support_hits %>%
  semi_join(
    links_main %>% select(FeatureID, BGC),
    by = c("FeatureID", "BGC")
  )

#Compute signed hierarchy distance (Δ)
support_hits <- support_hits %>%
  left_join(
    links_main %>% select(FeatureID, BGC, main_level_rank),
    by = c("FeatureID", "BGC")
  ) %>%
  mutate(
    delta_level = level_rank - main_level_rank
  )

#Direction-aware hierarchy weighting
support_hits <- support_hits %>%
  mutate(
    hierarchy_weight = case_when(
      delta_level <= -2 ~ 1.20,
      delta_level == -1 ~ 1.10,
      delta_level ==  0 ~ 1.00,
      delta_level ==  1 ~ 0.70,
      delta_level ==  2 ~ 0.45,
      delta_level >=  3 ~ 0.25
    )
  )

#Compute weighted support evidence per class
support_hits <- support_hits %>%
  mutate(
    evidence_raw = P_feat * P_bgc,
    evidence_weighted = pmin(evidence_raw * hierarchy_weight, 1)
  )

#Aggregate support per Feature–BGC (saturating model)
support_scores <- support_hits %>%
  group_by(FeatureID, BGC) %>%
  summarise(
    SupportScore = 1 - prod(1 - evidence_weighted),
    n_supporting_classes = n(),
    .groups = "drop"
  )

#-----------------------------
# 6. Compute a final composite score
#-----------------------------
###Add feature-side probability for the main class to the final score

#Extract feature-side probability for the matched class
feature_probs_long <- features %>%
  select(
    FeatureID,
    starts_with("ClassyFire_")
  ) %>%
  select(
    -ClassyFire_all_classifications
  ) %>%
  pivot_longer(
    cols = -FeatureID,
    names_to = "raw_name",
    values_to = "raw_value"
  ) %>%
  mutate(
    level = str_remove(raw_name, "_probability$"),
    type  = if_else(str_ends(raw_name, "_probability"),
                    "prob", "class")
  ) %>%
  select(-raw_name) %>%
  pivot_wider(
    names_from = type,
    values_from = raw_value
  ) %>%
  rename(
    classyfire_term = class,
    P_feat = prob
  ) %>%
  filter(!is.na(classyfire_term))


#Join feature-side probability to the links
links_main <- links_main %>%
  left_join(
    feature_probs_long,
    by = c(
      "FeatureID" = "FeatureID",
      "matched_class" = "classyfire_term",
      "matched_level_fixed" = "level"
    )
  )

#Incorporate feature probability into the main score
links_main <- links_main %>%
  mutate(
    P_bgc_main = as.numeric(P_bgc_main),
    P_feat     = as.numeric(P_feat),
    main_level_rank = as.numeric(main_level_rank)
  )

max_rank <- max(links_main$main_level_rank, na.rm = TRUE)

links_main <- links_main %>%
  mutate(
    hierarchy_weight = main_level_rank / max_rank,
    main_score = P_bgc_main * P_feat * hierarchy_weight
  )

#Final composite score
links_final <- links_main %>%
  left_join(support_scores, by = c("FeatureID", "BGC")) %>%
  mutate(
    SupportScore = replace_na(SupportScore, 0),
    final_score = main_score + 0.35 * SupportScore
  )

# Write output TSV with filtered matches
write_tsv(links_final, "Feature-BGC_links_scored_alpha035.tsv")

# Write output TSV with top scoring matches
links_final_0.8 <- links_final %>%
  filter(
    final_score > 0.8
  )

write_tsv(links_final_0.8, "Feature-BGC_TopScoringLinks.tsv")

#--------------------------------------------
# Stability-based α tuning
#--------------------------------------------

alpha_grid <- seq(0, 1, by = 0.05)

stability <- map_dfr(alpha_grid, function(alpha) {
  
  scores <- links_final %>%
    mutate(score = main_score + alpha * SupportScore)
  
  tibble(
    alpha = alpha,
    spearman = cor(
      scores$main_score,
      scores$score,
      method = "spearman",
      use = "complete.obs"
    ),
    n_above_1 = sum(scores$score > 1),
    score_sd = sd(scores$score, na.rm = TRUE)
  )
})

#How to pick α
#Choose:
# - the largest α before Spearman drops sharply
# - where n_above_1 starts increasing rapidly
# - where score variance inflates disproportionately
##### According to this, 0.35 is appropriate #####