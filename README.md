# SM-MG-LinkScoring
Chemical Ontology based Scoring Framework for Specialized Metabolism Metabologenomic Association

### Objective

To quantify the strength of potential associations between **MS/MS metabolomic features** and **biosynthetic gene clusters (BGCs)** using **hierarchy-aware chemical class evidence**, integrating:

- main ClassyFire class matches
- supporting class matches across hierarchy levels
- classification probabilities from both the MS/MS feature side and the BGC side

### In a nutshell
The scoring framework integrates probabilistic chemical-class annotations from metabolomic features and BGCs to prioritize potential BGC–metabolite associations. It combines prediction confidence with the specificity of shared chemical classes, while incorporating additional annotation evidence to derive an overall association score.
