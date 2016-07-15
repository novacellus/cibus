# List of word whose definition in LS contains the word "food" etc.
food_LS <- gloss_nostop_wide[which(gloss_nostop_wide[,food_NN > 0])]$L1
feast_LS <- gloss_nostop_wide[which(gloss_nostop_wide[,feast_NN > 0])]$L1
