setwd("Z:/bilder_smagnagerfotobokser/Data/formatted_data/V_rodents_cameratraps_automatic_image_classification_lemming_blocks/automatic_image_classification")

dir()

ko16 <- read.table("V_rodents_cameratraps_automatic_image_classification_lemming_blocks_komagdalen_2016.txt", header=T, sep=";")
vj20 <- read.table("V_rodents_cameratraps_automatic_image_classification_lemming_blocks_vestre_jakobselv_2020.txt", header=T, sep=";")
