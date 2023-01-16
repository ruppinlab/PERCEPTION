rm(list = ls())

# --- Functions ---
readSig <- function(x) {
  sig <- unique(scan(x, what = character(), skip = 2))
  return(sig)
}

# --- Code ---
# Apoptosis
apoptosis_up <- unique(c(readSig("data-raw/MSigDB/HERNANDEZ_MITOTIC_ARREST_BY_DOCETAXEL_1_UP.txt"),
                         readSig("data-raw/MSigDB/HERNANDEZ_MITOTIC_ARREST_BY_DOCETAXEL_2_UP.txt")))
apoptosis_down <- unique(c(readSig("data-raw/MSigDB/HERNANDEZ_MITOTIC_ARREST_BY_DOCETAXEL_1_DN.txt"),
                           readSig("data-raw/MSigDB/HERNANDEZ_MITOTIC_ARREST_BY_DOCETAXEL_2_DN.txt")))

sig_apoptosis_HERNANDEZ <- list(up = apoptosis_up, down = apoptosis_down)

# Senescence
senescence <- "data-raw/MSigDB/FRIDMAN_SENESCENCE_"
sig_senescence_FRIDMAN <- list(up = readSig(paste0(senescence, "UP.txt")),
                               down = readSig(paste0(senescence, "DN.txt")))

# Cell cycle
ccycle <- c(proliferation = "data-raw/MSigDB/WHITFIELD_CELL_CYCLE_LITERATURE.txt",
            G1S = "data-raw/MSigDB/WHITFIELD_CELL_CYCLE_G1_S.txt",
            G2M = "data-raw/MSigDB/WHITFIELD_CELL_CYCLE_G2_M.txt",
            G2 = "data-raw/MSigDB/WHITFIELD_CELL_CYCLE_G2.txt",
            MG1 = "data-raw/MSigDB/WHITFIELD_CELL_CYCLE_M_G1.txt",
            S = "data-raw/MSigDB/WHITFIELD_CELL_CYCLE_S.txt")

sig_proliferation_WHITFIELD <- list(up = readSig(ccycle["proliferation"]))
sig_G1S_WHITFIELD <- list(up = readSig(ccycle["G1S"]))
sig_G2M_WHITFIELD <- list(up = readSig(ccycle["G2M"]))
sig_G2_WHITFIELD <- list(up = readSig(ccycle["G2"]))
sig_MG1_WHITFIELD <- list(up = readSig(ccycle["MG1"]))
sig_S_WHITFIELD <- list(up = readSig(ccycle["S"]))

# EMT
EMT <- c(Alonso = "data-raw/MSigDB/ALONSO_METASTASIS_EMT_",
         Hallmarks = "data-raw/MSigDB/HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION.txt",
         GO = "data-raw/MSigDB/GOBP_EPITHELIAL_TO_MESENCHYMAL_TRANSITION.txt",
         negative_GO = "data-raw/MSigDB/GOBP_NEGATIVE_REGULATION_OF_EPITHELIAL_TO_MESENCHYMAL_TRANSITION.txt",
         positive_GO = "data-raw/MSigDB/GOBP_POSITIVE_REGULATION_OF_EPITHELIAL_TO_MESENCHYMAL_TRANSITION.txt",
         regulation_GO = "data-raw/MSigDB/GOBP_REGULATION_OF_EPITHELIAL_TO_MESENCHYMAL_TRANSITION.txt",
         Groger = "data-raw/Groger_2012/GROGER_")

sig_EMT_ALONSO <- list(up = readSig(paste0(EMT["Alonso"], "UP.txt")),
                       down = readSig(paste0(EMT["Alonso"], "DN.txt")))
sig_EMT_HALLMARKS <- list(up = readSig(EMT["Hallmarks"]))
sig_EMT_GO <- list(up = readSig(EMT["GO"]))
sig_negative_regulation_EMT_GO <- list(up = readSig(EMT["negative_GO"]))
sig_positive_regulation_EMT_GO <- list(up = readSig(EMT["positive_GO"]))
sig_regulation_EMT <- list(up = readSig(EMT["regulation_GO"]))

Groger_up <- unique(scan(paste0(EMT["Groger"], "UP.txt"), what = character()))
Groger_down <- unique(scan(paste0(EMT["Groger"], "DN.txt"), what = character()))
sig_EMT_GROGER_2012 <- list(up = Groger_up, down = Groger_down)

# All signatures
allsigs <- ls(pattern = "^sig_*")
pathways <- lapply(allsigs, get)
names(pathways) <- allsigs

# Save
usethis::use_data(pathways, internal = TRUE, overwrite = TRUE)
