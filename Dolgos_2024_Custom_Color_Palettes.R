library(paletteer)

# Sample_Description
pal_Sample_Name <- c(
  "P22-36_Org_Mat" = "#7EC3E5FF",    
  "P23-36_Org_ECMf" = "#E89242FF",  
  "P23-36_Org_Mat" = "#499A55",   
  "P23-51_Org_ECMf" = "#990F0FFF",  
  "P23-53_Org_ECMf" = "#9370DBFF",
  "P24-06_PAR_Tum" = "#3F0036FF", 
  "P24-06_Org_ECMf" = "#FFC44CFF",  
  "P24-06_Org_Mat" = "#E762D7FF",   
  "P24-12_Org_Mat" = "#686868",
  "P24-14_PAR_Tum" = "#00008BFF", 
  "P24-14_Org_ECMf" = "#A99F1F")  

# cell_cycle_phase
pal_cell_cycle_phase = c(
  "S" ="#1F77B4FF",
  "G1" = "#F8766D",  
  "G2M" = "#2CA02CFF")

# 
pal_walktrap_cluster = c(
  "1" = "#F38400",
  "2" = "#A1CAF1", 
  "3" = "#BE0032", 
  "4" = "#C2B280", 
  "5" = "#848482", 
  "6" = "#008856", 
  "7" = "#E68FAC", 
  "8" = "#0067A5", 
  "9" = "#F99379", 
  "10" = "#604E97", 
  "11"  =  "#F6A600", 
  "12" = "#B3446C", 
  "13" = "#DCD300", 
  "14" = "#882D17", 
  "15" =  "#8DB600", 
  "16" = "#654522", 
  "17" = "#E25822",
  "18" = "#2B3D26",
  "19" = "#4B0082")

pal_walktrap_cluster_old <- c("1" = "#729FCF", 
                              "2" = "#FF9E4A", 
                              "3" = "#67BF5B", 
                              "4" = "#EC665E", 
                              "5" = "#EC96C9", 
                              "6" = "#A8786E", 
                              "7" = "#AD8ACA", 
                              "8" = "#A1A1A1"
)


pal_line_type = c("Early_Passage_PDOs_Basal/Hillock" ="#3CB22D",
                  "Early_Passage_PDOs_Club" = "#FFCC00",
                  "Early_Passage_PDOs_Transitioning" = "#FF6E00",
                  "Early_Passage_PDOs_Tumor" = "#A60021",
                  "Long_Term_PDOs_Tumor" = "#662700")

pal_culture_condition = c(
  "Matrigel" = "#FF1493",
  "ECMf" = "#008B8B" ,
  "Short_Term_Lines_ECMf" = "#008B8B",
  "Short_Term_Lines_Matrigel" = "#FF1493",
  "Long_Term_Lines_ECMf" = "#008B8B",
  "Long_Term_Lines_Matrigel" ="#FF1493",
  "Huang_2023_Matrigel" = "#FF1493",  
  "Song_2022_Matrigel" = "#FF1493",    
  "McCray_2019_Matrigel" = "#FF1493",
  "early_passage_ECMf" = "#008B8B",
  "early_passage_Matrigel" = "#FF1493",
  "parental_tumor_Tissue" = "#B66DFF",
  "Tissue" = "#B66DFF")   

pal_data_type_details = c(
  "Tissue" = "#B66DFF",
  "PDOs" = "#007395",
  "PDOs_Early_Passage_Matrigel" = "#E5B17EFF",
  "PDOs_Long_Term_ECMf" = "#662700",
  "PDOs_Long_Term_Matrigel" = "#662700",
  "PDOs_Early_Passage_ECMf" = "#E5B17EFF",
  "Long_Term" = "#662700",
  "Early_Passage" = "#E5B17EFF",
  "Early_Passage_" = "#E5B17EFF"
)

pal_org_type= c(
  "Early_Passage_PDOs" = "#E5B17EFF",
  "Long_Term_PDOs" = "#662700",
  "Early_Passage_Previously_Published" = "#B66DFF")  

pal_sub = c(
  "Early_Passage_PDOs" = "#E5B17EFF",
  "Long_Term_PDOs" = "#662700",
  "Early_Passage_Previously_Published" = "#B66DFF",
  "Huang_2023" = "#FF7F0EFF",  
  "Song_2022" = "#007395",
  "McCray_2019" = "#EAE32C") 


pal_dataset = c(
  "Gao_Tang" = "purple" , 
  "Puca" = "cyan",
  "Chen_2022" = "#007395",
  "Hirz_2023" = "#EAE32C",
  "Dolgos_2024" = "#3CB22D",
  "Huang_2023" = "#B254A5",  
  "Song_2022" = "#4B0082",
  "Song_2022_Tumor" = "#B254A5",
  "Song_2022_Normal" = "#FFB6DB",
  "McCray_2019" = "#FFB6DB",
  "Chen_2021" = "red",
  "Dolgos_2024_Early_ECMf" = "#E5B17EFF",
  "Dolgos_2024_Early_Matrigel" = "#E5B17EFF",
  "Dolgos_2024_Long_Term" = "#662700",
  "Henry_2018" = "#00008BFF",
  "Tabula_Sapiens" = "#604E97",
  "Dolgos_2024_Early_Passage_Matrigel" = "#FF1493",
  "Dolgos_2024_Early_Passage_ECMf" = "#008B8B"
)  

pal_cell_type <- c(
  "Matrigel_PDOs" = "#FF1493",
  "ECMf_PDOs" = "#008B8B",
  "Basal" = "#3CB22D",
  "Club" = "#FFCC00",
  "Luminal" = "#BE0032",
  "Hillock" = "#358747",
  "Parental_Tumor_cells" = "red",
  "Other Epithelial" = "grey"
  
)

pal_dataset_cell_type <- c(
  "Dolgos_2024_Early_Passage_Matrigel_PDOs" = "#FF1493",
  "Dolgos_2024_Early_Passage_ECMf_PDOs" = "#008B8B",
  "Henry_2018_Basal" = "#3CB22D",
  "Henry_2018_Club" = "#FFCC00",
  "Henry_2018_Luminal" = "#BE0032",
  "Henry_2018_Hillock" = "#358747",
  "Tabula_Sapiens_Basal" = "#3CB22D",
  "Tabula_Sapiens_Other Epithelial" = "grey",  # Default for unmatched type
  "Tabula_Sapiens_Club" = "#FFCC00",
  "Tabula_Sapiens_Luminal" = "#BE0032",
  "Tabula_Sapiens_Hillock" = "#358747",
  "Dolgos_Matrigel_Basal/Hillock" = "#3CB22D",      # Basal/Hillock -> Basal
  "Dolgos_Matrigel_Club" = "#FFCC00",               # Club -> Club
  "Dolgos_Matrigel_Transitioning" = "#FF6E00",         # Transitioning -> Unspecified
  "Dolgos_ECMf_Basal/Hillock" = "#3CB22D",          # Basal/Hillock -> Basal
  "Dolgos_ECMf_Club" = "#FFCC00",                   # Club -> Club
  "Dolgos_ECMf_Tumor" = "#FF1493",                  # Tumor -> Tumor
  "Dolgos_ECMf_Transitioning" = "#FF6E00",             # Transitioning -> Unspecified
  "Dolgos_Matrigel_Tumor" = "#FF1493",              # Tumor -> Tumor
  "Song_Matrigel_Transitioning" = "#FF6E00",           # Transitioning -> Unspecified
  "Song_Matrigel_Basal/Hillock" = "#3CB22D",        # Basal/Hillock -> Basal
  "Song_Matrigel_Other" = "grey",                   # Other -> Unspecified
  "Song_Matrigel_Club" = "#FFCC00",                 # Club -> Club
  "Song_Matrigel_Tumor" = "#FF1493",                # Tumor -> Tumor
  "Song_Parental_Luminal" = "#BE0032",              # Luminal -> Luminal
  "Song_Parental_Basal" = "#3CB22D",                # Basal -> Basal
  "Song_Parental_ERGpos_Tumor" = "#FF1493",         # ERGpos Tumor -> Tumor
  "Song_Parental_ERGneg_Tumor" = "#FF1493",         # ERGneg Tumor -> Tumor
  "Song_Parental_Club" = "#FFCC00"                  # Club -> Club                   # "Club" matches with Club
)

pal_cell_type = c(
  "Basal/Hillock" = "#3CB22D",
  "Club" = "#FFCC00",
  "Transitioning" = "#FF6E00",
  "Tumor" = "#A60020",
  "Parental_Tumor" = "red",
  "Immune cells" = "#9467BDFF",
  "Fibroblasts" = "#1F77B4FF",
  "Endothelial cells" = "#C2B280",
  "Keratinocytes" = "#662700" ,
  "Adipocytes" = "#FFB6DB",
  "Other" = "grey")

pal_tissue_grade <- c(
  "Adjacent_Benign" = "#008600", 
  "ISUP_2" = "#86FF86",         
  "ISUP_3" = "#FFBCFF",         
  "ISUP_4" = "#FF51FF",         
  "ISUP_5" = "#860086",         
  "Not_Applicable" = "#924900FF", 
  "Metastasis" = "#924900FF", 
  "Not_Available" = "#E2E2E2")

pal_tissue_source <- c(
  "Radical_Prostatectomy" = "#75D3FF",   # Rouge vif
  "Metastasis_Resection" = "#99540FFF",    # Orange foncé
  "Prostate_Biopsy" = "#FF3D3D",         # Vert
  "Transurethral_Resection" = "#422CB2FF", # Bleu foncé
  "Not_Available" = "#E2E2E2"            # Violet foncé
)  

pal_cnv_treshold = c("CNV_pos" = "#F69541", "CNV_neg" = "#72A6CE")

pal_tissue_site <- c(
  "Prostate" = "#1F77B4", # Blue
  "Bone" = "#FF7F0E",     # Orange
  "Lung" = "#2CA02C"      # Green
)

pal_tissue_type = c(
  "Primary" = "#75D3FF",
  "Metastasis" = "#99540FFF"
)

pal_previous_treatment = c(
  "Untreated" = "#FFBCFF",
  "Pretreated" ="#860086",
  "Not_Available" = "#E2E2E2"
)

pal_org_origin = c(
  "Patient_Derived" = "#C2B280",
  "PDX_Derived" = "#99540FFF")

pal_zscore <- c(
  "#26456EFF", "#244C7CFF", "#21538BFF", "#1C5A99FF", "#1C63A1FF", "#1C6CAAFF", "#1F74B1FF", 
  "#2B7BB4FF", "#3482B6FF", "#3F8BBAFF", "#4F98C4FF", "#5EA5CEFF", "#78B1D3FF",
  "#9CBBCFFF", "#FFFFFF", "#D6BFBBFF", "#E9A79DFF", "#F78E80FF", "#F6796AFF", 
  "#EC6857FF", "#E25644FF", "#DC4636FF", "#D73529FF", "#D21E1CFF", "#CB1618FF", 
  "#C51517FF", "#BE1316FF", "#B3101BFF", "#A70C20FF", "#9C0824FF"
)

pal_walktrap_condition <- c(
  "Cluster-5-ECMf" = "#68D65E",
  "Cluster-5-Matrigel" = "#3CB22D",
  "Cluster-7-ECMf" = "#68D65E",
  "Cluster-7-Matrigel" = "#3CB22D",
  "Cluster-6-ECMf" = "#68D65E",
  "Cluster-6-Matrigel" = "#3CB22D",
  "Cluster-2-ECMf" = "#68D65E",
  "Cluster-2-Matrigel" = "#3CB22D",
  
  "Cluster-1-ECMf" = "#E06677",
  "Cluster-1-Matrigel" = "#A60021",
  "Cluster-8-ECMf" = "#E06677",
  "Cluster-8-Matrigel" = "#A60021",
  
  "Cluster-3-ECMf" = "#FF9540",
  "Cluster-3-Matrigel" = "#FF6E00",
  
  "Cluster-4-ECMf" = "#FFE066",
  "Cluster-4-Matrigel" = "#FFCC00"
)
