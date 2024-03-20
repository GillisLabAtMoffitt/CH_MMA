# MMA in all tumor types

# Import library
library(tidyverse)

################################################################################################# Load data
# path1 <- fs::path("","Volumes","Gillis_Research","Christelle Colin-Leitzinger", "CHIP in Avatar",
#                   "MM avatar", "Jamie")
# Clinical_linkage <- read.delim(paste0(path1, "/wes_somatic_mutations_metadata_v0.4.5.txt")) %>% 
#   select(c("subject", SLID_germline, SLID_tumor, moffittSampleId, 
#            moffittSampleId_tumor, moffittSampleId_germline,
#            "ClinicalSpecimenLinkage_AgeAtSpecimenCollection",
#            ClinicalSpecimenLinkage_HistologyBehavior, SpecimenDetail_DiseaseType)) %>% 
#   arrange(subject, ClinicalSpecimenLinkage_AgeAtSpecimenCollection)
# Samples dates and Ids for v04.5
path2 <- fs::path("","Volumes","Gillis_Research","Christelle Colin-Leitzinger", "CHIP in Avatar",
                  "CH all tumor types", "raw data", "M2GEN")
sample_data <- 
  readxl::read_xlsx(paste0(path2, "/Yifen data/Avatar SLIDs and 06Sdatesv2.xlsx"),
                    na = "NULL") %>% 
  select(-DateOfCollection,
         tumor_anatomic_site = "tumor anatomic site", germline_anatomic_site = "germline anatomic site")

path <- fs::path("","Volumes","Gillis_Research","Christelle Colin-Leitzinger", "CHIP in Avatar",
                  "CH all tumor types")
WES_v4.7 <- read_csv(paste0(path, "/processed WES ids list/matched germline tumor samples ids all tumor type v04.7.csv"))

IDs_list <- readxl::read_xlsx(here::here("MMA_avatar_lab_report.xlsx"), 
                              sheet = "IDs list") %>%
  mutate(PATIENT_ID = as.character(PATIENT_ID)) %>%
  mutate(MRN = as.character(MRN))
MMA_2 <- readxl::read_xlsx(here::here("MMA_avatar_lab_report.xlsx"), 
                           sheet = "Lab_MMA") %>%
  left_join(IDs_list, .,
            by = "PATIENT_ID") %>%
  filter(!is.na(LAB_RESULT))


################################################################################################# Data mining

WES_v4.5 <- inner_join(sample_data, Clinical_linkage,
                       by = c("subject", "SLID_germline", "SLID_tumor", 
                              "moffittSampleId", "moffittSampleId_tumor", 
                              "moffittSampleId_germline")) %>% 
  select(-moffittSampleId, -ClinicalSpecimenLinkage_AgeAtSpecimenCollection) %>% 
  rename(avatar_key = subject,
         tumor_site_collection = tumor_anatomic_site,
         histology = ClinicalSpecimenLinkage_HistologyBehavior) %>% 
  filter(germline_anatomic_site == "Blood")

WES_v4.7 <- WES_v4.7 %>% 
  select(-germline_sample_family_id, dob, date_of_diagnosis, age_at_diagnosis) %>% 
  rename(SLID_tumor = dna_sequencing_library_id,
         moffittSampleId_germline = "germline_orien_id",
         DateOfCollection_germline = "germline_collection_date",
         moffittSampleId_tumor = "tumor_orien_id", 
         SLID_RNA_tumor = "rna_sequencing_library_id",
         DateOfCollection_tumor = "tumor_collection_date")

WES_v4.7 <- WES_v4.7 %>% 
  mutate(interval_germ_tumor = interval(start = DateOfCollection_germline, 
                                        end = DateOfCollection_tumor)/
           duration(n=1, units = "days"))

WES_all_version <- bind_rows(WES_v4.5, 
                             WES_v4.7,
                             .id = "version")

Germline_available <- WES_all_version %>% 
  filter(!is.na(SLID_germline))

MMA_in_germline <- inner_join(Germline_available %>% 
                                mutate(mrn = as.character(mrn)), 
                              MMA_2,
                              by= c("mrn" = "MRN"))

non_HEM <- # MMA_in_germline %>% filter(histology != "97323 Multiple myeloma") %>% 
  # distinct(avatar_key, SLID_germline, .keep_all = TRUE)
  read_csv(paste0(path, "/processed WES ids list/matched germline tumor samples ids non-HEM v04.7.csv")) %>% 
  mutate(mrn = as.character(mrn)) %>% 
  inner_join(. , 
             MMA_2 %>% 
               distinct(MRN, .keep_all = TRUE),
             by= c("mrn" = "MRN"))
MM_MMA <- # MMA_in_germline %>% filter(histology == "97323 Multiple myeloma") %>% 
  # distinct(avatar_key, SLID_germline, .keep_all = TRUE)
  read_csv(paste0(path, "/processed WES ids list/matched germline tumor samples ids Multiple Myeloma v04.7.csv")) %>% 
  mutate(mrn = as.character(mrn)) %>% 
  inner_join(. , 
             MMA_2 %>% 
               distinct(MRN, .keep_all = TRUE),
             by= c("mrn" = "MRN"))
MGUS_MMA <- # MMA_in_germline %>% filter(histology == "97323 Multiple myeloma") %>% 
  # distinct(avatar_key, SLID_germline, .keep_all = TRUE)
  read_csv(paste0(path, "/processed WES ids list/matched germline tumor samples ids MGUS v04.7.csv")) %>% 
  mutate(mrn = as.character(mrn)) %>% 
  inner_join(. , 
             MMA_2 %>% 
               distinct(MRN, .keep_all = TRUE),
             by= c("mrn" = "MRN"))


write_csv(non_HEM %>% select("avatar_key" : interval_germ_tumor), 
          "sample ids non-HEM v407 with MMA data.csv")
write_csv(MGUS_MMA %>% select("avatar_key" : interval_germ_tumor), 
          "sample ids MGUS v407 with MMA data.csv")
write_csv(MM_MMA %>% select("avatar_key" : interval_germ_tumor), 
          "sample ids in multiple myeloma patients v407 with MMA data.csv")



write_csv(sample_data %>% 
            distinct(subject), 
          "ids in v04.5 all tumor types data.csv")


