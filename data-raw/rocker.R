rocker = create_protein_dat(path_u = "inst/extdata/Rocker_ref_sequences_filtered.txt.gz",
                            path_l = "inst/extdata/Rocker_sel_sequences_filtered.txt.gz")

usethis::use_data(rocker,overwrite = T)
