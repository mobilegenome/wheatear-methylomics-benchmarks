
methylBenchmark <- new.env()


methylBenchmark$nested_methylKit_to_df <- function(mk_obj, sample_names, subsampling_frac_list) {
    df_frac <- imap(sample_names, function (samplename, sample_idx) {
            map(
                paste("frac", subsampling_frac_list, sep="-"), 
                function (frac) {
                    methylKit::getData(mk_obj[[frac]][[sample_idx]]) |>
                    #head(100) |>
                    dplyr::rename("seqnames" = "chr") |>
                    dplyr::mutate(
                        sample = samplename,
                        subset_frac = str_remove(frac, "`") 
                    ) 
                }
            ) |> 
            compact() |>
            list_rbind() 
        }
    ) |> 
    list_rbind()

    return(df_frac) # Ensure the function returns the result
}

methylBenchmark$get_relative_error <- function(data_wide, x, y) {
    
    data_wide |>
        left_join(
            wheatearcommons::df_genome |> 
                dplyr::filter(genome == "OenMel1.1") |> 
                select(seqnames, seqlength)
        ) |>
        select(seqnames, seqlength, start, sample, !!sym(x), !!sym(y)) |>
        group_by(seqnames, sample, seqlength) |>
        mutate(diff = !!sym(x) - !!sym(y)) |>
        mutate(is_FN = if_else(!!sym(x) == 0 & !!sym(y) > 0, 1, 0))
}

methylBenchmark$get_correlation <- function(data_wide, x_lab, y_lab) {

    # calculate correlation coefficient and p-value between the two fractions
    cor.test(
        data_wide[[x_lab]],
        data_wide[[y_lab]]) |>
        broom::tidy() |>
        dplyr::select(estimate, p.value) 
}

methylBenchmark$widen_data <- function(data) {
    data_wide <-
        data |>
        # calculate mCpG
        mutate(mCpG = numCs / coverage) |>
        select(seqnames, start, coverage, mCpG, sample, subset_frac) |>
        pivot_wider(id_cols = c(seqnames, start, sample), names_from = subset_frac, values_from = mCpG)
    return(data_wide)

}

#' Plots the fraction correlation between two variables for a specific sample and sequence ----
#'
#' Assumes chromosome/contig in column "seqnames"
#' 
#' @param .data The data frame containing the relevant data.
#' @param samplename The name of the sample to filter the data for.
#' @param seqname The name of the sequence to filter the data for.
#' @param x The name of the variable to plot on the x-axis.
#' @param y The name of the variable to plot on the y-axis.
methylBenchmark$plot_frac_correlation <- function(data_wide,
    samplename,
    seqname,
    x,
    y) {

    p <- ggplot(data_wide,
        aes(
            x=!!sym(x),
            y=!!sym(y))
        ) + 
        geom_point() + 
        geom_smooth(method="lm") + 
        ggtitle(glue("{samplename} - {y} vs {x} - {seqname}"))

        gc() 

        return(p)
}


# get list of sample names
methylBenchmark$get_samplenames_from_mkobj <- function(mk_obj) {
    sample_names <-
        map(mk_obj, ~ methylKit::getSampleID(.x)) |> 
        unique() |>
        unlist()
}


#' data in long format
methylBenchmark$calculate_conservation_score <- function(data) {
    data |>
        mutate(mCpG = numCs/coverage) |>
        complete(sample, subset_frac, seqnames, start, fill = list(mCpG = NA)) |> # Ensure all sample-site combinations exist
        mutate(call_present = !is.na(mCpG)) |>
        group_by(seqnames, start, subset_frac) |>
        summarize(conservation_score = mean(call_present)) 
}

methylBenchmark$plot_boxplot_conservation <- function(df_cons) {
    ggplot(df_cons, aes(
        x=subset_frac,
        y=conservation_score*6,
        fill=subset_frac)) + 
        geom_violin() + 
        geom_boxplot(fill="white", alpha=0.5, width=0.2) + 
        theme_minimal() + 
        theme(
            axis.text.x = element_blank(),
            axis.title.x = element_blank(),
            axis.ticks.x = element_blank()
        )
}


methylBenchmark$calculate_table_of_sites_present <- function(df_cons) {
    df_cons |>
        group_by(subset_frac) |>
        # bin conservation scores in as many bins as there are samples, add 0.01 to achieve this
        mutate(conservation_bin = cut_width(conservation_score, width = 1/6+0.01, boundary = 0, labels=FALSE)) |> 
        group_by(subset_frac, conservation_bin) |>
        # calculate percentage per subset_frac
        summarize(n = n()) |>
        pivot_wider(names_from = conservation_bin, values_from = n) |>
        janitor::adorn_totals("col") |>
        knitr::kable()
}


# calculate number of missing for all pairwise comparisons between fractions
# and calculate the percentage of missing values
methylBenchmark$print_table_of_missing_sites <- function(data_wide) {
    data_wide |> 
        mutate(
            missing = case_when(
                is.na(`frac-0.125`) & !is.na(`frac-0.5`) ~ "comp-0.125-0.5",
                is.na(`frac-0.125`) & !is.na(`frac-0.25`) ~ "comp-0.125-0.25",
                is.na(`frac-0.25`) & !is.na(`frac-0.5`) ~ "comp-0.25-0.5",
                TRUE ~ "no-missing")) |>
        group_by(missing) |>
        summarize(n = n()) |>
        mutate(perc = round(100*(n/sum(n)),2)) |> 
        janitor::adorn_totals("row") |>
        knitr::kable()
}


methylBenchmark$plot_mean_mCpG_by_bin_along_chr <- function(data) {
    p <- data |> 
        arrange(start) |>
        mutate(mCpG = numCs/coverage) |>
        group_by(pos = consecutive_id(seqnames, start)) |>
        select(-c(end, strand)) |>
        ungroup() |>
        mutate(bin = cut_width(pos, width=1e3, boundary = 0, labels=FALSE)) |> 
        group_by(subset_frac, bin) |>
        mutate(
            mean_mCpG = mean(mCpG, na.rm=F),
            sd_mCpG = sd(mCpG)
        ) |>
        ggplot(aes(x=as.factor(bin), y=mean_mCpG, color=subset_frac)) + 
        geom_line(aes(group=subset_frac)) + 
        theme_minimal() + 
        theme(
            axis.text.x = element_blank(),
            axis.title.x = element_blank(),
            axis.ticks.x = element_blank()
        )
    
    return(p)
}

methylBenchmark$calculate_methylation_bins <- function(data) {
    data |> 
        select(-c(end, strand)) |>
        arrange(start) |>
        mutate(mCpG = numCs/coverage) |>
        group_by(pos = consecutive_id(seqnames, start)) |>
        ungroup() |>
        pivot_wider(id_cols = c(sample, seqnames, start),
            names_from = subset_frac,
            values_from = mCpG) |> 
        mutate(bin = cut_width(`frac-0.5`, width=0.05, boundary = 0, labels=FALSE)) 
}

methylBenchmark$calculate_na_counts <- function(data) {
    data |>  
        group_by(sample) |>
        summarise(
            n = n(),
            n_frac125 = sum(!is.na(`frac-0.125`)),
            n_frac25 = sum(!is.na(`frac-0.25`)),
            n_frac5 = sum(!is.na(`frac-0.5`)),
            # NA counts
            na_frac125 = sum(is.na(`frac-0.125`)),
            na_frac25 = sum(is.na(`frac-0.25`)),
            na_frac5 = sum(is.na(`frac-0.5`)),
            # NA fractions
            perc_na_frac125 = round(100*na_frac125/n, 2),
            perc_na_frac25 = round(100*na_frac25/n, 2),
            perc_na_frac5 = round(100*na_frac5/n, 2)
        ) 
}

methylBenchmark$plot_boxplot_with_na_counts <- function(data, na_count) {

    data_median_x_axis <- data |>
        group_by(sample, bin) |>
        summarize(median = round(median(`frac-0.5`, na.rm=T),2)) |> 
        drop_na() |>  
        select(bin,sample,  median)

    data_long <- 
        data |>
        select(-`frac-0.5`) |>
        drop_na() |>
        # pivot_long again 
        pivot_longer(
            cols = contains("frac"),
            names_to = "subset_frac",
            values_to = "mCpG"
        ) 

    ggplot() + 
    geom_boxplot(
        data = data_long,
        aes(x=as.factor(bin), y=mCpG, fill=subset_frac)
    ) + 
    scale_fill_manual(values=
        c("frac-0.125" = "lightgray",
            "frac-0.25" = "gray")) +
    geom_point(
        data = data_median_x_axis, 
        aes(x=bin, y=median), shape = 4, size = 3, color="red") +
        facet_wrap(~sample) +
    ggtitle("", glue("dropped {na_count} % of sites not present at lower coverage"))
}

# Plot by position
methylBenchmark$plot_frac_correlation_by_pos <- function(
    data_long,
    samplename,
    seqname,
    x,
    y,
    filename = "results/{samplename}-fraccor_pos-{seqname}-{fname_part}.png"
    ) {
    
    seqlength <- wheatearcommons::df_genome |> 
        dplyr::filter(genome == "OenMel1.1", 
        seqnames == seqname) |> 
        dplyr::pull(seqlength)



    p <- data_long |>
        dplyr::filter(sample == samplename) |>
        dplyr::filter(str_starts(seqnames, "chr")) |>
        left_join(
            wheatearcommons::df_genome |> 
                dplyr::filter(genome == "OenMel1.1") |> 
                select(seqnames, seqlength)
        ) |>
        #dplyr::filter(str_detect(chr, seqname)) |>
        # calculate mCpG
        mutate(
            mCpG = numCs / coverage,
            #start = round(start/seqlength*100),5)
            ) |> 
        select(seqnames,seqlength, start, coverage, mCpG, sample, subset_frac) |>
        filter(subset_frac %in% c(x, y)) |>
        group_by(seqnames, subset_frac) |>
        mutate(bin = cut_width(start, width=1e5, labels=FALSE)) |> 
        group_by(bin,seqnames,seqlength, subset_frac) |>
        summarize(mCpG = mean(mCpG)) |> 
        pivot_wider(id_cols = c(seqnames, bin, seqlength), names_from = subset_frac, values_from = mCpG) |>
        mutate(rel_error = (!!sym(x) - !!sym(y)) / !!sym(x)) |>
        # caluclate difference between fractions 
        ggplot(aes(
            x=bin,
            y=reorder(seqnames, seqlength),
            fill=rel_error)
        ) + 
        colorspace::scale_fill_continuous_diverging() + 
        geom_tile() + 
        ggtitle(samplename) 

        fname_part <- paste(str_extract(x, "[0-9\\.]+"), "-", str_extract(y, "[0-9\\.]+"), sep="")

        return(p)
}
