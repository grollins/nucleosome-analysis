library(tidyverse)

df <- read_delim("matrix_nuc_signal_spcyt.txt", "\t", escape_double = FALSE,
                trim_ws = TRUE, skip = 2)

df$idx <- seq(1, nrow(df))
df$gene_group <- ""
df$gene_group[1:10] <- "bam.off_48.off_72.on"
df$gene_group[11:459] <- "bam.off_48.on_72.on"
df$gene_group[460:479] <- "bam.on_48.on_72.up"
df$gene_group[480:1680] <- "bam.on_48_72.up"

table(df$gene_group)

df2 <- df %>%
    gather(sample_name_raw, value, -idx, -gene_group)

parse_sample_name <- function(x) {
    return(stringr::str_split(x, '\\.', simplify = T)[[1]])
}

parse_sample_num <- function(x) {
    stringr::str_split(x, '_') %>%
        unlist() %>% 
        tail(n = 1) %>%
        as.integer()
}

df2$sample_name <- map_chr(df2$sample_name_raw, parse_sample_name)
df2$sample_num <- map_chr(df2$sample_name_raw, parse_sample_num)

df3 <- df2 %>%
    mutate(sample_num = ifelse(is.na(sample_num), 0, sample_num)) %>%
    mutate(sample_num = as.integer(sample_num + 1))

df3 %>%
    filter(gene_group == 'bam.on_48_72.up') %>%
    # filter(sample_name == 'aly') %>%
    mutate(value2 = ifelse(value > 2, 2, value)) %>%
    mutate(value2 = ifelse(value2 < -2, -2, value2)) %>%
    ggplot(aes(sample_num, idx)) +
        geom_raster(aes(fill = value2)) +
        facet_wrap(~sample_name) +
        scale_y_reverse()

# df3 %>% sample_frac(0.1) %$% sample_num %>% table(useNA = 'ifany') %>% prop.table()
# df3 %>% sample_frac(0.5) %>% filter(is.na(value)) %>% nrow()

df3_mean <- df3 %>%
    filter(idx < 10) %>%
    filter(sample_name == 'can') %>%
    group_by(idx) %>%
    summarize(mu = mean(value))

df3 %>%
    filter(idx < 10) %>%
    filter(sample_name == 'can') %>%
    left_join(df3_mean, by = 'idx') %>%
    ggplot(aes(sample_num, value)) +
        geom_line() +
        geom_hline(aes(yintercept = mu), size = 0.25) +
        facet_wrap(~idx)

find_peaks <- function(df) {
    pp_matrix <- df %>%
        select(idx, value, sample_num) %>%
        spread(idx, value) %>%
        select(-sample_num) %>%
        as.matrix() %>%
        peakPick::detect.spikes(roi = c(20, 50), winlen = 10, spike.min.sd = 6)
    
    pp_df <- pp_matrix %>%
        as_tibble() %>%
        mutate(sample_num = seq(1, nrow(.))) %>%
        gather(idx, is_peak, -sample_num) %>%
        mutate(idx = as.integer(idx))
    
    df_with_peaks <- df %>%
        left_join(pp_df, by = c('idx', 'sample_num'))
    
    return(df_with_peaks)
}

df4 <- df3 %>%
    filter(idx < 10) %>%
    split(.$sample_name) %>%
    map(find_peaks) %>%
    bind_rows()


df4 %>%
    ggplot(aes(sample_num, value, color = sample_name)) +
        geom_line() +
        facet_wrap(~idx, scales = 'free_y') +
        geom_point(data = filter(df4, is_peak == T), aes(sample_num, value,
                   color = sample_name)) +
        ggthemes::theme_few()

df4 %>%
    filter(idx == 1) %>%
    filter(sample_name == 'bam') %>%
    filter(is_peak == T) %$% value %>%
    max()

df_max_peak_values <- df4 %>%
    filter(is_peak == T) %>%
    group_by(idx, sample_name) %>%
    summarize(m = max(value))

df5 <- df4 %>%
    left_join(df_max_peak_values, by = c('idx', 'sample_name')) %>%
    mutate(is_peak2 = ifelse(abs(m - value) < 1e-5, T, F))

df5 %>%
    ggplot(aes(sample_num, value, color = sample_name)) +
        geom_line() +
        facet_wrap(~idx, scales = 'free_y') +
        geom_point(data = filter(df5, is_peak == T), aes(sample_num, value,
                                                     color = sample_name)) +
        geom_point(data = filter(df5, is_peak2 == T), aes(sample_num, value),
                                                      color = 'black') +
        ggthemes::theme_few()

df6 <- df5 %>%
    split(.$idx) %>%
    map(function(x) {
        bam <- x %>%
            filter(sample_name == 'bam') %>%
            filter(is_peak2 == T)
        
        bam_idx = ifelse(nrow(bam) == 0, 0, bam$sample_num)
        
        y <- x %>%
            mutate(delta = sample_num - bam_idx)
        return(y)
    }) %>%
    bind_rows()

df6 %>%
    filter(is_peak2 == T) %>%
    filter(sample_name != 'bam') %>%
    mutate(idx_fct = as.factor(idx)) %>%
    select(idx_fct, sample_name, value, delta) %>%
    ggplot(aes(value, delta, color = sample_name)) +
        geom_point() +
        facet_wrap(~idx_fct) +
        ggthemes::theme_few()