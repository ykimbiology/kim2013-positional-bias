

plot_distr <- function(normpos_list, main_text) {
    hist(normpos_list, n=5, main=main_text, xlab='normpos', ylab='count')
}

get_fraction_position_bp_b <- function(pos_start, pos_end, seq_length) {
    #'''Second version suggested by bpeters.
    #   pos_start uses zero-index system.'''
    peptide_length = pos_end - pos_start + 1
    pos_fraction = pos_start/(seq_length - peptide_length + 1)
    
        #print 'error pos fraction', pos_start, pos_end, seq_length, pos_fraction
    pos_fraction
}

sample_normpos <- function(dtab) {
    # [pos_norm, peptide_len, protein_len]
    kmax <- 10
    normpos_list <- rep(-1,length(dtab$normpos)*kmax)
    index <- 1
    for (k in 1:kmax) {
    for (i in 1:length(dtab$normpos)) {
        protein_len <- dtab$protein_len[3]
        peptide_len <- dtab$peptide_len[2]
        protein_len_sub <- protein_len - peptide_len + 1  # for 9mer 9 - 9 + 1 = 1
        pos_start <- sample.int(protein_len_sub, size=1, replace=T)
        pos_end <- pos_start + peptide_len - 1
        normpos <- get_fraction_position_bp_b(pos_start-1, pos_end-1, protein_len)
        normpos_list[index] <- normpos
        index <- index + 1
    }
    }
    normpos_list
}

dtab_a <- read.table('normpos.bruteforce.all.txt',header=T)
dtab_b <- read.table('normpos.blast_a.all.txt',header=T)
normpos_list_theoretical <- sample_normpos(dtab_a)

png(filename='plot_distr_normpos.png', width=1200, height=600)
par(mfrow=c(1,3), pty='s', cex=1.5)
plot_distr(dtab_a$normpos, 'bruteforce')
plot_distr(dtab_b$normpos, 'blast_a')
plot_distr(normpos_list_theoretical, 'theoretical')
dev.off()
