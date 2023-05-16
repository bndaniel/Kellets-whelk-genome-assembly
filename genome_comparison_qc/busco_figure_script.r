library(ggplot2);
library(ggtext);
# library(gridExtra);
library(scales);
library(rjson);

# quast and busco analysis; CA
figure.data <- data.frame();

# extract data from quast
setwd('~/Desktop/whelk_busco_figure/quast_summaries');
for (quast.report in dir()){
    quast.data <- read.delim(
        file = quast.report, 
        sep = '\t', 
        header = FALSE,
        col.names = c('statistic', 'value')
        );
    assembly <- paste(
            strsplit(
                quast.data[quast.data$statistic == "Assembly", 2], 
                split = ".", 
                fixed = TRUE)[[1]][c(1, 2)], 
            collapse = '.'
            );
    assembly <- paste(
            strsplit(assembly, split = "_", fixed = TRUE)[[1]][c(1, 2)], 
            collapse = '_'
            );
    assembly <- gsub(pattern = '_NA', replacement = '', x = assembly);
    total.length <- quast.data[quast.data$statistic == "Total length (>= 0 bp)", 2];
    num.contigs <- quast.data[quast.data$statistic == "# contigs (>= 0 bp)", 2];
    largest.contig <- quast.data[quast.data$statistic == "Largest contig", 2];
    GC.content <- quast.data[quast.data$statistic == "GC (%)", 2];
    N50 <- quast.data[quast.data$statistic == "N50", 2];
    L50 <- quast.data[quast.data$statistic == "L50", 2];
    figure.data <- rbind(
        figure.data, 
        c(
        assembly, 
        total.length, 
        num.contigs,
        largest.contig, 
        GC.content, 
        N50, 
        L50
        )
    );
};
colnames(figure.data) <- c(
    "assembly", 
    "total.length", 
    "num.contigs",
    "largest.contig",
    "GC.content",
    "N50",
    "L50"
    );


# extract data from busco
readShortSummary <- function(figure.df, directory, dataset) {
    setwd(directory);
    figure.df$busco.data = NA;
    for (short.summary in dir()){
        file_data <- fromJSON(
            file = short.summary
            );
        completeness <- round(file_data$C / as.numeric(file_data$dataset_total_buscos) * 100, digits = 1);
        accession <- paste(
            strsplit(short.summary, split = ".", fixed = TRUE)[[1]][c(4, 5)], 
            collapse = '.'
            );
        accession <- gsub(pattern = '_busco_output', replacement = '', x = accession); # had to do this for the one exception, assembly.fasta
        accession <- paste(
            strsplit(accession, split = "_", fixed = TRUE)[[1]][c(1, 2)], 
            collapse = '_'
            );
        accession <- gsub(pattern = '_NA', replacement = '', x = accession);
        figure.df$busco.data[figure.df$assembly == accession] <- completeness;
        } ;
    names(figure.df)[names(figure.df) == "busco.data"] <- dataset;
    return(figure.df);
    };

figure.data = readShortSummary(
    figure.df = figure.data, 
    directory = '~/Desktop/whelk_busco_figure/busco_short_summaries', 
    dataset = "metazoaBUSCO"
    );

# color coding for the plot
taxonomy.group <- c(
    'GCA_011634625.1' = 'Buccindae',
    'GCA_018398815.1' = 'Negastropoda',
    'GCA_017654935.1' = 'Negastropoda',
    'GCA_009936545.1' = 'Negastropoda',
    'GCA_016801955.1' = 'Negastropoda',
    'GCA_001262575.1' = 'Negastropoda',
    'GCA_004193615.1' = 'Negastropoda',
    'GCA_018857735.1' = 'Ceanogastropoda',
    'GCA_018292915.1' = 'Ceanogastropoda',
    'GCA_004794575.1' = 'Ceanogastropoda',
    'GCA_004794325.1' = 'Ceanogastropoda',
    'GCA_004794655.1' = 'Ceanogastropoda',
    'GCF_003073045.1' = 'Ceanogastropoda',
    'GCA_944989445.1' = 'Gastropoda',
    'GCA_936450465.1' = 'Gastropoda',
    'GCA_932274485.1' = 'Gastropoda',
    'GCA_022045235.1' = 'Gastropoda',
    'GCA_012295275.1' = 'Gastropoda',
    'GCF_000327385.1' = 'Gastropoda',
    'GCF_000002075.1' = 'Gastropoda',
    'GCA_019648995.1' = 'Gastropoda',
    'GCA_025434175.1' = 'Gastropoda',
    'GCF_016097555.1' = 'Gastropoda',
    'assembly.fasta' = 'Flye',
    'pb1.dmo' = 'SMARTdenovo',
    'assembly.scaffolds' = 'MaSuRCA'
    );

figure.data$group <- taxonomy.group[match(
    x = figure.data$assembly,
    table = names(taxonomy.group)
    )];

figure.data$group <- factor(figure.data$group, levels = c("Flye", "MaSuRCA", "SMARTdenovo", "Gastropoda", "Negastropoda", "Ceanogastropoda" , "Buccindae"));
figure.data$num.contigs <- as.numeric(figure.data$num.contigs); 
figure.data$N50 <- as.numeric(figure.data$N50); 
figure.data$total.length <- as.numeric(figure.data$total.length); 

figure.data$scaled.contigs <- figure.data$total.length / figure.data$num.contigs;
figure.data$scaled.N50 <- figure.data$total.length / figure.data$N50;

# contigs plot
contig.plot <- ggplot(figure.data) +
  
  geom_point(aes(x = metazoaBUSCO, y = num.contigs, color = group, shape = group), size = 5) +

    scale_color_manual(name = expression(underline("Genomes")), 
                     labels = c(
                        expression(paste(italic("K. kelletii"), "  (Flye)")), 
                        expression(paste(italic("K. kelletii"), "  (MaSuRCA)")), 
                        expression(paste(italic("K. kelletii"), "  (SMARTdenovo)")), 
                        "Gastropoda", "Negastropoda", "Ceanogastropoda" , "Buccindae"), 
                     values = c("blue", "#0096FF", "#88ACE0", "#FFC000", "#380474", "#228B22", 'slategrey')) +
  
  scale_shape_manual(name = expression(underline("Genomes")), 
                     labels = c(
                        expression(paste(italic("K. kelletii"), "  (Flye)")), 
                        expression(paste(italic("K. kelletii"), "  (MaSuRCA)")), 
                        expression(paste(italic("K. kelletii"), "  (SMARTdenovo)")), 
                        "Gastropoda", "Negastropoda", "Ceanogastropoda" , "Buccindae"),
                     values = c(17, 17, 17, 19, 19, 19, 19)) +
  
  scale_y_log10(limits = c(10, 1e+06), breaks = c(1e+01, 1e+02, 1e+03, 1e+04, 1e+05, 1e+06), labels = trans_format(trans = "log10", format = math_format(10^.x))) +
  
  scale_x_continuous(breaks = seq(from = 40, to = 100, by = 10), labels = paste(seq(from = 40, to = 100, by = 10), "%", sep = ""), limits = c(40, 100)) +
  
  labs(x = "BUSCO Score", y = "Scaffolds and Contigs") +
  
  theme_classic() +
  
  theme(legend.box.background = element_rect(color = "black", size = 1), 
        legend.text.align = 0,
        legend.text = element_text(size = 18),
        legend.title = element_text(size = 20),
        axis.text.x = element_text(size = 18, color = "black"), 
        axis.text.y = element_text(size = 18, color = "black"), 
        axis.title.x = element_text(size = 20), 
        axis.title.y = element_text(size = 20))

# scaled contig plot
scaled.contig.plot <- ggplot(figure.data) +
  
  geom_point(aes(x = metazoaBUSCO, y = scaled.contigs, color = group, shape = group), size = 5) +

    scale_color_manual(name = expression(underline("Genomes")), 
                     labels = c(
                        expression(paste(italic("K. kelletii"), "  (Flye)")), 
                        expression(paste(italic("K. kelletii"), "  (MaSuRCA)")), 
                        expression(paste(italic("K. kelletii"), "  (SMARTdenovo)")), 
                        "Gastropoda", "Negastropoda", "Ceanogastropoda" , "Buccindae"), 
                     values = c("blue", "#0096FF", "#88ACE0", "#FFC000", "#380474", "#228B22", 'slategrey')) +
  
  scale_shape_manual(name = expression(underline("Genomes")), 
                     labels = c(
                        expression(paste(italic("K. kelletii"), "  (Flye)")), 
                        expression(paste(italic("K. kelletii"), "  (MaSuRCA)")), 
                        expression(paste(italic("K. kelletii"), "  (SMARTdenovo)")), 
                        "Gastropoda", "Negastropoda", "Ceanogastropoda" , "Buccindae"),
                     values = c(17, 17, 17, 19, 19, 19, 19)) +
  
  scale_y_log10(limits = c(10, 1e+06), breaks = c(1e+01, 1e+02, 1e+03, 1e+04, 1e+05, 1e+06), labels = trans_format(trans = "log10", format = math_format(10^.x))) +
  
  scale_x_continuous(breaks = seq(from = 40, to = 100, by = 10), labels = paste(seq(from = 40, to = 100, by = 10), "%", sep = ""), limits = c(40, 100)) +
  
  labs(x = "BUSCO Score", y = "Genome Size / Scaffolds and Contigs") +
  
  theme_classic() +
  
  theme(legend.box.background = element_rect(color = "black", size = 1), 
        legend.text.align = 0,
        legend.text = element_text(size = 18),
        legend.title = element_text(size = 20),
        axis.text.x = element_text(size = 18, color = "black"), 
        axis.text.y = element_text(size = 18, color = "black"), 
        axis.title.x = element_text(size = 20), 
        axis.title.y = element_text(size = 20))


setwd('/Users/cassidyandrasz/Desktop/whelk_busco_figure');
ggsave(
    filename = "whelk_contig_plot_v5.png",
    plot = contig.plot,
    device = 'png',
    width = 11, 
    height = 7,
    units = 'in'
    );

# N50 plot

N50.plot <- ggplot(figure.data) +
  
  eom_point(aes(x = metazoaBUSCO, y = N50, color = group, shape = group), size = 5) +

scale_color_manual(name = expression(underline("Genomes")), 
                     labels = c(
                        expression(paste(italic("K. kelletii"), "  (Flye)")), 
                        expression(paste(italic("K. kelletii"), "  (MaSuRCA)")), 
                        expression(paste(italic("K. kelletii"), "  (SMARTdenovo)")), 
                        "Gastropoda", "Negastropoda", "Ceanogastropoda" , "Buccindae"), 
                     values = c("blue", "#0096FF", "#88ACE0", "#FFC000", "#380474", "#228B22", 'slategrey')) +
  
  scale_shape_manual(name = expression(underline("Genomes")), 
                     labels = c(
                        expression(paste(italic("K. kelletii"), "  (Flye)")), 
                        expression(paste(italic("K. kelletii"), "  (MaSuRCA)")), 
                        expression(paste(italic("K. kelletii"), "  (SMARTdenovo)")), 
                        "Gastropoda", "Negastropoda", "Ceanogastropoda" , "Buccindae"),
                     values = c(17, 17, 17, 19, 19, 19, 19)) +
  
  scale_y_log10(limits = c(10, 1e+09), breaks = c(1e+01, 1e+02, 1e+03, 1e+04, 1e+05, 1e+06, 1e+07, 1e+08, 1e+09), labels = trans_format(trans = "log10", format = math_format(10^.x))) +
  
  scale_x_continuous(breaks = seq(from = 40, to = 100, by = 10), labels = paste(seq(from = 40, to = 100, by = 10), "%", sep = ""), limits = c(40, 100)) +
  
  labs(x = "BUSCO Score", y = "N50") +
  
  theme_classic() +
  
  theme(legend.box.background = element_rect(color = "black", size = 1), 
        legend.text.align = 0,
        legend.text = element_text(size = 18),
        legend.title = element_text(size = 20),
        axis.text.x = element_text(size = 18, color = "black"), 
        axis.text.y = element_text(size = 18, color = "black"), 
        axis.title.x = element_text(size = 20), 
        axis.title.y = element_text(size = 20))

# scaled N50 plot

scaled.N50.plot <- ggplot(figure.data) +
  
geom_point(aes(x = metazoaBUSCO, y = scaled.N50, color = group, shape = group), size = 5) +

scale_color_manual(name = expression(underline("Genomes")), 
                     labels = c(
                        expression(paste(italic("K. kelletii"), "  (Flye)")), 
                        expression(paste(italic("K. kelletii"), "  (MaSuRCA)")), 
                        expression(paste(italic("K. kelletii"), "  (SMARTdenovo)")), 
                        "Gastropoda", "Negastropoda", "Ceanogastropoda" , "Buccindae"), 
                     values = c("blue", "#0096FF", "#88ACE0", "#FFC000", "#380474", "#228B22", 'slategrey')) +
  
  scale_shape_manual(name = expression(underline("Genomes")), 
                     labels = c(
                        expression(paste(italic("K. kelletii"), "  (Flye)")), 
                        expression(paste(italic("K. kelletii"), "  (MaSuRCA)")), 
                        expression(paste(italic("K. kelletii"), "  (SMARTdenovo)")), 
                        "Gastropoda", "Negastropoda", "Ceanogastropoda" , "Buccindae"),
                     values = c(17, 17, 17, 19, 19, 19, 19)) +
  
  scale_y_log10(limits = c(10, 1e+09), breaks = c(1e+01, 1e+02, 1e+03, 1e+04, 1e+05, 1e+06, 1e+07, 1e+08, 1e+09), labels = trans_format(trans = "log10", format = math_format(10^.x))) +
  
  scale_x_continuous(breaks = seq(from = 40, to = 100, by = 10), labels = paste(seq(from = 40, to = 100, by = 10), "%", sep = ""), limits = c(40, 100)) +
  
  labs(x = "BUSCO Score", y = "Genome Size / N50") +
  
  theme_classic() +
  
  theme(legend.box.background = element_rect(color = "black", size = 1), 
        legend.text.align = 0,
        legend.text = element_text(size = 18),
        legend.title = element_text(size = 20),
        axis.text.x = element_text(size = 18, color = "black"), 
        axis.text.y = element_text(size = 18, color = "black"), 
        axis.title.x = element_text(size = 20), 
        axis.title.y = element_text(size = 20))

ggsave(
    filename = "whelk_N50_plot_v5.png",
    plot = N50.plot,
    device = 'png',
    width = 11, 
    height = 7,
    units = 'in'
    );

write.table(
    x = figure.data,
    file = 'whelk_busco_figure_data.tsv',
    sep = '\t',
    col.names = TRUE,
    row.names = FALSE
    );
