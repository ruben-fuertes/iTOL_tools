#!/usr/bin/Rscript
args = commandArgs(trailingOnly = T)

# test if there is at least one argument: if not, return an error
if (length(args) == 0) {
        stop("At least one argument must be supplied (input file).n",
             call. = FALSE)
} else if (length(args) == 2) {
        # default output file
        args[3] = "tree/binary.txt"
}

library(biomformat)
#Create the dataset for colouring the nodes depending on the
#functional groups


#load the data
biom2 = read_biom(args[1])
tax = read.csv(args[2], sep = '\t', stringsAsFactors = F)[-3]

seq_prune = function(b, tax) {
        #Prune the sequences keeping only one seq that correspond to each
        #taxonomy assignment
        #extract the otu table
        columns = b$columns
        rows = b$rows
        data = b$data
        #loop
        d = data.frame(stringsAsFactors = F)
        for (i in 1:length(rows)) {
                otu_id = rows[[i]][[1]]
                taxa = tax[tax == otu_id, 'taxonomy']
                count = sum(data[[i]])
                if (count < length(columns)) {
                        next
                }
                if (!(taxa %in% rownames(d))) {
                        d[taxa, 'id'] = otu_id
                        d[taxa, 'count'] = count
                        next
                }
                if (count > d[taxa, 'count']) {
                        d[taxa, 'id'] = otu_id
                        d[taxa, 'count'] = count
                }

        }
        return(d)
}
binary = function(d) {
        #Functional groups load
        f = read.csv(
                'input/filtros.csv',
                stringsAsFactors = F
        )
        d2 = data.frame(
                id = rep(NA, length(rownames(d))),
                ACCD = NA,
                Copiotroficos = NA,
                Fosforo = NA,
                MicroNitroactivos = NA,
                Nitroactivos = NA,
                PGPR = NA
        )
        rownames(d2) = rownames(d)
        d2$id = d$id
        for (i in 1:length(d$id)) {
                mo = rownames(d)[i]
                for (filt in levels(as.factor(f$filter))) {
                        a = (grep(
                                paste(
                                        paste0('__', paste(f[f$filter == filt, 'organism'],
                                                           collapse = ';|__'), ';'),
                                        paste0('__', paste(f[f$filter == filt, 'organism'],
                                                           collapse = '$|__'), '$'),
                                        sep = '|'
                                ),
                                ignore.case = T,
                                rownames(d)
                        ))

                        d2[a, filt] = 1
                        d2[-a, filt] = -1
                }
        }
        return(d2)
}

d2 = binary(seq_prune(biom2, tax))

label = colnames(d2)[-1]
color = list('#E4D00A',
             '#CB6D51',
             '#0047AB',
             '#FFB200',
             '#FF9933',
             '#2F4F4F')

a = paste0(
        'DATASET_BINARY

        SEPARATOR COMMA

        DATASET_LABEL,filters

        COLOR,#ff0000

        FIELD_LABELS,',
        label[[1]],
        ',',
        label[[2]],
        ',',
        label[[3]],
        ',',
        label[[4]],
        ',',
        label[[5]],
        ',',
        label[[6]],
        '\n',
        'FIELD_COLORS,',
        color[[1]],
        ',',
        color[[2]],
        ',',
        color[[3]],
        ',',
        color[[4]],
        ',',
        color[[5]],
        ',',
        color[[6]],
        '\n',

        'FIELD_SHAPES,1,1,1,1,1,1

        SHOW_INTERNAL,0

        MARGIN,0

        HEIGHT_FACTOR,1.8

        SYMBOL_SPACING,10

        DATA
        '
)

#Print the final file
write(a, file = args[3])

write.table(
        d2,
        args[3],
        row.names = F,
        sep = ',',
        col.names = F,
        quote = F,
        append = T
)
