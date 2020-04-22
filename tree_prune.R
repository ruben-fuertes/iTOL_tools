#!/usr/bin/Rscript
args = commandArgs(trailingOnly = T)

# test if there is at least one argument: if not, return an error
if (length(args) == 0) {
        stop("At least one argument must be supplied (input file).n",
             call. = FALSE)
} else if (length(args) == 2) {
        # default output file
        args[3] = "tree/retain.tsv"
}

library(biomformat )

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

write.table(
        seq_prune(biom2, tax)$id,
        file = args[3],
        sep = '\t',
        row.names = F,
        col.names = 'id'
)
