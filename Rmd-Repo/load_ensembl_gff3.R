loadEnsemblGff3 =  function(filename, type="gene") {
    gff = read.table(filename,
                     sep='\t', row.names=NULL, header=FALSE,
                     comment.char='#', quote='',
                     check.names=FALSE, stringsAsFactors=FALSE)
    colnames(gff) = c(
        'seqname',
        'source',
        'feature',
        'start',
        'end',
        'score',
        'strand',
        'frame',
        'attribute'
    )
    gff$gene_id = gsub('^.*=gene:(.*?)(;|$).*', '\\1', gff$attribute)
    gff$gene_name = gsub('^.*Name=(.*?)(;|$).*', '\\1', gff$attribute)
    gff[grepl('ID=', gff$gene_name), 'gene_name'] = NA
    gff$biotype = gsub('^.*biotype=(.*?)(;|$).*', '\\1', gff$attribute)
    if (tolower(type) %in% c('gene', 'genes')) {
        genes = gff[gff$feature == 'gene', ]
        rownames(genes) = genes$gene_id
        return(genes)
    } else if (tolower(type) %in% c('transcript', 'transcripts', 'trx')) {
        trx = gff[grepl('RNA|transcript', gff$feature) &
                  !grepl('RNA_gene', gff$feature), ]
        rownames(trx) = gsub('^.*=transcript:(.*?)(;|$).*', '\\1', trx$attribute)
        return(trx)
    }
    return(gff)
}
