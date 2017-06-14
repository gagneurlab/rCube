#' A function to create constitutive features (i.e. exons or introns) from an
#' annotation file. All exonic/intronic bases are considered constitutive, when
#' they belong to all transcript isoform of the corresponding gene.
#' 
#' @param granges A GRanges annotation file with exons/introns. Needs columns 
#' "gene_id" and "transcript_id"
#' @param ncores Number of cores available for parallel computation
#' 
#' @import BiocParallel
#'
#' @return Returns a GRanges object with constitutive exons.
#' @author Carina Demel
#' 
#' @examples
#' # Gencode annotation of MYC gene
#' data(example.exons)
#' constitutive.exons = createConstitutiveFeaturesGRangesFromGRanges(example.exons, 1)
#' @export
createConstitutiveFeaturesGRangesFromGRanges = function(granges, ncores=2){

    #TODO which Param?
    snowparam <- BiocParallel::SnowParam(workers = ncores, type = "SOCK")
    BiocParallel::register(snowparam, default = TRUE)
    BiocParallel::registered()
    
    # build constitutive features (exons/introns):
    # set of (exonic/intronic) bases that belong to each isoform of the gene
    getConstitutiveFeatures <- function(gid){
        ranges <- granges[which(granges$gene_id == gid)]
        tr_ids <- unique(ranges$transcript_id) # transcript ids
        tr_num <- length(tr_ids) # number of transcript isoforms
        
        if(tr_num > 1){
            cov <- coverage(ranges)[[1]]
            w <- which(runValue(cov) == tr_num)
            runLengthsSum <- cumsum(runLength(cov))
            start <- runLengthsSum[w-1] + 1
         
            if(length(start)>0){
                end <- runLengthsSum[w]
                const.feat <- GRanges(seqnames = seqnames(ranges[1]),
                                      ranges = IRanges(start = start, end = end),
                                      strand = strand(ranges[1]),
                                      type = factor("constitutive feature"),
                                      gene_id = rep(gid,length(start)))
            }else{
                const.feat <- GRanges()
            }
        }else{
            const.feat <- ranges
            elementMetadata(const.feat) <- data.frame(type = factor("constitutive feature"),
                                                      gene_id = gid)
        }
        
        return(const.feat)
    }
    
    gene_ids = unique(granges$gene_id)
    const.feature.list <- BiocParallel::bplapply(gene_ids, getConstitutiveFeatures)
    const.feature.list <- GRangesList(const.feature.list)
    constitutive.features <- do.call("c", const.feature.list)
   
    feature.num.formatted <- sprintf("%05d", 1:length(constitutive.features))
    names(constitutive.features) <- paste("CF", feature.num.formatted, sep="")
    
    return(constitutive.features)
}
