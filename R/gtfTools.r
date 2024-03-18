require(data.table)
require(GenomicRanges)
require(IRanges)
require(Biostrings)

infoExtract <- function(gtf, infoColName, key, type) {

    n <- length(key)
    if (n != length(type)) stop("key & type must be of same length")

    # parse info into dataframe keeping track of feature index
    info <- strsplit(gtf[[infoColName]], ";")
    infodt <- data.table(
        info = unlist(info),
        idx = unlist(mapply(function(v, i) rep(i, length(v)), info, 1:length(info)))
    )
    infodt$info <- gsub("\"", "", infodt$info)
    infodt$info <- trimws(infodt$info)

    # format keys for matching
    fkey <- ifelse(type == "tag", paste("tag", key), paste0("^", key, " "))

    # match
    for (i in 1:n) {

        if (type[i] == "tag") {

            r <- getTag(infodt, fkey[i], nrow(gtf))

        } else {

            r <- getMulti(infodt, fkey[i], nrow(gtf))

        }
        gtf[[ key[i] ]] <- r

    }

    return(gtf)

}

getMulti <- function(infodt, multi, nfeatures) {

    # get values for the queried multilevel key for each feature in gtf
    r <- rep(NA_character_, nfeatures)
    ismulti <- grepl(multi, infodt$info)
    r[infodt$idx[ismulti]] <- gsub(multi, "", infodt$info[ismulti])

    return(r)

}

getTag <- function(infodt, tag, nfeatures) {

    # get values for the queried binary key for each feature in gtf
    r <- rep(0L, nfeatures)
    r[infodt$idx[infodt$info == tag]] <- 1L

    return(r)

}

# gets features that appear in the mRNA between the start and stop codon (inclusive) in transcription order
cdsfy <- function(gtf) {

    gtfFilter <- gtf[Feature == "CDS" | Feature == "start_codon" | Feature == "stop_codon"]
    gtfStrand <- split(gtfFilter, gtfFilter$Strand)
    gtfStrand <- lapply(names(gtfStrand), function(s) {
        .dt <- split(gtfStrand[[s]], gtfStrand[[s]]$Chromosome)
        if (s == "+") .dt <- lapply(.dt, function(x) x[order(Start_Position)])
        if (s == "-") .dt <- lapply(.dt, function(x) x[order(Start_Position, decreasing = TRUE)])
        rbindlist(.dt)
    })

    gtfMrna <- rbindlist(gtfStrand)
    return(gtfMrna)

}

# gets features that get transcribed to RNA in transcription order
rnaify <- function(gtf, mrnaOnly = FALSE) {
    
    if (mrnaOnly) {
    
        gtfFilter <- gtf[Feature == "CDS" | Feature == "UTR"]

    } else {

        gtfFilter <- gtf[Feature == "exon"]

    }
    

    gtfStrand <- split(gtfFilter, gtfFilter$Strand)
    gtfStrand <- lapply(names(gtfStrand), function(s) {
        .dt <- split(gtfStrand[[s]], gtfStrand[[s]]$Chromosome)
        if (s == "+") .dt <- lapply(.dt, function(x) x[order(Start_Position)])
        if (s == "-") .dt <- lapply(.dt, function(x) x[order(Start_Position, decreasing = TRUE)])
        rbindlist(.dt)
    })

    gtfMrna <- rbindlist(gtfStrand)
    return(gtfMrna)

}

# gets the 5'UTR of a mRNA gtf
get5utr <- function(mrnaGtf) {

    # split gtf by strand
    strandGtf <- split(mrnaGtf, mrnaGtf$Strand)

    utr5Gtf <- lapply(
        names(strandGtf),
        function(s) {
            cdsGtf <- strandGtf[[s]][Feature == "CDS"]
            utrGtf <- strandGtf[[s]][Feature == "UTR"]
            if (s == "+") {
                firstcds <- tapply(cdsGtf$Start_Position, cdsGtf$transcript_id, min)
                return(utrGtf[utrGtf$End_Position < firstcds[utrGtf$transcript_id]])
            }
            firstcds <- tapply(cdsGtf$End_Position, cdsGtf$transcript_id, max)
            return(utrGtf[utrGtf$Start_Position > firstcds[utrGtf$transcript_id]])
        }
    )

    return(rbindlist(utr5Gtf))

}

# gets the 5'UTR of a mRNA gtf
get3utr <- function(mrnaGtf) {

    # split gtf by strand
    strandGtf <- split(mrnaGtf, mrnaGtf$Strand)

    utr3Gtf <- lapply(
        names(strandGtf),
        function(s) {
            cdsGtf <- strandGtf[[s]][Feature == "CDS"]
            utrGtf <- strandGtf[[s]][Feature == "UTR"]
            if (s == "+") {
                lastcds <- tapply(cdsGtf$End_Position, cdsGtf$transcript_id, max)
                return(utrGtf[utrGtf$Start_Position > lastcds[utrGtf$transcript_id]])
            }
            lastcds <- tapply(cdsGtf$Start_Position, cdsGtf$transcript_id, min)
            return(utrGtf[utrGtf$End_Position < lastcds[utrGtf$transcript_id]])
        }
    )

    return(rbindlist(utr3Gtf))

}

# gets the introns of a mRNA gtf
getIntron <- function(mrnaGtf) {

    mrnaGtf$transcript_id <- factor(mrnaGtf$transcript_id, levels = unique(mrnaGtf$transcript_id))
    exonSeq <- paste(mrnaGtf$Chromosome, mrnaGtf$transcript_id, sep = "_")
    exon <- GRanges(exonSeq, IRanges(mrnaGtf$Start_Position, mrnaGtf$End_Position))
    trStart <- tapply(mrnaGtf$Start_Position, mrnaGtf$transcript_id, min)
    trEnd <- tapply(mrnaGtf$End_Position, mrnaGtf$transcript_id, max)
    tr <- GRanges(unique(exonSeq), IRanges(trStart, trEnd))
    intron <- setdiff(tr, exon)
    intron <- data.table(
        Chromosome = sapply(strsplit(as.character(seqnames(intron)), "_"), "[", 1),
        Start_Position = start(intron),
        End_Position = end(intron)
    )
    return(intron)

}

# collapse all exons and introns of transcripts to 1 annotation
collapse2tx <- function(gtf, trcol) {

    trStart <- tapply(gtf$Start_Position, gtf[[trcol]], min)
    trEnd <- tapply(gtf$End_Position, gtf[[trcol]], max)
    trIdx <- match(names(trStart), gtf[[trcol]])
    r <- data.table(
        Chromosome = gtf$Chromosome[trIdx],
        Start_Position = trStart,
        End_Position = trEnd
    )
    r[[trcol]] <- gtf[[trcol]][trIdx]
    carrycols <- names(gtf)[!(names(gtf) %in% c("Chromosome", "Start_Position", "End_Position", "Feature"))]
    r <- cbind(r, gtf[trIdx, ..carrycols])
    r$Feature <- "transcript"

    return(r)

}

gtf2dna <- function(gtf, genome) {

    gtfRanges <- GenomicRanges::GRanges(
        gtf$Chromosome,
        IRanges::IRanges(gtf$Start_Position, gtf$End_Position)
    )
    dna <- genome[gtfRanges]
    isminus <- gtf$Strand == "-"
    dna[isminus] <- Biostrings::reverseComplement(dna[isminus])

    return(dna)

}

utr5width <- function(gtfMrna) {

    firstCDSidx <- sapply(gtfMrna, function(.dt) min(which(.dt$Feature == "CDS")))
    utr5w <- rep(0L, length(firstCDSidx))
    hasutr <- firstCDSidx > 1
    utr5w[hasutr] <- mapply(
        function(.dt, i) {
            ii <- 1:(i - 1)
            sum(.dt$End_Position[ii] - .dt$Start_Position[ii] + 1)
        },
        gtfMrna[hasutr],
        firstCDSidx[hasutr]
    )

    return(utr5w)

}

cdsWidth <- function(gtfMrna) {

    cdsidx <- lapply(gtfMrna, function(.dt) which(.dt$Feature == "CDS"))
    hascds <- sapply(cdsidx, length) > 0
    cdsw <- rep(0L, length(cdsidx))
    cdsw[hascds] <- mapply(
        function(.dt, ii) sum(.dt$End_Position[ii] - .dt$Start_Position[ii] + 1),
        gtfMrna[hascds],
        cdsidx[hascds]
    )

    return(cdsw)

}

filterOvTx <- function(exonGtf) {

    # find overlaps between exons 
    exonRanges <- GRanges(exonGtf$Chromosome,  IRanges(exonGtf$Start_Position, exonGtf$End_Position))
    exonOv <- findOverlaps(exonRanges, exonRanges)

    # count overlaps with own gene or other genes for each tx
    ovdt <- data.table(
        tx = exonGtf$transcript_id[queryHits(exonOv)],
        ov = ifelse(
            exonGtf$gene_id[queryHits(exonOv)] == exonGtf$gene_id[subjectHits(exonOv)],
            "own",
            "other"
        )
    )
    count <- table(
        tx = factor(ovdt$tx, unique(exonGtf$transcript_id)),
        ov = factor(ovdt$ov, c("own", "other"))
    )

    # filter tx that overlap with tx from other genes
    filterGtf <- exonGtf[transcript_id %in% rownames(count)[count[, "other"] == 0L]]

    return(filterGtf)

}

snpAnnotate <- function(snpdt, exonGtf) {

    # annotation vector
    r <- rep(NA_character_, nrow(snpdt))

    # compute gtf-snp overlap
    snpRanges <- GRanges(snpdt$Chromosome, IRanges(snpdt$Start_Position, snpdt$Start_Position))
    gtfRanges <- GRanges(exonGtf$Chromosome, IRanges(exonGtf$Start_Position, exonGtf$End_Position))
    gtfSnpOv <- findOverlaps(gtfRanges, snpRanges)
    if (length(gtfSnpOv) == 0) return(r)
    
    # get number of snps per tx
    txdt <- as.data.table(table(tx = exonGtf$transcript_id[queryHits(gtfSnpOv)]))

    # add info that will help decide which transcript to use per gene
    exonGtf[, "width" := width(gtfRanges)]
    exonGtf[
        exonGtf[, list("txwidth" = sum(width)), by = "transcript_id"],
        "txwidth" := i.txwidth,
        on = "transcript_id"
    ]
    txdt[exonGtf, ':=' ("gene" = i.gene_id, "level" = i.level, "width" = i.txwidth), on = c("tx" = "transcript_id")]
    
    # get best snp number, evidence level and length across tx of each gene
    bestdt <- txdt[, list("maxN" = max(N), "maxLevel" = min(level), "maxWidth" = max(width)), by = "gene"]
    txdt[bestdt, ':=' ("maxN" = i.maxN, "maxLevel" = i.maxLevel, "maxWidth" = i.maxWidth), on = "gene"]

    # choose 1 tx per gene
    chosedt <- txdt[N == maxN & level == maxLevel & width == maxWidth, list("tx" = head(tx, 1)), by = "gene"]
    finalGtf <- exonGtf[transcript_id %in% chosedt$tx]

    # assign a tx to each snp
    finalRanges <-  GRanges(finalGtf$Chromosome, IRanges(finalGtf$Start_Position, finalGtf$End_Position))
    snpGtfOv <- findOverlaps(snpRanges, finalRanges)
    r[queryHits(snpGtfOv)] <- finalGtf$transcript_id[subjectHits(snpGtfOv)]

    return(r)

}

# given a gtf with the neccessary columns, it gets the variant classification of every possible mutation in the gtf's transcripts
pvarAnnotate <- function(gtf, genome) {

    # get variant classification for protein coding sites
    cdsgtf <- gtf[Feature %in% c("CDS",  "stop_codon")]
    cdspmutdt <- pvarCDSannotate(cdsgtf, genome)
    cdspmutdt[, c("relpos", "frame", "codon", "wtAA", "mutCodon", "mutAA") := NULL]
    
    # get variant classification for 5'UTR sites
    mrnagtf <- rnaify(gtf, TRUE)
    utr5gtf <- get5utr(mrnagtf)
    utr5pmutdt <- pvarNCannotate(utr5gtf, genome)
    utr5pmutdt$type <- rep("utr5", nrow(utr5pmutdt))

    # get variant classification for 3'UTR sites
    utr3gtf <- get3utr(mrnagtf)
    utr3pmutdt <- pvarNCannotate(utr3gtf, genome)
    utr3pmutdt$type <- rep("utr3", nrow(utr3pmutdt))

    # get variant classification for non-coding sites
    ncgtf <- gtf[Feature == "exon" & gene_type != "protein_coding"]
    ncpmutdt <- pvarNCannotate(ncgtf, genome)
    ncpmutdt$type <- ncgtf$gene_type[match(ncpmutdt$transcript_id, ncgtf$transcript_id)]

    return(rbind(cdspmutdt, utr5pmutdt, utr3pmutdt, ncpmutdt))

}

pvarNCannotate <- function(ncgtf, genome) {

    # get coding sites, their transcript of origin and strand
    gtfRanges <- GRanges(ncgtf$Chromosome, IRanges(ncgtf$Start_Position, ncgtf$End_Position))
    siteRanges <- unlist(slidingWindows(gtfRanges, 1))
    siteRanges$transcript_id <- rep(ncgtf$transcript_id, width(gtfRanges))
    siteRanges$codingStrand <- rep(ncgtf$Strand, width(gtfRanges))

    # get nucleotides in the + strand
    siteRanges$ref <- genome[siteRanges]

    # get nucleotides in the coding strand
    siteRanges$codingRef <- siteRanges$ref
    reverse <- siteRanges$codingStrand == "-"
    siteRanges$codingRef[reverse] <- reverseComplement(siteRanges$codingRef[reverse])

    # determine pyrimidine strand and get pyrimidine
    siteRanges$pyrimidineStrand <- ifelse(siteRanges$ref %in% c("C", "T"), "+", "-")
    siteRanges$pyrimidine <- siteRanges$ref
    reverse <- siteRanges$pyrimidineStrand == "-"
    siteRanges$pyrimidine[reverse] <- reverseComplement(siteRanges$pyrimidine[reverse])

    # get all possible mutations in the examined sited oriented by the pyrimidine
    pmuts <- list(C = c("A", "G", "T"), T = c("A", "C", "G"))[as.character(siteRanges$pyrimidine)]
    pmutRanges <- siteRanges[rep(1:length(siteRanges), each = 3)]
    pmutRanges$pyrimidineMut <- DNAStringSet(unlist(pmuts))

    # get the mutation in the coding strand according to pyrimidine strand
    reverse <- pmutRanges$pyrimidineStrand != pmutRanges$codingStrand
    pmutRanges$codingMut <- pmutRanges$pyrimidineMut
    pmutRanges$codingMut[reverse] <- reverseComplement(pmutRanges$codingMut[reverse])

    # convert to dta.table
    names(pmutRanges) <- 1:length(pmutRanges)
    pmutdt <- data.table(as.data.frame(pmutRanges))
    pmutdt[, c("width", "strand", "end") := NULL]
    names(pmutdt)[2] <- "position"

    return(pmutdt)

}

# assumes all transcripts have been checked to be divisible by 3!
pvarCDSannotate <- function(cdsgtf, genome) {

    # ensure order of this chromosome
    cdsgtf <- cdsgtf[order(Start_Position)]

    # get coding sites, their transcript of origin and strand
    gtfRanges <- GRanges(cdsgtf$Chromosome, IRanges(cdsgtf$Start_Position, cdsgtf$End_Position))
    siteRanges <- unlist(slidingWindows(gtfRanges, 1))
    siteRanges$transcript_id <- rep(cdsgtf$transcript_id, width(gtfRanges))
    siteRanges$codingStrand <- rep(cdsgtf$Strand, width(gtfRanges))

    # split sites by transcript, if they are in negative strand reverse order so that everything is in transcription order
    siteRanges <- split(siteRanges, siteRanges$transcript_id)
    reverse <- sapply(siteRanges, function(gr) gr$codingStrand[1]) == "-"
    siteRanges[reverse] <- endoapply(siteRanges[reverse], rev)

    # get relative position of each site within its transcript CDS
    relpos <- lapply(elementNROWS(siteRanges), function(l) 1:l)
    siteRanges <- unlist(siteRanges)
    siteRanges$relpos <- unlist(relpos)

    # get nucleotides in the + strand
    siteRanges$ref <- genome[siteRanges]

    # get nucleotides in the coding strand
    siteRanges$codingRef <- siteRanges$ref
    reverse <- siteRanges$codingStrand == "-"
    siteRanges$codingRef[reverse] <- reverseComplement(siteRanges$codingRef[reverse])

    # determine pyrimidine strand and get pyrimidine
    siteRanges$pyrimidineStrand <- ifelse(siteRanges$ref %in% c("C", "T"), "+", "-")
    siteRanges$pyrimidine <- siteRanges$ref
    reverse <- siteRanges$pyrimidineStrand == "-"
    siteRanges$pyrimidine[reverse] <- reverseComplement(siteRanges$pyrimidine[reverse])

    # get the codon and aminoacid that each site influences
    siteRanges$frame <- ((siteRanges$relpos + 2L) %% 3L) + 1L # gives 1 if in frame, 2 for one base out of frame and 3 for two bases out of frame
    siteRanges$codon <- rep(
        xscat(
            siteRanges$codingRef[siteRanges$frame == 1],
            siteRanges$codingRef[siteRanges$frame == 2],
            siteRanges$codingRef[siteRanges$frame == 3]
        ),
        each = 3
    )
    siteRanges$wtAA <- translate(siteRanges$codon, no.init.codon = TRUE)
    
    # get all possible mutations in the examined sites oriented by the pyrimidine
    pmuts <- list(C = c("A", "G", "T"), T = c("A", "C", "G"))[as.character(siteRanges$pyrimidine)]
    pmutRanges <- siteRanges[rep(1:length(siteRanges), each = 3)]
    pmutRanges$pyrimidineMut <- DNAStringSet(unlist(pmuts))

    # get the mutation in the coding strand according to pyrimidine strand
    reverse <- pmutRanges$pyrimidineStrand != pmutRanges$codingStrand
    pmutRanges$codingMut <- pmutRanges$pyrimidineMut
    pmutRanges$codingMut[reverse] <- reverseComplement(pmutRanges$codingMut[reverse])

    # get the mutated codon and mutated aminoacid according to mutation in coding strand
    pmutRanges$mutCodon <- pmutRanges$codon
    subseq(pmutRanges$mutCodon, pmutRanges$frame, pmutRanges$frame) <- pmutRanges$codingMut
    pmutRanges$mutAA <- translate(pmutRanges$mutCodon, no.init.codon = TRUE)

    # determine the type of mutation
    pmutRanges$type <- rep("missense", length(siteRanges))
    pmutRanges$type[pmutRanges$wtAA == pmutRanges$mutAA] <- "syn"
    pmutRanges$type[pmutRanges$wtAA != "*" & pmutRanges$mutAA == "*"] <- "nonsense"
    pmutRanges$type[pmutRanges$wtAA == "*" & pmutRanges$mutAA != "*"] <- "nonstop"

    # convert to dta.table
    names(pmutRanges) <- 1:length(pmutRanges)
    pmutdt <- data.table(as.data.frame(pmutRanges))
    pmutdt[, c("width", "strand", "end") := NULL]
    names(pmutdt)[2] <- "position"

    return(pmutdt)

}
