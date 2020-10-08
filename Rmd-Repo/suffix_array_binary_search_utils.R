 ## arguments to pwmProbScore:
 ## - pwm is position weight matrix;
 ##   rows should be named by characters from sequence alphabet
 ## - kmer can be either single string or vector of k single characters
pwmProbScore = function(pwm, kmer) {
    if (is.character(kmer) && length(kmer) == 1) {
         ## split kmer up into vector of individual bases:
        kmer = strsplit(kmer, "")[[1]]
    }
     ## use lapply to select the prob value from the i-th column of pwm
     ## for the character kmer[[i]] as the i-th element of ouput vector p:
    p = lapply(1:ncol(pwm), function(i) {pwm[kmer[[i]], i]})
     ## would rather have p as vector (for prob function below) than list;
     ## can use R's unlist function to accomplish this:
    p = unlist(p)
     ## use prod function to multiply p[[1]] * p[[2]] ... * p[[k]],
     ## where k is the number of characters in kmer:
    return(prod(p))
}


bestKmer = function(pwm) {
     ## pick out maximum probability entry from each column:
    kmer = apply(pwm, 2, which.max)
     ## convert kmer from integer representation to character vector:
    kmer = rownames(pwm)[kmer]
     ## paste together with arg collapse="" into single string and return:
    return(paste(kmer, collapse=""))
}


## arguments to nextKmerPrefix:
## - kmer is represented as *integer* vector for efficiency; for example,
##   "TCGAGG" would be represented by c(4, 2, 3, 1, 3, 3)
## - activeIndex is single integer indicating that kmer[1:activeIndex]
##   is the actual prefix under current consideration in enumerateKmers;
##   (the suffix kmer[(activeIndex+1):length(kmer)] is ignored)
nextKmerPrefix = function(kmer, activeIndex) {
    while (activeIndex > 0 && kmer[[activeIndex]] == 4) {
         ## there are only 4 possible nucleotides, so if
         ## kmer[[activeIndex]] is 4, need to reset kmer[[activeIndex]]
         ## to 1 (encoding lexicographically first base A)...
        kmer[[activeIndex]] = 1
         ## and then shift activeIndex back to the left one position:
        activeIndex = activeIndex - 1
    }
    if (activeIndex == 0) {
         ## in this case we've exhausted all possible kmers and are done
        return(NULL)
    }
     ## now advance the base at position activeIndex by 1
     ## (from A to C, C to G, or G to T, depending on current value):
    kmer[[activeIndex]] = kmer[[activeIndex]] + 1
     ## need to return both the updated kmer and activeIndex values:
    return(list(kmer=kmer, activeIndex=activeIndex))
}


enumerateKmers = function(pwm, threshold) {
     ## represent kmer being considered as *integer* vector:
    kmer = rep(1, ncol(pwm))
    nucs = c("A", "C", "G", "T")  ## so nucs[kmer] is actual base sequence
    activeIndex = 1  ## index position in kmer currently under consideration
     ## p[[i]] tracks probability from PWM for kmer base at index i;
     ## initialize for kmer consisting of ncol(pwm) consecutive As:
    p = pwm[1, ]     ## (since A is nucs[[1]])
    out = character(0)   ## good kmer matches will be added to out
    done = FALSE     ## will be reset to TRUE inside while loop when done
    while (!done) {
         ## some passes through loop tell us there cannot be any more kmers
         ## starting with kmer[1:activeIndex] which could be a good match:
         ## in these cases we'll call nextKmerPrefix function to figure out
         ## next kmer prefix to consider, for now assume we will do so:
        callNextKmerPrefix = TRUE  ## but reset to FALSE below if necessary
        if (prod(p[1:activeIndex]) >= threshold) {
             ## then kmer[1:activeIndex] could be prefix of good match(es):
            if (activeIndex == ncol(pwm)) {  ## we've confirmed right-most
                 ## position, thus kmer is a good match to pwm, so paste
                 ## nucs[kmer] together with collapse="" and add to out:
                out[[length(out)+1]] = paste(nucs[kmer], collapse="")
            } else {  ## haven't checked full kmer yet; check next position:
                activeIndex = activeIndex + 1
                 ## and make sure p[[activeIndex]] correctly set according
                 ## to whatever base is encoded by kmer[[activeIndex]]:
                p[[activeIndex]] = pwm[kmer[[activeIndex]], activeIndex]
                 ## still working on extending this kmer, so:
                callNextKmerPrefix = FALSE
            }
        }
        if (callNextKmerPrefix) {
            nextInfo = nextKmerPrefix(kmer, activeIndex)
            if (length(nextInfo) == 0) {  ## signals we're done:
                done = TRUE
            } else {
                kmer = nextInfo$kmer
                activeIndex = nextInfo$activeIndex
                 ## make sure p[[activeIndex]] correctly set according
                 ## to whatever base is encoded by kmer[[activeIndex]]:
                p[[activeIndex]] = pwm[kmer[[activeIndex]], activeIndex]
            }
        }
    }
    return(out)
}


 ## arguments to binarySearchSuffixArray:
 ## - kmer: the kmer we want to find the locations of in the:
 ## - sequence: the single string to search for kmer in
 ## - suffixArray: integer vector suffix array for sequence
 ## - comparator: should be either `<` or `<=`; function will
 ##               find lower end of suffix array block if `<`
 ##               or upper end of suffix array block if `<=`
binarySearchSuffixArray = function(kmer, sequence, suffixArray,
                                   comparator=`<`) {
    k = nchar(kmer)
    lower = 0
    upper = length(suffixArray) + 1
    middle = round(mean(c(lower, upper)))
    while (lower < (upper-1)) {
        kmMid = substr(sequence,
                       suffixArray[[middle]], suffixArray[[middle]]+k-1)
        if (comparator(kmMid, kmer)) {
            lower = middle
        } else {
            upper = middle
        }
        middle = round(mean(c(lower, upper)))
    }
     ## need to figure out whether to return lower or upper:
     ## lower <= lower but not lower < lower, so use:
    if (comparator(lower, lower)) {  ## to test what comparator is
        return(lower)  ## in this case comparator is <=, so return lower
    } else {
        return(upper)  ## in this case comparator is <. so return upper
    }
}


findKmer = function(kmer, sequence, suffixArray) {
    lower = binarySearchSuffixArray(kmer, sequence, suffixArray, `<`)
    upper = binarySearchSuffixArray(kmer, sequence, suffixArray, `<=`)
    if (upper < lower) {    ## happens when no match found
        return(integer(0))  ## return empty numeric vector
    }
     ## now just need to use suffix array to get "spatial" positions of
     ## suffixes beginning with kmer---i.e., locations of kmer in sequence:
    return(suffixArray[lower:upper])
}



## -----------------------------------------------------------------------------
sdPwm = rbind(
    c(0.7, 0.1, 0.1, 0.7, 0.1, 0.1),
    c(0.1, 0.1, 0.1, 0.1, 0.1, 0.1),
    c(0.1, 0.7, 0.7, 0.1, 0.7, 0.7),
    c(0.1, 0.1, 0.1, 0.1, 0.1, 0.1)
)
rownames(sdPwm) = c("A", "C", "G", "T")
