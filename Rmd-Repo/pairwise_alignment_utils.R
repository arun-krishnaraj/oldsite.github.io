bestFinalStep = function(v, w, cost,
                         subCost=1, insertCost=1, deleteCost=1,
                         matchCost=0, maxCost=Inf) {
    finalPositionMismatched = (v[[length(v)]] != w[[length(w)]])
    finalPositionMatched = !finalPositionMismatched
    i = length(v) + 1
    j = length(w) + 1
    costOptions = c(
        M = cost[i-1, j-1] +
            finalPositionMismatched * subCost +
            finalPositionMatched * matchCost,
        D = cost[i-1, j] + deleteCost,
        I = cost[i, j-1] + insertCost
    )
    bestEdit = names(which.min(costOptions))
    return(list(
        edit = bestEdit,
        cost = min(costOptions[[bestEdit]], maxCost)
    ))
}


buildCostMatrix = function(v, w,
                           subCost=1, indelCost=1,
                           matchCost=0, trimCost=1, maxCost=Inf) {
    cost = matrix(0, nrow=length(v)+1, ncol=length(w)+1)
    priorEdit = matrix("", nrow=length(v)+1, ncol=length(w)+1)    
    priorEdit[2:nrow(priorEdit), 1] = "D"
    cost[ , 1] = pmin(0:length(v) * trimCost, maxCost)
    priorEdit[1, 2:ncol(priorEdit)] = "I"
    cost[1, ] = pmin(0:length(w) * trimCost, maxCost)
    for (i in 2:nrow(cost)) {
        for (j in 2:ncol(cost)) {
            insertCost = deleteCost = indelCost
            if (i == nrow(cost)) {insertCost = trimCost}
            if (j == ncol(cost)) {deleteCost = trimCost}
            priorStep = bestFinalStep(v[1:(i-1)], w[1:(j-1)], cost,
                                      subCost, insertCost, deleteCost,
                                      matchCost, maxCost)
            cost[i, j] = priorStep[["cost"]]
            priorEdit[i, j] = priorStep[["edit"]]
        }
    }
    return(list(
        cost = cost,
        priorEdit = priorEdit
    ))
}


traceSteps = function(cost, priorEdit, maxCost=Inf) {
    i = nrow(cost)
    j = ncol(cost)
    prefix = suffix = character(0)
    if (maxCost < Inf) {
        ij = max(which(cost == min(cost)))
        i = 1 + ((ij-1) %% nrow(cost))
        j = ceiling(ij / nrow(cost))
        if (i < nrow(cost)) {
            suffix = c(suffix, rep("D", nrow(cost)-i))
        }
        if (j < ncol(cost)) {
            suffix = c(suffix, rep("I", ncol(cost)-j))
        }
    }
    bestPathBack = character(0)
    localDone = FALSE
    while (i > 1 && j > 1 && !localDone) {
        if (priorEdit[i, j] == "M") {
            i = i-1
            j = j-1
            bestPathBack = c(bestPathBack, "M")
        } else if (priorEdit[i, j] == "D") {
            i = i-1
            bestPathBack = c(bestPathBack, "D")
        } else if (priorEdit[i, j] == "I") {
            j = j-1
            bestPathBack = c(bestPathBack, "I")
        }
        localDone = (cost[i, j] >= maxCost)
    }
    if (i > 1) {
        prefix = c(prefix, rep("D", i-1))
    }
    if (j > 1) {
        prefix = c(prefix, rep("I", j-1))
    }
    return(c(prefix, rev(bestPathBack), suffix))
}


cigar = function(cost, priorEdit, maxCost=Inf) {
    forwardSteps = traceSteps(cost, priorEdit, maxCost=maxCost)
    cigar = ""
    activeStep = forwardSteps[[1]]
    nConsecutiveSame = 1
    for (i in 2:length(forwardSteps)) {
        if (forwardSteps[[i]] != activeStep) {
            cigar = paste0(cigar, nConsecutiveSame, activeStep)
            activeStep = forwardSteps[[i]]
            nConsecutiveSame = 1
        } else {
            nConsecutiveSame = nConsecutiveSame + 1
        }
    }
    cigar = paste0(cigar, nConsecutiveSame, activeStep)
    return(cigar)
}


editDistanceAlign = function(v, w,
                             subCost=1, indelCost=1,
                             matchCost=0, trimCost=1, maxCost=Inf) {
    vVec = strsplit(v, split="")[[1]]
    wVec = strsplit(w, split="")[[1]]
    costAndPriorEdits = buildCostMatrix(
        vVec, wVec,
        subCost, indelCost, matchCost, trimCost, maxCost
    )
    i = nrow(costAndPriorEdits$cost)
    j = ncol(costAndPriorEdits$cost)
    stepsChar = do.call(traceSteps,
                        c(costAndPriorEdits, list(maxCost=maxCost)))
    i = j = 1
    steps = list(c(i, j))
    for (stepChar in stepsChar) {
        if (stepChar == "M") {
            i = i+1; j = j+1
        } else if (stepChar == "D") {
            i = i+1
        } else if (stepChar == "I") {
            j = j+1
        }
        steps = c(steps, list(c(i, j)))
    }
    steps = do.call(rbind, steps)    
    return(list(
        v = v,
        w = w,
        cost = costAndPriorEdits$cost,
        priorEdit = costAndPriorEdits$priorEdit,
        cigar = do.call(cigar, c(costAndPriorEdits, list(maxCost=maxCost))),
        steps = steps
    ))
}


plotAlignmentGrid = function(alignment,
                             costs=TRUE, numbers=FALSE, local=FALSE,
                             alpha=1, vjust=2.75, hjust=-4.25) {
    require(ggplot2)
    require(tidyr)
    v = alignment$v
    w = alignment$w
    priorV = matrix(1, nrow=nrow(alignment$cost), ncol=ncol(alignment$cost))
    priorW = matrix(1, nrow=nrow(alignment$cost), ncol=ncol(alignment$cost))
    for (i in 1:nrow(priorV)) {
        for (j in 1:ncol(priorW)) {
            if (i != 1 || j != 1) {
                if (alignment$priorEdit[i, j] == "M") {
                    priorV[i, j] = i-1
                    priorW[i, j] = j-1
                } else if (alignment$priorEdit[i, j] == "D") {
                    priorV[i, j] = i-1
                    priorW[i, j] = j
                } else if (alignment$priorEdit[i, j] == "I") {
                    priorV[i, j] = i
                    priorW[i, j] = j-1
                }
            }
        }
    }
    steps = alignment$steps
    if (local) {
        while (steps[1, 1] == steps[2, 1] ||
               steps[1, 2] == steps[2, 2]) {
            steps = steps[2:nrow(steps), ]
        }
        while (steps[nrow(steps)-1, 1] == steps[nrow(steps), 1] ||
               steps[nrow(steps)-1, 2] == steps[nrow(steps), 2]) {
            steps = steps[1:(nrow(steps)-1), ]
        }
    }
    ggd = rbind(
        data.frame(i=1:(nchar(v)+1), priorV,
                   type='i*', check.names=FALSE) %>%
            pivot_longer(c(-i, -type), names_to='j', values_to='parent') %>%
            as.data.frame,
        data.frame(i=1:(nchar(v)+1), priorW,
                   type='j*', check.names=FALSE) %>%
            pivot_longer(c(-i, -type), names_to='j', values_to='parent') %>%
            as.data.frame
    ) %>% pivot_wider(names_from=type, values_from=parent) %>% as.data.frame
    for (column in colnames(ggd)) {ggd[[column]] = as.integer(ggd[[column]])}
    ggd$best = FALSE
    ggd$match = ''
    for (s in 1:(nrow(steps)-1)) {
        ggd[ggd$'i*' == steps[s, 1] & ggd$'j*' == steps[s, 2] &
            ggd$i == steps[s+1, 1] & ggd$j == steps[s+1, 2],
            'best'] = TRUE
        if (all(steps[s+1, ] - steps[s, ] > 0) &&
            (substr(v, steps[s+1, 1]-1, steps[s+1, 1]-1) ==
             substr(w, steps[s+1, 2]-1, steps[s+1, 2]-1))) {
            ggd[ggd$'i*' == steps[s, 1] & ggd$'j*' == steps[s, 2] &
                ggd$i == steps[s+1, 1] & ggd$j == steps[s+1, 2],
                'match'] = substr(v, steps[s+1, 1]-1, steps[s+1, 1]-1)
        }
    }
    for (column in c('i', 'j', 'i*', 'j*')) {ggd[[column]] = ggd[[column]] - 0.5}
    ggd = ggd[ggd$i != 0.5 | ggd$j != 0.5, ]
    ggo = ggplot(ggd, aes(x=`i*`, y=`j*`, color=best, label=match))    
    if (costs) {
	    costMelt = data.frame(`i*` = 1:nrow(alignment$cost),
                                  alignment$cost,
                                  check.names=FALSE) %>%
	               pivot_longer(-`i*`, names_to='j*', values_to='cost') %>%
	               as.data.frame
	    costMelt$'i*' = costMelt$'i*' - 0.5
	    costMelt$'j*' = as.numeric(as.character(costMelt$'j*')) - 0.5
        ggo = ggo + geom_raster(aes(x=`i*`, y=`j*`, fill=cost),
                                data=costMelt, inherit.aes=FALSE)
        ggo = ggo + scale_fill_gradientn(colors=c('white',
                                                  'lightgoldenrod',
                                                  'darkgoldenrod'))
    }
    ggo = ggo + geom_text(
        data = data.frame('i*'=1:nchar(v), match=strsplit(v, '')[[1]],
                          'j*'=0, best=FALSE, check.names=FALSE),
        color = 'black'
    ) + geom_text(
        data = data.frame('j*'=1:nchar(w), match=strsplit(w, '')[[1]],
                          'i*'=0, best=FALSE, check.names=FALSE),
        color = 'black'
    )            
    ggo = ggo + geom_text(
        data = data.frame('j*'=1:nchar(w), match=strsplit(w, '')[[1]],
                          'i*'=0, best=FALSE, check.names=FALSE),
        color = 'black'
    )
    ggo = ggo + geom_segment(aes(xend=i, yend=j), alpha=alpha,
                             arrow=arrow(length=unit(0.075, units='inches')))
    if (length(vjust) * length(hjust) > 0) {
        ggo = ggo + geom_text(vjust=vjust, hjust=hjust, alpha=alpha)
    }
    if (costs && numbers) {
        ggo = ggo + geom_text(aes(x=`i*`, y=`j*`, label=cost),
                              data=costMelt, inherit.aes=FALSE)
    }
    ggo = ggo + xlab(expression(italic('i'))) + ylab(expression(italic('j')))
    ggo = ggo + scale_y_reverse()
    ggo = ggo + theme_bw()
    ggo = ggo + theme(axis.ticks.x = element_blank(),
                      axis.ticks.y = element_blank())
    if (costs) {
        ggo = ggo + scale_color_manual(values=c('darkgray', 'firebrick'), guide=FALSE)
    } else {
        ggo = ggo + scale_color_manual(values=c('darkgray', 'firebrick'), guide=FALSE)
    }
    invisible(ggo)   
}
