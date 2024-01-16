#' @title f
#' @description Computes the F_beta score given an intenger number of True Positives (TP), True Negatives (TN). It is optimized for speed and n is thus not the total number of events
#' @param TP Number of true positive events
#' @param TN Number of true negative events
#' @param n beta^2*(TP+FN)+TN+FP
#' @param beta2 squared-beta to weight precision (low beta) or recall (high beta) more
f<-function(TP,TN,n,beta2=1){
    (1+beta2)*TP/(n+TP-TN) ##n equals beta2*(TP+FN)+(TN+FP)=beta2*nT+nF
}

#' @title fill_FNTN_matrix
#' @description fill_FNTN_matrix Used for assessing whether an expansion move is possible
#' @param xp_FN Expression matrix of False Negative events
#' @param xp_TN Expression matrix of True Negative events
#' @param B_FN Boolean matrix of FN events
#' @param B_TN Boolean matrix of TN events
#' @param par Current hyper-rectangle parametrization
fill_FNTN_matrix<-function(xp_FN,xp_TN,B_FN,B_TN,par){
    FN=nrow(xp_FN)
    TN=nrow(xp_TN)
    N=FN+TN
    xp=rbind(xp_FN,xp_TN)
    B=rbind(B_FN,B_TN)

    res=matrix(nrow=N,ncol=FN,data=FALSE)

    if(ncol(B)<=.Machine$double.digits){
        summer=matrix(nrow=ncol(B),ncol=1,data=2^(1:ncol(B)-1))
        hashes=(!B)%*%summer
        ##Split events according to the state for which channels they are FALSE
        B_hash=split(1:N,hashes)
        B_FN_hash=split(1:FN,hashes[1:FN])

        hash_table=matrix(nrow=length(B_hash),ncol=ncol(B),dimnames=list(names(B_hash),colnames(B)),data=NA)
        for(hash in names(B_hash)){
            hash_table[hash,]=B[B_hash[[hash]][1],]
        }
    } else {

        hashes=apply(!B,1,function(x)paste(x,collapse="/")) ##For every active channel, identify which are in FALSE state.
        unique_hashes=unique(hashes)
        hash_table=matrix(nrow=length(unique_hashes),ncol=ncol(B),dimnames=list(unique_hashes,colnames(B)),data=TRUE)
        for(hash in unique_hashes){
            hash_table[hash,as.logical(strsplit(hash,split="/")[[1]])]=FALSE
        }
        ##hash_table=unique(B)
        ##rownames(hash_table)=apply(hash_table,1,function(x)paste(which(!x),collapse="/"))

        ##Split events according to the state for which channels they are FALSE
        B_hash=split(1:N,hashes)
        B_FN_hash=split(1:FN,hashes[1:FN])
    }

    for(hash in names(B_FN_hash)){
        true_for_hash=hash_table[hash,] ##Check for a given hash which channels are TRUE
        compatible_candidates=rownames(hash_table)[rowSums(hash_table[,true_for_hash,drop=FALSE])==sum(true_for_hash)] ##Only other hashes for which only a subset of FALSE channels are FALSE are compatible (for the others there are guaranteed to not be implicated

        for(other_hash in compatible_candidates){
            others=B_hash[[other_hash]] ##Compatible events
            xp.tmp=xp[others,,drop=FALSE] ##Subsetting the matrix of events to only keep compatible events
            res.tmp=matrix(nrow=length(others),ncol=length(B_FN_hash[[hash]]))
            i=0
            for(event in B_FN_hash[[hash]]){
                i=i+1
                affected_parameters=xp[event,]<par
                affected_parameters.values=xp[event,affected_parameters]
                sum_affected=sum(affected_parameters) ##Number of affected parameters
                res.tmp[,i]=(rowSums(xp.tmp[,affected_parameters,drop=FALSE]>=matrix(nrow=length(others),ncol=sum_affected,data=affected_parameters.values,byrow=TRUE))==rep(sum_affected,nrow(xp.tmp))) ##If all affected parameters are higher than the values for the event screened, the other event is implicated by the event screened
            }
            res[others,B_FN_hash[[hash]]]=res.tmp
        }
    }
    res
}

#' @title plot_gating_strategy
#' @description Plot a hypergate return
#' @param gate A hypergate object (produced by hypergate())
#' @param xp The expression matrix from which the 'gate' parameter originates
#' @param gate_vector Categorical data from which the 'gate' parameter originates
#' @param level Level of gate_vector identifying the population of interest
#' @param highlight color of the positive population when plotting
#' @param path Where png files will be produced
#' @param cex size of dots
#' @param ... passed to png
#' @examples
#' data(Samusik_01_subset)
#' xp=Samusik_01_subset$xp_src[,Samusik_01_subset$regular_channels]
#' gate_vector=Samusik_01_subset$labels
#' hg=hypergate(xp=xp,gate_vector=gate_vector,level=23,delta_add=0.01)
#' par(mfrow=c(1,ceiling(length(hg$active_channels)/2)))
#' plot_gating_strategy(gate=hg,xp=xp,gate_vector=gate_vector,level=23,highlight="red")
#' @export

plot_gating_strategy<-function(gate,xp,gate_vector,level,cex=0.5,highlight="black",path="./",...){
    if(missing(path)){
        warning("path argument is missing, output won't be saved to file")
    }
    w=!is.na(gate_vector)
    xp=xp[w,,drop=F]
    gate_vector=gate_vector[w]
    truce=gate_vector==level

    parameters=gate$pars.history
    active_parameters=gate$active_channels##apply(parameters,2,function(x){x[length(x)]!=x[1]})
    parameters=parameters[,active_parameters,drop=FALSE]
    if (nrow(parameters) > 1) {
        parameters_order=
            apply(parameters,2,function(x)min(which(x!=x[1])))
        parameters=parameters[,order(parameters_order,decreasing=FALSE),drop=FALSE]
    }
    parameters=setNames(parameters[nrow(parameters),,drop=TRUE],colnames(parameters))
    channels=sub("_max","",names(parameters))
    channels=sub("_min","",channels)

    ranges.global=apply(xp[,channels,drop=F],2,range)
    rownames(ranges.global)=c("min","max")
    
    cols=rep("black",nrow(xp))
    cols[gate_vector==level]=highlight

    n=length(parameters)

    active_events=rep(T,nrow(xp))
    iter=0

    ##All parameters that are "_min" take value 2, the "_max" take value 1
    direction=rep(2,length(parameters))
    direction[grep("_max",names(parameters))]=1

    ##Loop over pairs of consecutive parameters
    for(i in seq(1,n,by=2)){
        if((i+1)<=n){
            iter=iter+1
            if(!missing(path)){
                png(paste(path,iter,".png",sep=""),...)
            }
            chan1=channels[i]
            chan2=channels[i+1]
            plot(
                xp[active_events,chan1],
                xp[active_events,chan2],
                xlab=chan1,
                ylab=chan2,
                xlim=ranges.global[,chan1],
                ylim=ranges.global[,chan2],
                bty="l",
                pch=16,
                cex=cex,
                col=cols[active_events]
            )
            segments(
                x0=parameters[i],
                y0=parameters[i+1],
                x1=ranges.global[direction[i],chan1],
                col="red"
            )
            segments(
                x0=parameters[i],
                y0=parameters[i+1],
                y1=ranges.global[direction[i+1],chan2],
                col="red"
            )

            ##Updating active_events
            if(direction[i]==2){
                test1=xp[,chan1]>=parameters[i] ##If _min, events above parameter are selected
            } else {
                test1=xp[,chan1]<=parameters[i] ##Else events above parameter below
            }
            if(direction[i+1]==2){
                test2=xp[,chan2]>=parameters[i+1]
            } else {
                test2=xp[,chan2]<=parameters[i+1]
            }
            active_events=active_events&test1&test2
            title(main=paste(paste(channels[1:(i+1)],ifelse(direction[1:(i+1)]==2,"+","-"),sep=""),collapse=", "))
            title(sub=paste("F=",signif(F_beta(truce,active_events),4),sep=""))
            if(!missing(path)){
                dev.off()
            }
        }
    }
    ##Single last parameter if n is odd
    if(n%%2==1){
        iter=iter+1
        chan1=channels[i]
        if(!missing(path)){
            png(paste(path,iter,".png",sep=""),...)
        }
        plot(xp[active_events,chan1],main="",
             ylab=channels[i],
             xlab="Events index",
             ylim=ranges.global[,chan1],
             bty="l",
             pch=16,
             cex=cex,
             col=cols[active_events]
             )
        abline(h=parameters[i],col="red")
        if(direction[i]==2){
            test1=xp[,chan1]>=parameters[i] ##If _min, events above parameter are selected
        } else {
            test1=xp[,chan1]<=parameters[i] ##Else events above parameter below
        }
        active_events=active_events&test1
        title(main=paste(paste(channels[1:(i)],ifelse(direction[1:(i)]==2,"+","-"),sep=""),collapse=", "))
        title(sub=paste("F=",signif(F_beta(truce,active_events),4),sep=""))
        if(!missing(path)){
            dev.off()
        }
    }
}

#' @title hgate_info
#' @description Extract information about a hypergate return: the channels of
#'   the phenotype, the sign of the channels, the sign of the comparison, the
#'   thresholds. The function could also compute the Fscores if the xp,
#'   gate_vector and level parameters are given.
#' @param hgate A hypergate object (produced by hypergate())
#' @param xp The expression matrix from which the 'hgate' parameter originates,
#'   needed for Fscore computation
#' @param gate_vector Categorical data from which the 'hgate' parameter
#'   originates, needed for Fscore computation
#' @param level Level of gate_vector identifying the population of interest,
#'   needed for Fscore computation
#' @param beta Beta to weight purity (low beta) or yield (high beta) more,
#'   needed for Fscore computation
#' @return A data.frame with channel, sign, comp and threshold columns, and
#'   optionnally deltaF (score deterioration when parameter is ignored),Fscore1d (F_value when using only this parameter) and Fscore (F score when all parameters up to this one are included). Fscores are computed if xp, gate_vector
#'   and level are passed to the function.
#' @seealso \code{hg_pheno}, \code{hg_rule}
#' @examples
#' data(Samusik_01_subset)
#' xp=Samusik_01_subset$xp_src[,Samusik_01_subset$regular_channels]
#' gate_vector=Samusik_01_subset$labels
#' hg=hypergate(xp=xp,gate_vector=gate_vector,level=23,delta_add=0.01)
#' hgate_info(hgate=hg)
#' hgate_pheno(hgate=hg)
#' hgate_rule(hgate=hg)
#' @export

hgate_info <- function(hgate, xp, gate_vector, level, beta = 1) {
    miss = missing(xp) + missing(level) + missing(level)
    if (miss == 1 || miss == 2) {
        warning("at least one parameter is missing in order to compute scores.")
    }
                                        # retrieve threshold
    pars = hgate$pars.history
    active_pars = hgate$active_channels
    pars = pars[, active_pars, drop = FALSE]
    if (nrow(pars) > 1) {
        pars_order = apply(pars, 2, function(x) min(which(x != x[1])))
        pars = pars[, order(pars_order, decreasing = FALSE), drop = FALSE]
    }
    pars = setNames(pars[nrow(pars), , drop = TRUE], colnames(pars))
                                        # get channel names
    channels = sub("_max", "", names(pars))
    channels = sub("_min", "", channels)
                                        # phenotype sign
    dir.sign = rep('+', length(pars))
    dir.sign[grep("_max", names(pars))] = '-'
                                        # comparison sign
    dir.comp = rep(' >= ', length(pars))
    dir.comp[grep("_max", names(pars))] = ' <= '
                                        # all together
    res = data.frame(
        channels, sign = dir.sign, comp = dir.comp, threshold = pars
    )
                                        # scores
    if (miss == 0) {
        
        w = !is.na(gate_vector)
        xp = xp[w,,drop=F]
        gate_vector = gate_vector[w]
        truce = gate_vector==level
        
        ##Loop over parameters
        active_events = rep(TRUE, nrow(xp))
        Fscore = Fscore1D = c()
        for(i in seq(length(pars))) {
                                        # events for the current 1D gate
            if(dir.comp[i] == ' >= '){
                test1D = xp[,channels[i]] >= pars[i]
            } else {
                test1D = xp[,channels[i]] <= pars[i]
            }
            Fscore1D = c(Fscore1D, signif(F_beta(truce, test1D, beta = beta), 4))
                                        # update active events
            active_events = active_events & test1D
            Fscore = c(Fscore, signif(F_beta(truce, active_events, beta = beta), 4))
        }
        res = cbind(res, deltaF=channels_contributions(hgate, xp, gate_vector, level, beta),Fscore1D, Fscore)
    }
    res
}

#' @title hgate_pheno
#' @description Build a human readable phenotype, i.e. a combination of channels
#'   and sign (+ or -) from a hypergate return.
#' @param hgate A hypergate object (produced by hypergate())
#' @param collapse A character string to separate the markers.
#' @return A string representing the phenotype.
#' @seealso \code{hg_rule}, \code{hg_info}
#' @examples
#' ## See hgate_info
#' @export

hgate_pheno <- function(hgate, collapse = ", ") {
    with(hgate_info(hgate), paste0(channels, sign, collapse = collapse))
}

#' @title hgate_rule
#' @description Build a human readable rule i.e. a combination of channels, sign
#'   of comparison and threshold.
#' @param hgate A hypergate object (produced by hypergate())
#' @param collapse A character string to separate the markers.
#' @param digits An integer that specifies the decimal part when rounding.
#' @return A data.frame with channel, sign, comp and threshold columns
#' @seealso \code{hg_pheno}, \code{hg_rule}
#' @examples
#' ## See hgate_info
#' @export

hgate_rule <- function(hgate, collapse = ", ", digits = 2) {
    with(hgate_info(hgate), paste0(channels, comp, round(threshold, digits), collapse = collapse))
}

#' @title hgate_sample
#' @description Downsample the data in order to fasten the computation and
#'   reduce the memory usage.
#' @param gate_vector A Categorical vector whose length equals the number of
#'   rows of the matrix to sample (nrow(xp))
#' @param level A level of gate_vector so that gate_vector == level will produce
#'   a boolean vector identifying events of interest
#' @param size An integer specifying the maximum number of events of interest to
#'   retain. If the count of events of interest is lower than \code{size}, than
#'   \code{size} will be set to that count.
#' @param method A string specifying the method to balance the count of events.
#'   \code{"prop"} means proportionnality: if events of interest are sampled in
#'   a 1/10 ratio, then all others events are sampled by the same ratio.
#'   \code{"10x"} means a balance of 10 between the count events of interest and
#'   the count all others events. \code{"ceil"} means a uniform sampling no more
#'   than the specified size for each level of the gate_vector. \code{level} is
#'   unused in that method.
#' @return A logical vector with TRUE correspond to the events being sampled, ie
#'   kept to further analysis
#' @note No replacement is applied. If there are less events in one group or the
#'   alternate than the algorithm requires, then all available events are
#'   returned. NA values in gate_vector are not sampled, ie ignored.
#' @examples
#' # Standard procedure with downsampling
#' data(Samusik_01_subset)
#' xp <- Samusik_01_subset$xp_src[,Samusik_01_subset$regular_channels]
#' gate_vector <- Samusik_01_subset$labels
#' sampled <- hgate_sample(gate_vector, level=8, 100)
#' table(sampled)
#' table(gate_vector[sampled])
#' xp_sampled <- xp[sampled, ]
#' gate_vector_sampled <- gate_vector[sampled]
#' hg <- hypergate(xp_sampled, gate_vector_sampled, level=8, delta_add=0.01)
#' # cluster 8 consists in 122 events
#' table(gate_vector)
#' # Downsampling
#' table(gate_vector[hgate_sample(gate_vector, level=8, 100)])
#' # Downsampling reduces the alternate events
#' table(gate_vector[hgate_sample(gate_vector, level=8, 100, "10x")])
#' # Downsampling is limited to the maximum number of events of interest
#' table(gate_vector[hgate_sample(gate_vector, level=8, 150)])
#' # Downsampling is limited to the maximum number of events of interest, and
#' # the alternate events are downsampled to a total of 10 times
#' table(gate_vector[hgate_sample(gate_vector, level=8, 150, "10x")])
#' # More details about sampling
#' # Convert -1 to NA, NA are not sampled
#' gate_vector[gate_vector==-1] = NA
#' gate_vector = factor(gate_vector)
#' table(gate_vector, useNA = "alw")
#' #
#' # target size = 100 whereas initial freq is 122 for pop 8
#' smp.prop = hgate_sample(gate_vector, level = 8, size = 100, method = "prop")
#' smp.10x  = hgate_sample(gate_vector, level = 8, size = 100, method = "10x")
#' smp.ceil = hgate_sample(gate_vector, size = 10, method = "ceil")
#' table(smp.prop)
#' table(smp.10x)
#' table(smp.ceil)
#' rbind(raw = table(gate_vector),
#'       prop = table(gate_vector[smp.prop]),
#'       `10x` = table(gate_vector[smp.10x]),
#'       ceil = table(gate_vector[smp.ceil]))
#' #
#' # target size = 30 whereas initial freq is 25 for pop 14
#' smp.prop = hgate_sample(gate_vector, level = 14, size = 30, method = "prop")
#' smp.10x  = hgate_sample(gate_vector, level = 14, size = 30, method = "10x")
#' table(smp.prop)
#' table(smp.10x)
#' rbind(raw = table(gate_vector),
#'       prop = table(gate_vector[smp.prop]),
#'       `10x` = table(gate_vector[smp.10x]))
#' # prop returns original data, because target size ids larger than initial freq
#' # 10x  returns sampled data according to initial freq, such as the total amount
#' # of other events equals 10x initial freq of pop 14
#' @export

hgate_sample <- function(gate_vector, level, size = 1000, method = "prop") {
    ## Where gate_vector is the vector of clusters and level the population of interest)
    subsample <- rep(FALSE, length(gate_vector))
                                        # multi-class methods
    if (method == "ceil") {
        if (!missing(level))
            warning(sprintf("level is ignored when method is %s", method))
        for (level in unique(gate_vector)) {
            if (is.na(level)) next()
            nna_pop <- !is.na(gate_vector)
            pos_pop <- nna_pop & (gate_vector==level)
            sum_pos <- sum(pos_pop, na.rm = TRUE)
            if (sum_pos <= size) {
                subsample[pos_pop] = TRUE
            } else {
                subsample[pos_pop][sample.int(sum_pos, size)] = TRUE
            }
        }
        return(subsample)
    }
                                        # pos vs neg methods
    nna_pop <- !is.na(gate_vector)
    pos_pop <- nna_pop & (gate_vector==level)
    sum_pos <- sum(pos_pop)
                                        # downsample positive population
    if (sum_pos <= size) {
        subsample[pos_pop] = TRUE
        pos_size = sum_pos
    } else {
        subsample[pos_pop][sample.int(sum_pos, size)] = TRUE
        pos_size = size
    }
                                        # downsample positive population
    sum_neg <- sum(nna_pop) - sum_pos
    if (method == "prop") {
        neg_size = round(pos_size / sum_pos * sum_neg)
    } else if (method == "10x") {
        neg_size = 10 * pos_size
    } else {
        stop(sprintf("method \"\" is not implemented.", method))
    }
    neg_pop <- nna_pop & (gate_vector!=level)
    if (neg_size < sum_neg) {
        idx <- sample.int(sum_neg, neg_size)
        subsample[neg_pop][idx] <- TRUE
    } else {
        subsample[neg_pop] <- TRUE
    }
    subsample
}

#' @title subset_matrix_hg
#' @description Returns a boolean vector whose TRUE elements correspond to events inside the hyperrectangle
#' @param gate a return from hypergate
#' @param xp Expression matrix used for gate
#' @examples
#' data(Samusik_01_subset)
#' xp=Samusik_01_subset$xp_src[,Samusik_01_subset$regular_channels]
#' gate_vector=Samusik_01_subset$labels
#' hg=hypergate(xp=xp,gate_vector=gate_vector,level=23,delta_add=0.01)
#' gating_state=subset_matrix_hg(hg,xp)
#' gating_state=ifelse(gating_state,"Gated in","Gated out")
#' target=ifelse(gate_vector==23,"Target events","Others")
#' table(gating_state,target)
#' @export
subset_matrix_hg<-function(gate,xp){
    if(length(gate$active_channels)>0){
        state=rep(TRUE,nrow(xp))
        pars=tail(gate$pars.history,1)[1,]
        for(chan in gate$active_channels){
            chan_real=substr(chan,1,nchar(chan)-4)
            if(substr(chan,nchar(chan)-3,nchar(chan))=="_min"){
                state[state][xp[state,chan_real]<pars[chan]]=FALSE
            } else {
                state[state][xp[state,chan_real]>pars[chan]]=FALSE
            }        
        }
        return(state)
    } else {
        return(rep(TRUE,nrow(xp)))
    }
}

#' @title contract
#' @description Test (some) possible contractions of the hyperrectangle
#' @param par Current parametrization of the hyperrectangle
#' @param xp_pos Expression matrix for positive events
#' @param state_pos State vector of the positive events
#' @param xp_neg Expression matrix for negative events
#' @param state_neg State vector of the negative events
#' @param n passed to f
#' @param TP integer: current number of TP
#' @param TN integer: current number of TN
#' @param beta Passed from the top-level function
#' @param envir Current environment of the optimization
contract<-function(
                   par=par,
                   xp_pos=envir$xp_pos,
                   state_pos=envir$state_pos,
                   xp_neg=envir$xp_neg,
                   state_neg=envir$state_neg,
                   n=envir$n,
                   TP=envir$TP,
                   TN=envir$TN,
                   beta=envir$beta2,
                   envir=parent.frame()
                   ){
    ##Evaluate what happens when we increase (strictly) the parameter to fit a TP
    TP_par=sort(xp_pos[state_pos,par]) ##Sorted values, possibly duplicates
    FP_par=sort(xp_neg[state_neg,par]) ##Sorted values, possibly duplicates

    cc=unique(TP_par) #Contraction candidates.
    if(length(cc)<1){
        return(
            list(
                newF=envir$f_current,
                par_best=setNames(min(envir$hg.env$xp[,par]),par),
                flag_add=FALSE,
                dTP_max=0,
                dTN_max=0
            )
        )
    }
    i_TP=1L ##Iterator for TP vector
    i_FP=1L ##Iterator for FP vector
    dTP=0L ##Keep track of changes in TP for every TP cutoff value
    dTN=0L ##Keep track of changes in TN for every TP cutoff value

    ncc=length(cc)
    dTPvec=rep(0L,ncc)
    dTNvec=rep(0L,ncc)
    j=1L
    for(level in cc){
        ##How many TP we lose if we cut just below this level
        while(TP_par[i_TP]<level&i_TP<=length(TP_par)){
            dTP=dTP-1L
            i_TP=i_TP+1L
        }
        ##How many TN we gain if we cut just below this level
        while(FP_par[i_FP]<level&i_FP<=length(FP_par)){
            dTN=dTN+1L
            i_FP=i_FP+1L
        }

        dTPvec[j]=dTP
        dTNvec[j]=dTN
        j=j+1L
    }
    
    newF=f(rep(TP,ncc)+dTPvec,rep(TN,ncc)+dTNvec,n,beta2=beta)
    
    w=which.max(newF)
    
    list(
        newF=newF[w],
        par_best=setNames(cc[w],par),
        flag_add=FALSE,
        dTP_max=dTPvec[w],
        dTN_max=dTNvec[w]
    )
}

##From current state and update parameters, returns updated state
#' @title contract.update
#' @description Update the hyperrectangle to the best contraction move found
#' @param contract_object output of the contract function
#' @param pars Current parametrization of the hyperrectangle
#' @param active_channels vector of currently-used parameters
#' @param b_pos boolean matrix of positive events
#' @param b_neg boolean matrix of negative events
#' @param xp_pos Expression matrix for positive events
#' @param state_pos State vector of the positive events
#' @param xp_neg Expression matrix for negative events
#' @param state_neg State vector of the negative events
#' @param envir Current environment of the optimization
#' @param TP integer: current number of TP
#' @param TN integer: current number of TN
contract.update<-function(
                          contract_object,
                          pars=envir$pars,
                          active_channels=envir$active_channels,
                          b_pos=envir$b_pos,
                          b_neg=envir$b_neg,
                          state_pos=envir$state_pos,
                          state_neg=envir$state_neg,
                          TN=envir$TN,
                          TP=envir$TP,
                          xp_pos=envir$xp_pos,
                          xp_neg=envir$xp_neg,
                          envir=parent.frame()
                          ){
    pars[names(contract_object$par_best)]=contract_object$par_best
    if(!names(contract_object$par_best)%in%active_channels){
        active_channels=c(active_channels,names(contract_object$par_best))
    }

    ##Values in B that may be affected (previously TRUE)
    w_pos=b_pos[,names(contract_object$par_best)]
    w_neg=b_neg[,names(contract_object$par_best)]

    ##Those that are indeed affected in B
    w_state_pos=xp_pos[w_pos,names(contract_object$par_best)]<contract_object$par_best
    w_state_neg=xp_neg[w_neg,names(contract_object$par_best)]<contract_object$par_best

    ##We update the state
    state_pos[w_pos][w_state_pos]=FALSE
    state_neg[w_neg][w_state_neg]=FALSE

    ##We perform the changes in B
    b_pos[w_pos,names(contract_object$par_best)][w_state_pos]=FALSE
    b_neg[w_neg,names(contract_object$par_best)][w_state_neg]=FALSE

    list(
        pars=pars,
        active_channels=active_channels,
        state_pos=state_pos,
        state_neg=state_neg,
        b_pos=b_pos,
        b_neg=b_neg,
        TP=TP+contract_object$dTP,
        TN=TN+contract_object$dTN
    )
}

#' @title expand
#' @description Test (some) possible expansions of the hyperrectangle
#' @param n passed to f
#' @param TP integer: current number of TP
#' @param TN integer: current number of TN
#' @param FN integer: current number of FP
#' @param beta Passed from the top-level function 
#' @param envir Coreloop environment
#' @param FNTN_matrix Boolean matrix of dim (FN, FN + TN), where Mij is TRUE if and only if expanding to include the ith FN in the gate would lead to the inclusion of the jth column event
expand<-function(
                 FN=envir$FN,
                 FNTN_matrix=envir$FNTN_matrix,
                 TP=envir$TP,
                 TN=envir$TN,
                 n=envir$n,
                 beta=envir$beta2,
                 envir=parent.frame()
                 ){
    ##Expanding only concerns events that are currently !state. Events for which "state" is positive are guaranteed to stay "state==TRUE".
    ##FNTN_matrix tries to create expansion that match a given event, and flag !state events that will be included along with it

    dTP=colSums(FNTN_matrix[1:FN,1:FN,drop=FALSE])
    dTN=-colSums(FNTN_matrix[(FN+1):nrow(FNTN_matrix),1:FN,drop=FALSE])
    TP.vector=rep(TP,FN)
    TN.vector=rep(TN,FN)
    n.vector=rep(n,FN)
    dF_expansion=f(TP.vector+dTP,TN.vector+dTN,n.vector,beta2=beta)
    w=which.max(dF_expansion)
    newF=dF_expansion[w]
    list(
        newF=newF,
        dTP_max=dTP[w],
        dTN_max=dTN[w],
        which_expansion=w,
        FN=FN,
        FNTN_matrix=FNTN_matrix
    )
}

#' @title expand.update
#' @description Update the hyperrectangle to the best expansion move found
#' @param expand.object output of the expand function
#' @param pars Current parametrization of the hyperrectangle
#' @param xp_pos Expression matrix for positive events
#' @param xp_neg Expression matrix for negative events
#' @param state_pos State vector of the positive events
#' @param state_neg State vector of the negative events
#' @param b_pos boolean matrix of positive events
#' @param b_neg boolean matrix of negative events
#' @param n passed to f
#' @param TP integer: current number of TP
#' @param TN integer: current number of TN
#' @param envir Current environment of the optimization
expand.update<-function(
                        expand.object,
                        pars=envir$pars,
                        xp_pos=envir$xp_pos,
                        xp_neg=envir$xp_neg,
                        state_pos=envir$state_pos,
                        state_neg=envir$state_neg,
                        b_pos=envir$b_pos,
                        b_neg=envir$b_neg,
                        n=envir$n,
                        TP=envir$TP,
                        TN=envir$TN,
                        envir=parent.frame()
                        ){

    f_best=expand.object$newF
    expansion_parameters=setNames(xp_pos[!state_pos,,drop=FALSE][expand.object$which_expansion,],names(pars))
    affected_expansion_parameters=expansion_parameters<pars
    par_best=expansion_parameters[affected_expansion_parameters]

    pars[names(par_best)]=par_best

    ##For recycling
    B_FN_old=b_pos
    B_TN_old=b_neg
    
    for(p in which(affected_expansion_parameters)){
        b_pos[,p]=xp_pos[,p]>=pars[p]
        b_neg[,p]=xp_neg[,p]>=pars[p]
    }
    
    w_FN=expand.object$FNTN_matrix[1:expand.object$FN,expand.object$which_expansion]
    w_TN=expand.object$FNTN_matrix[(expand.object$FN+1):nrow(expand.object$FNTN_matrix),expand.object$which_expansion]
    
    state_pos[!state_pos][w_FN]=TRUE ##FN converted to TP
    state_neg[!state_neg][w_TN]=TRUE ##TN converted to FP

    ##For recycling
    B_FN_old=B_FN_old[!state_pos,,drop=FALSE]
    B_TN_old=B_TN_old[!state_neg,,drop=FALSE]
    B_FN_new=b_pos[!state_pos,,drop=FALSE]
    B_TN_new=b_neg[!state_neg,,drop=FALSE]
    xp_FN=xp_pos[!state_pos,,drop=FALSE]
    xp_TN=xp_neg[!state_neg,,drop=FALSE]
    
    TN=TN+expand.object$dTN_max
    TP=TP+expand.object$dTP_max

    ##f_current=f(TP,TN,n,beta2=beta)

    flag_expansion=TP<length(state_pos) ##If the expansion included all positive events, cannot try to expand more. Otherwise it's TRUE (so we try expanding again)

    FNTN_matrix=FNTN_matrix.recycle(
        expand.object$FNTN_matrix[!c(w_FN,w_TN),!w_FN,drop=FALSE],
        B_FN_old,
        B_TN_old,
        B_FN_new,
        B_TN_new,
        xp_FN,
        xp_TN,
        pars[envir$active_channels]
    )
    ## ##TMP failcheck
    ## if(TN!=sum(!state_neg)){
    ##     stop("Problem TN mismatch")
    ## }
    ## if(TP!=sum(state_pos)){
    ##     stop("Problem TP mismatch")
    ## }
    ## if(f_current!=f_best){
    ##     stop("f_current != f_best in contraction")
    ## }
    ## if(any(state_pos!=(rowSums(b_pos)==ncol(b_pos)))){
    ##     stop("Mismatch state_pos b_pos")
    ## }
    ## if(any(state_neg!=(rowSums(b_neg)==ncol(b_neg)))){
    ##     stop("Mismatch state_neg b_neg")
    ## }

    return(
        list(
            pars=pars,
            state_pos=state_pos,
            state_neg=state_neg,
            b_pos=b_pos,
            b_neg=b_neg,
            flag_expansion=flag_expansion,
            TN=TN,
            TP=TP,
            par_best=par_best,
            FNTN_matrix=FNTN_matrix
        )
    )
}

#' @title FNTN_matrix.recycle
#' @description Recycle an expansion matrix
#' @param FNTN_matrix Expansion matrix to recycle
#' @param xp_FN Expression matrix of False Negative events
#' @param xp_TN Expression matrix of True Negative events
#' @param B_FN_old Boolean matrix of FN events before the last expansion
#' @param B_TN_old Boolean matrix of TN events before the last expansion
#' @param B_FN_new Boolean matrix of FN events after the last expansion
#' @param B_TN_new Boolean matrix of TN events after the last expansion
#' @param par Current hyper-rectangle parametrization
FNTN_matrix.recycle<-function(
                              FNTN_matrix, ##element to recycle
                              B_FN_old,
                              B_TN_old,
                              B_FN_new,
                              B_TN_new,
                              xp_FN,
                              xp_TN,
                              par
                              ){
    
    ##Identify elements that changed
    ##w_FN_changed=rowSums(xor(B_FN_old,B_FN_new))>0
    w_FN_changed=rep(TRUE,ncol(FNTN_matrix)) ##We have to include every FN whether they changed or not
    w_TN_changed=rowSums(xor(B_TN_old[,names(par)],B_TN_new[,names(par)]))>0

    if(any(w_FN_changed)){
        FNTN_matrix[c(w_FN_changed,w_TN_changed),w_FN_changed]=fill_FNTN_matrix(xp_FN[w_FN_changed,names(par),drop=FALSE],xp_TN[w_TN_changed,names(par),drop=FALSE],B_FN_new[w_FN_changed,names(par),drop=FALSE],B_TN_new[w_TN_changed,names(par),drop=FALSE],par)
    }
    FNTN_matrix
}

#' @title coreloop
#' @param par Current parametrization of the hyperrectangle
#' @param hg.env Environment where the main execution of hypergate takes place
#' @description Core optimization loop of hypergate
coreloop<-function(par,hg.env=hg.env$hg.env){
    loop.env=environment()
    pars=hg.env$par
    active_channels=hg.env$active_channels
    TP=hg.env$TP
    state_pos=hg.env$state_pos
    state_neg=hg.env$state_neg
    b_pos=hg.env$b_pos
    b_neg=hg.env$b_neg
    
    
    ##Step 1 ADD NEW CHANNEL
    contractions=contract(loop.env$par,envir=loop.env$hg.env)
    ##best_contraction=which.max(sapply(contractions,function(x)x$newF))
    if(contractions$newF>loop.env$hg.env$f_current){
        contraction=contract.update(contractions,envir=loop.env$hg.env)
        
        sapply(names(contraction),function(x)assign(x,contraction[[x]],envir=loop.env))
        f_current=contractions$newF
        
        ##Update add channel
        if(hg.env$verbose){
            print(paste("Trying to add",names(contractions$par_best),":",contractions$par_best,",F=",signif(f_current,4L)))
        }
        f_res=c(loop.env$hg.env$f_res,f_current)
        pars.history.rank=rbind(loop.env$hg.env$pars.history.rank,pars)
        flag_expansion=TRUE
        flag_contraction=TRUE
    } else {
        if(hg.env$verbose){
            print("No improvement")
        }
        assign("f_current",hg.env$f_current,envir=loop.env)
        return(loop.env)
    }
    
    ##Once we have added a channel
    ##We expand as much as possible, and then try to contract
    ##We stop only if we haven't done anything (so no contraction or expansion is possible)
    while(flag_expansion|flag_contraction){
        ##EXPANSION
        flag_expansion=length(active_channels)>1&TP<length(state_pos) ##We only try to expand if we already have two channels... and its not possible if we don't have a FN
        flag_contraction=length(active_channels)>1
        flag_recycling=FALSE
        while(flag_expansion){
            FN=length(state_pos)-TP
            ##Expanding only concerns events that are currently !state. Events for which "state" is positive are guaranteed to stay "state==TRUE".
            ##FNTN_matrix tries to create expansion that match a given event, and flag !state events that will be included along with it
            if(!flag_recycling){
                FNTN_matrix=fill_FNTN_matrix(loop.env$hg.env$xp_pos[!state_pos,active_channels,drop=FALSE],loop.env$hg.env$xp_neg[!state_neg,active_channels,drop=FALSE],b_pos[!state_pos,active_channels,drop=FALSE],b_neg[!state_neg,active_channels,drop=FALSE],pars[active_channels])
            } else {
                FNTN_matrix=expansion$FNTN_matrix
                ## FNTN_matrix.raw=fill_FNTN_matrix(xp_pos[!state_pos,active_channels,drop=FALSE],xp_neg[!state_neg,active_channels,drop=FALSE],b_pos[!state_pos,active_channels,drop=FALSE],b_neg[!state_neg,active_channels,drop=FALSE],pars[active_channels])
                ## if(any(FNTN_matrix.raw!=expansion$FNTN_matrix)){
                ##     recover()
                ## }
            }
            expansions=expand(envir=loop.env,n=hg.env$n,beta=hg.env$beta2)
            if(expansions$newF>f_current){
                which_expansion=expansions$which_expansion
                expansion=expand.update(expansions,envir=loop.env,xp_pos=loop.env$hg.env$xp_pos,xp_neg=loop.env$hg.env$xp_neg)
                flag_expansion=TP<length(state_pos) ##If the expansion included all positive events, cannot try to expand more.
                sapply(names(expansion),function(x)assign(x,expansion[[x]],envir=loop.env))
                f_current=expansions$newF
                f_res=c(f_res,f_current)
                pars.history.rank=rbind(pars.history.rank,pars)
                
                if(hg.env$verbose){
                    print(paste("Expansion of",paste(names(expansion$par_best),":",expansion$par_best,", ",collapse=""),"F=",signif(f_current,4)))
                }
                flag_recycling=TRUE
            } else {
                flag_expansion=FALSE
            }
        }

        ##CONTRACTION
        contractions=sapply(active_channels,contract,simplify=FALSE,envir=loop.env,xp_pos=hg.env$xp_pos,xp_neg=hg.env$xp_neg,n=hg.env$n,beta=hg.env$beta2)
        best_contraction=which.max(sapply(contractions,function(x)x$newF))
        if(contractions[[best_contraction]]$newF>f_current){
            contraction=contract.update(contractions[[best_contraction]],envir=loop.env,xp_pos=hg.env$xp_pos,xp_neg=hg.env$xp_neg)
            sapply(names(contraction),function(x)assign(x,contraction[[x]],envir=loop.env))
            f_current=contractions[[best_contraction]]$newF
            if(hg.env$verbose){
                print(paste("Contracting ",names(contractions[[best_contraction]]$par_best),":",contractions[[best_contraction]]$par_best,",F=",signif(f_current,4L)))
            }
            f_res=c(f_res,f_current)
            pars.history.rank=rbind(pars.history.rank,pars)
            flag_contraction=TRUE
        } else {
            flag_contraction=FALSE
        }
    }
    return(loop.env)
}

#' @title hypergate
#' @description Finds a hyperrectangle gating around a population of interest
#' @param xp an Expression matrix
#' @param gate_vector A Categorical vector of length nrow(xp)
#' @param level A level of gate_vector so that gate_vector == level will produce a boolean vector identifying events of interest
#' @param delta_add If the increase in F after an optimization loop is lower than delta_add, the optimization will stop (may save computation time)
#' @param beta Purity / Yield trade-off
#' @param verbose Boolean. Whether to print information about the optimization status.
#' @examples
#' data(Samusik_01_subset)
#' xp=Samusik_01_subset$xp_src[,Samusik_01_subset$regular_channels]
#' gate_vector=Samusik_01_subset$labels
#' hg=hypergate(xp=xp,gate_vector=gate_vector,level=23,delta_add=0.01)
#' @seealso \code{\link{channels_contributions}} for ranking parameters within the output, \code{\link{reoptimize_strategy}} for reoptimizing a output on a subset of the markers, \code{\link{plot_gating_strategy}} for plotting an output, \code{\link{subset_matrix_hg}} to apply the output to another input matrix, \code{\link{boolmat}} to obtain a boolean matrix stating which events are filtered out because of which markers
#' @export

hypergate<-function(xp,gate_vector,level,delta_add=0,beta=1,verbose=FALSE){
    beta2=beta^2
    if(is.null(rownames(xp))){
        rownames(xp)=1:nrow(xp)
    }
    r=setNames(rownames(xp),1:nrow(xp))
    rownames(xp)=1:nrow(xp)

    xp_src=xp

    ##Sorting values so that cutoffs are very easy to compute. Rank = 1 is small, rank = n is big
    xp=apply(xp,2,rank,ties.method="min")
    rownames(xp)=1:nrow(xp)

    ##Each parameter is duplicated (for upper and lower cutoffs)
    ##n=nrow(xp)
    colnames=paste(rep(colnames(xp),each=2),c("min","max"),sep="_")

    xp.2=matrix(nrow=nrow(xp),ncol=2*ncol(xp),data=0,dimnames=list(rownames(xp),colnames))
    storage.mode(xp.2)="integer"

    ##For both mins (odd columns) and maxs (even columns), make sure than low rank = "extreme" values (likely to pop). We can then use min to select
    xp.2[,seq(1,ncol(xp.2),by=2)]=xp
    xp.2[,seq(2,ncol(xp.2),by=2)]=nrow(xp)-xp
    xp=xp.2
    rm(xp.2)

    truce=gate_vector==level
    xp_pos=xp[truce,,drop=FALSE] ##xp for positive elements
    state_pos=rep(T,nrow(xp_pos)) ##State vector (event selected or not in current gate)
    b_pos=matrix(nrow=nrow(xp_pos),ncol=ncol(xp_pos),dimnames=dimnames(xp_pos),data=TRUE) ##State vector for every channel. Keeps track of whether the event is above the cut-off value for this parameter (parameter = (channel,upper/lower bound)

    xp_neg=xp[!truce,,drop=FALSE] ##xp for negative elements
    state_neg=rep(T,nrow(xp_neg))
    b_neg=matrix(nrow=nrow(xp_neg),ncol=ncol(xp_neg),dimnames=dimnames(xp_neg),data=TRUE)

    ##Number of TN and TP to compute F
    TN=sum(!state_neg)
    TP=sum(state_pos)
    n=beta2*length(state_pos)+length(state_neg) ##This is only used in the function that computes F_beta.
    
    pars=setNames(as.integer(apply(xp,2,min)),colnames(xp)) ##Rank-parameters initialization
    pars.history.rank=matrix(nrow=1,ncol=length(pars),dimnames=list("0",names(pars)),data=pars) ##For reporting purposes we keep track of every optimization step
    storage.mode(pars.history.rank)="integer"

    active_channels=character()

    ##Optimization loop

    ##f value
    f_current=0 ##This is wrong but the optimization sometimes gets stuck really early if set otherwise (i.e. big cluster will naturally have "high" F if all==TRUE)
    f_res=f_current ##Just for reporting purposes

    hg.env=environment()

    while(length(active_channels)<=ncol(xp)){ ##This condition just makes sure that we can't add more channels than possible. The true ending condition is at the very end of the loop and enforced using "break"
        f_previous=f_current

        ##Core loop
        cycle=new.env()
        assign("f_current",f_current,envir=cycle)
        for(par in setdiff(colnames(xp),active_channels)){
            current_cycle=coreloop(par,hg.env=hg.env)
            if(current_cycle$f_current>cycle$f_current){
                cycle=current_cycle
                rm(current_cycle)
            }
        }
        
        if(cycle$f_current==f_previous){
            if(verbose){
                print("Found no channel to add. Exiting")
            }
            break
        }
        
        sapply(ls(envir=cycle),function(x){
            assign(x,envir=hg.env,value=get(x,envir=cycle))
        })
        
        ##Ending condition: last channel did not bring much.
        if((f_current-f_previous)<delta_add){
            if(verbose){
                print("Increase in f lower than delta_add, exiting")
            }
            if(length(active_channels)>1){ ##Unless its the only channel
                pop=which(pars.history.rank[,tail(active_channels,1)]!=pars.history.rank[1,tail(active_channels,1)])[1] ##Flag history from when we added this last channel...
                pars.history.rank=pars.history.rank[1:(pop-1),] ##... and remove it
                f_res=f_res[1:(pop-1)]
                active_channels=active_channels[-length(active_channels)]
            }
            break
        }
    }

    ##RETURN
    pars.history=pars.history.rank
    pars.history=matrix(nrow=nrow(pars.history.rank),ncol=ncol(pars.history.rank),data=NA,dimnames=dimnames(pars.history.rank))
    for(j in 1:ncol(pars.history)){
        pars.history[,j]=xp_src[match(pars.history.rank[,j],xp[,j]),j%%2+j%/%2]
    }
    return(list(pars.history=pars.history,pars.history.rank=pars.history.rank,f=f_res,active_channels=active_channels))
}

#' @title channels_contributions
#' @description Gives scores for the contribution of individual channels to a gating strategy
#' @param gate A return from hypergate
#' @param xp Expression matrix as in the hypergate call
#' @param gate_vector Categorical vector of length nrow(xp)
#' @param level A level of gate_vector that identifies the population of interest
#' @param beta, should be the same as for the hypergate object
#' @examples
#' data(Samusik_01_subset)
#' xp=Samusik_01_subset$xp_src[,Samusik_01_subset$regular_channels]
#' gate_vector=Samusik_01_subset$labels
#' hg=hypergate(xp=xp,gate_vector=gate_vector,level=23,delta_add=0)
#' contribs=channels_contributions(gate=hg,xp=xp,gate_vector=gate_vector,level=23,beta=1)
#' contribs
#' @export
channels_contributions<-function(gate,xp,gate_vector,level,beta=1){
    truce=gate_vector==level
    F_beta(subset_matrix_hg(gate,xp),truce,beta)-sapply(gate$active_channels,function(chan){
        subchans=setdiff(gate$active_channels,chan)
        gate$active_channels=subchans
        pred=subset_matrix_hg(gate,xp)
        F_beta(pred,truce,beta)
    })
}

#' @title boolmat
#' @description Convert an expression matrix and a gating strategy to a boolean matrix (whether each event is gated out by each channel)
#' @param gate A return from hypergate
#' @param xp Expression matrix as in the hypergate callxp=Samusik_01_subset$xp_src[,Samusik_01_subset$regular_channels]
#' @examples
#' data(Samusik_01_subset)
#' xp=Samusik_01_subset$xp_src
#' gate_vector=Samusik_01_subset$labels
#' hg=hypergate(xp=xp,gate_vector=gate_vector,level=23,delta_add=0.01)
#' head(boolmat(hg,xp))
#' @export
boolmat<-function(gate,xp){
    chans=gate$active_channels
    if(length(chans)>0){
        cols=colnames(xp)
        signs=rep(0,2*ncol(xp))
        signs[2*which(paste(cols,"_min",sep="")%in%chans)-1]=1
        signs[2*which(paste(cols,"_max",sep="")%in%chans)]=-1
        ui=matrix(nrow=2*ncol(xp),ncol=ncol(xp),data=0)
        for(i in 1:length(signs)){
            ui[i,ceiling(i/2)]=signs[i]
        }
        ci=tail(gate$pars.history,1)[1,]*-signs
        ci=matrix(nrow=ncol(xp)*2,ncol=nrow(xp),data=rep(ci,nrow(xp)),byrow=FALSE)
        res=t((ui%*%t(xp)+ci)>=0)
        colnames(res)=colnames(gate$pars.history)
        return(res[,chans])
    }
}

#' @title reoptimize_strategy
#' @description Optimize a gating strategy given a manual selection of channels
#' @param gate A return from hypergate
#' @param channels_subset Character vector identifying the channels that will be retained (others are ignored). The form is e.g. c("CD4_min","CD8_max")
#' @param xp Expression matrix as in the hypergate call
#' @param gate_vector Categorical vector as in the hypergate call
#' @param level Level of gate_vector identifying the population of interest
#' @param beta Yield / purity trade-off
#' @param verbose Whether to print information about optimization status
#' @examples
#' data(Samusik_01_subset)
#' xp=Samusik_01_subset$xp_src[,Samusik_01_subset$regular_channels]
#' gate_vector=Samusik_01_subset$labels
#' hg=hypergate(xp=xp,gate_vector=gate_vector,level=23,delta_add=0)
#' contribs=channels_contributions(gate=hg,xp=xp,gate_vector=gate_vector,level=23,beta=1)
#' significant_channels=names(contribs)[contribs>=0.01]
#' hg_reoptimized=reoptimize_strategy(gate=hg,channels_subset=significant_channels,xp,gate_vector,23)
#' @export
reoptimize_strategy<-function(gate,channels_subset,xp,gate_vector,level,beta=1,verbose=FALSE){
    beta2=beta^2
    gate$active_channels=channels_subset
    if(is.null(rownames(xp))){
        rownames(xp)=1:nrow(xp)
    }
    r=setNames(rownames(xp),1:nrow(xp))
    rownames(xp)=1:nrow(xp)
    
    xp_src=xp
    
    ##Sorting values so that cutoffs are very easy to compute. Rank = 1 is small, rank = n is big
    xp=apply(xp,2,rank,ties.method="min")
    rownames(xp)=1:nrow(xp)

    ##Each parameter is duplicated (for upper and lower cutoffs)
    ##n=nrow(xp)
    colnames=paste(rep(colnames(xp),each=2),c("min","max"),sep="_")
    
    xp.2=matrix(nrow=nrow(xp),ncol=2*ncol(xp),data=0,dimnames=list(rownames(xp),colnames))
    storage.mode(xp.2)="integer"

    ##For both mins (odd columns) and maxs (even columns), make sure than low rank = "extreme" values (likely to pop). We can then use min to select
    xp.2[,seq(1,ncol(xp.2),by=2)]=xp
    xp.2[,seq(2,ncol(xp.2),by=2)]=nrow(xp)-xp
    xp=xp.2
    rm(xp.2)

    b_src=matrix(nrow=nrow(xp),ncol=ncol(xp),data=TRUE,dimnames=dimnames(xp))
    b_src.tmp=boolmat(gate,xp_src)
    b_src[,colnames(b_src.tmp)]=b_src.tmp
    rm(b_src.tmp)
    
    truce=gate_vector==level
    xp_pos=xp[truce,,drop=FALSE] ##xp for positive elements
    b_pos=b_src[truce,,drop=FALSE]
    state_pos=rowSums(b_pos)==ncol(b_pos)
    
    xp_neg=xp[!truce,,drop=FALSE] ##xp for negative elements
    b_neg=b_src[!truce,,drop=FALSE]
    state_neg=rowSums(b_neg)==ncol(b_neg)

    
    ##Number of TN and TP to compute F
    TN=sum(!state_neg)
    TP=sum(state_pos)
    n=beta2*length(state_pos)+length(state_neg) ##This is only used in the function that computes F_beta.

    pars.history=gate$pars.history
    pars.history.tmp=tail(pars.history,1)
    pars.history.tmp[,setdiff(colnames(pars.history.tmp),gate$active_channels)]=pars.history[1,setdiff(colnames(pars.history.tmp),gate$active_channels)]
    pars.history=pars.history.tmp
    rm(pars.history.tmp)

    pars.history.rank=gate$pars.history.rank
    pars.history.rank.tmp=tail(pars.history.rank,1)
    pars.history.rank.tmp[,setdiff(colnames(pars.history.rank.tmp),gate$active_channels)]=pars.history.rank[1,setdiff(colnames(pars.history.rank.tmp),gate$active_channels)]
    pars.history.rank=pars.history.rank.tmp
    rm(pars.history.rank.tmp)
    storage.mode(pars.history.rank)="integer"

    pars=pars.history.rank[1,]    
    active_channels=gate$active_channels

    ##Optimization loop   
    ##f value
    f_current=f(TP,TN,n,beta2)
    f_res=f_current
    
    hg.env=environment()
    loop.env=environment()

    flag_expansion=TRUE
    flag_contraction=TRUE
    while(flag_expansion|flag_contraction){
        ##EXPANSION
        flag_recycling=FALSE
        flag_expansion=length(state_pos)>TP
        while(flag_expansion){
            ##Expanding only concerns events that are currently !state. Events for which "state" is positive are guaranteed to stay "state==TRUE".
            ##FNTN_matrix tries to create expansion that match a given event, and flag !state events that will be included along with it
            FN=length(state_pos)-TP
            if(!flag_recycling){
                FNTN_matrix=fill_FNTN_matrix(
                    loop.env$hg.env$xp_pos[!state_pos,active_channels,drop=FALSE],
                    loop.env$hg.env$xp_neg[!state_neg,active_channels,drop=FALSE],
                    b_pos[!state_pos,active_channels,drop=FALSE],
                    b_neg[!state_neg,active_channels,drop=FALSE],
                    pars[active_channels]
                )
            } else {
                FNTN_matrix=expansion$FNTN_matrix
            }
            expansions=expand(envir=loop.env,n=hg.env$n,beta=hg.env$beta2)
            if(expansions$newF>f_current){
                which_expansion=expansions$which_expansion
                expansion=expand.update(expansions,envir=loop.env,xp_pos=loop.env$hg.env$xp_pos,xp_neg=loop.env$hg.env$xp_neg)
                FN=length(state_pos)-TP
                flag_expansion=FN>0 ##If the expansion included all positive events, cannot try to expand more.
                sapply(names(expansion),function(x)assign(x,expansion[[x]],envir=loop.env))
                f_current=expansions$newF
                f_res=c(f_res,f_current)
                pars.history.rank=rbind(pars.history.rank,pars)
                
                if(hg.env$verbose){
                    print(paste("Expansion of",paste(names(expansion$par_best),":",expansion$par_best,", ",collapse=""),"F=",signif(f_current,4)))
                }
                flag_recycling=TRUE
            } else {
                flag_expansion=FALSE
            }
        }

        ##CONTRACTION
        contractions=sapply(active_channels,contract,simplify=FALSE,envir=loop.env,xp_pos=hg.env$xp_pos,xp_neg=hg.env$xp_neg,n=hg.env$n,beta=hg.env$beta2)
        best_contraction=which.max(sapply(contractions,function(x)x$newF))
        if(contractions[[best_contraction]]$newF>f_current){
            contraction=contract.update(contractions[[best_contraction]],envir=loop.env,xp_pos=hg.env$xp_pos,xp_neg=hg.env$xp_neg)
            sapply(names(contraction),function(x)assign(x,contraction[[x]],envir=loop.env))
            f_current=contractions[[best_contraction]]$newF
            if(hg.env$verbose){
                print(paste("Contracting ",names(contractions[[best_contraction]]$par_best),":",contractions[[best_contraction]]$par_best,",F=",signif(f_current,4L)))
            }
            f_res=c(f_res,f_current)
            pars.history.rank=rbind(pars.history.rank,pars)
            flag_contraction=TRUE
        } else {
            flag_contraction=FALSE
        }
    }

    ##RETURN
    pars.history=pars.history.rank
    pars.history=matrix(nrow=nrow(pars.history.rank),ncol=ncol(pars.history.rank),data=NA,dimnames=dimnames(pars.history.rank))
    for(j in 1:ncol(pars.history)){
        pars.history[,j]=xp_src[match(pars.history.rank[,j],xp[,j]),j%%2+j%/%2]
    }
    return(list(pars.history=pars.history,pars.history.rank=pars.history.rank,f=f_res,active_channels=active_channels))
}

#' @title F_beta
#' @description Compute a F_beta score comparing two boolean vectors
#' @param pred boolean vector of predicted values
#' @param truth boolean vector of true values
#' @param beta Weighting of yield as compared to precision. Increase beta so that the optimization favors yield, or decrease to favor purity.
#' @examples
#' data(Samusik_01_subset)
#' truth=c(rep(TRUE,40),rep(FALSE,60))
#' pred=rep(c(TRUE,FALSE),50)
#' table(pred,truth) ##40% purity, 50% yield
#' #' F_beta(pred=pred,truth=truth,beta=2) ##Closer to yield
#' F_beta(pred=pred,truth=truth,beta=1.5) ##Closer to yield
#' F_beta(pred=pred,truth=truth,beta=1) ##Harmonic mean
#' F_beta(pred=pred,truth=truth,beta=0.75) ##Closer to purity
#' F_beta(pred=pred,truth=truth,beta=0.5) ##Closer to purity
#' @export

F_beta=function(pred,truth,beta=1){
    TP=sum(truth&pred)
    if(TP==0){
        return(0)
    }
    FP=sum(!truth&pred)
    FN=sum(truth&!pred)
    yield=TP/(TP+FN) #aka recall
    purity=TP/(TP+FP) #aka precision
    if(is.nan(purity)|is.nan(yield)){
        return(0)
    }
    F=(1+beta^2)*(yield*purity)/(yield+beta^2*purity)
    F
}

#' @title gate_from_biplot
#' @description From a biplot let the user interactively draw polygons to create a "Gate" vector
#' @param matrix A matrix
#' @param x_axis character, colname of matrix used for x-axis in the biplot
#' @param y_axis character, colname of matrix used for y-axis in the biplot
#' @param sample Used to downsample the data in case there are too many events to plot quickly
#' @param bty passed to plot
#' @param cex passed to plot
#' @param pch passed to plot
#' @param ... passed to plot
#' @examples
#' if(interactive()){
#'     ##See the details section to see how this function works
#'     gate_from_biplot(matrix=Samusik_01_subset$tsne,x_axis="tSNE1",y_axis="tSNE2")
#' }
#' @export
#' @return A named vector of length nrow(matrix) and names rownames(matrix). Ungated events are set to NA
#' @details Data will be displayed as a bi-plot according to user-specified x_axis and y_axis arguments, then a call to locator() is made. The user can draw a polygon around parts of the plot that need gating. When done, 'right-click' or 'escape' (depending on the IDE) escapes locator() and closes the polygon. Then the user can press "n" to draw another polygon (that will define a new population), "c" to cancell and draw the last polygon again, or "s" to exit. When exiting, events that do not fall within any polygon are assigned NA, the others are assigned an integer value corresponding to the last polygon they lie into.

gate_from_biplot<-function(matrix,x_axis,y_axis,...,bty="l",pch=16,cex=0.5,sample=NULL)
{
    xp=matrix[,c(x_axis,y_axis)]
    if(!is.null(sample)){
        s=sort(sample(1:nrow(xp),sample))
    } else {
        s=1:nrow(xp)
    }

    gate_updated=rep(0,nrow(xp))
    color_biplot_by_discrete(xp[s,],gate_updated[s],bty=bty,pch=pch,cex=cex,...)

    input.message=" n : new gate, c : redo last, s : stop gating. "
    cat("\nPlease use the mouse pointer to draw")
    polygons=list()
    i=0
    u="n"
    while(u!="s"){
        if(!u%in%c("n","s","c")){
            u=readline(paste("Incorrect input.",input.message,sep=""))
        }
        if(u=="n"){
            gate=gate_updated
            i=i+1
            col=setNames(c("black",rainbow(i)),0:i)

            new.pol=en.locator()

            new.pol=polygon.clean(new.pol)
            polygons=c(polygons,list(new.pol))
            gate_updated=update_gate(xp,polygons[[i]],gate,i)

            color_biplot_by_discrete(xp[s,],gate_updated[s],bty=bty,pch=pch,cex=cex,colors=col,...)

        }
        if(u=="c"){
            gate_updated=gate
            color_biplot_by_discrete(xp[s,],gate_updated[s],bty=bty,pch=pch,cex=cex,colors=col,...)
            new.pol=en.locator()
            new.pol=polygon.clean(new.pol)
            polygons[[i]]=new.pol
            gate_updated=update_gate(xp,polygons[[i]],gate_updated,i)
            color_biplot_by_discrete(xp[s,],gate_updated[s],bty=bty,pch=pch,cex=cex,colors=col,...)
        }
        u=readline(paste("Input ?",input.message,"\n",sep=""))
    }
    
    gate=gate_updated

    gate[apply(xp,1,function(x)any(is.na(x)))]=NA

    setNames(gate,rownames(matrix))
}

#' Colors a biplot according to a vector with discrete values
#' @param matrix a two columns matrix
#' @param discrete_vector a vector of size nrow(matrix)
#' @param colors Palette to used named after the unique elements of discrete_vector. Generated from rainbow() if missing.
#' @param pch passed to plot
#' @param cex passed to plot
#' @param bty passed to plot
#' @param ... passed to plot
#' @examples
#' data(Samusik_01_subset)
#' levels=unique(sort(Samusik_01_subset$labels))
#' colors=setNames(colorRampPalette(palette())(length(levels)),sort(levels))
#' with(Samusik_01_subset,color_biplot_by_discrete(matrix=tsne,discrete_vector=labels,colors=colors))
#' @export

color_biplot_by_discrete<-function(matrix,discrete_vector,...,bty="l",pch=16,cex=0.5,colors=NULL){
    levels=unique(discrete_vector)
    if(missing(colors)){
        colors=setNames(c("black",rainbow(length(levels)-1)),levels)
    }
    plot(matrix,bty=bty,pch=pch,cex=cex,col=colors[as.character(discrete_vector)],...)
}

#' Wrapper to locator that plots segments on the fly
en.locator<-function(){
    input=TRUE
    x=vector()
    y=vector()
    while(!is.null(input)){
        input=locator(1)
        x=c(x,input$x)
        y=c(y,input$y)

        if(length(x)>1){
            segments(x0=x[length(x)-1],x1=x[length(x)],y0=y[length(y)-1],y1=y[length(y)],lty=2)
        }
    }
    segments(x0=x[1],x1=x[length(x)],y0=y[1],y1=y[length(y)],lty=2)
    list(x=x,y=y)
}

#' Remove self intersection in polygons
#' @param poly a polygon (list with two components x and y which are equal-length numerical vectors)
#' @return A polygon without overlapping edges and new vertices corresponding to non-inner points of intersection

polygon.clean<-function(poly){
    if(requireNamespace("sf", quietly = TRUE)) {
        ## Create a data frame from the input coordinates
        coords_df = data.frame(x = poly$x, y = poly$y)

        ## Close polygon
        coords_df = rbind(coords_df, coords_df[1, ])
        
        ## Convert to an sf object
        coords_sf = sf::st_as_sf(coords_df, coords = c("x", "y"))
        
        ## Convert points to line
        line = sf::st_cast(sf::st_combine(coords_sf), "MULTILINESTRING")
        line = sf::st_node(line)
        
        ## Convert line to polygon
        polygon = sf::st_polygonize(line)

        ## Simplify and clean the polygon
        ## Note: st_union is used here for a cleaning effect similar to gUnaryUnion
        cleaned_polygon = sf::st_union(polygon, by_feature = TRUE)

        ## Extract coordinates
        coords = sf::st_coordinates(sf::st_cast(cleaned_polygon, "POLYGON"))
        x = coords[, 1]
        y = coords[, 2]

        ## Return the cleaned coordinates
        return(list(x = x, y = y))

    } else {
        return(poly)
    }
}

#' Updates a gate vector
#' @param xp A two colums matrix
#' @param polygon A list with two components x and y of equal lenghts and numeric values
#' @param gate_vector a vector of length nrow(xp) with integer values
#' @param value The number that will be assigned to gate_vector, corresponding to points that lie in the polygon
#' @return The updated gate_vector

update_gate=function(xp,polygon,gate_vector=rep(0,nrow(xp)),value=1){
    if(requireNamespace("sp",quietly=FALSE)){
        gate_vector[sp::point.in.polygon(xp[,1],xp[,2],polygon$x,polygon$y)!=0]=value
    }
    gate_vector
}

#' 2000 events randomly sampled from the 'Samusik_01' dataset
#' @name Samusik_01_subset
#' @docType data
#' @format list with four elements: fs_src (a flowSet), xp_src (its expression matrix), labels (manual gates of the events) and tsne (a tSNE projection of the dataset)
#' @references https://flowrepository.org/id/FR-FCM-ZZPH
"Samusik_01_subset"

#' @importFrom grDevices col2rgb
NULL
#' @importFrom grDevices dev.off
NULL
#' @importFrom grDevices png
NULL
#' @importFrom grDevices rainbow
NULL
#' @importFrom grDevices rgb
NULL
#' @importFrom graphics abline
NULL
#' @importFrom graphics locator
NULL
#' @importFrom graphics plot
NULL
#' @importFrom graphics segments
NULL
#' @importFrom graphics text
NULL
#' @importFrom graphics title
NULL
#' @importFrom stats setNames
NULL
#' @importFrom utils tail
NULL
