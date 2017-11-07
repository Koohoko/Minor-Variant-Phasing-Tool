# load the library
library(Rsamtools)
library(seqinr)
library(reshape2)
library(ggplot2)
options(stringsAsFactors = F)

setwd('c:/Users/Haogao/Google Drive/work/2017Sem3/2017-10-19_Pacbio_Phasing/')
readBAM <- function(bamFile){
    
    bam <- scanBam(bamFile)
    
    # A function for collapsing the list of lists into a single list
    # as per the Rsamtools vignette
    .unlist <- function (x){
        x1 <- x[[1L]]
        if (is.factor(x1)){
            structure(unlist(x), class = "factor", levels = levels(x1))
        } else {
            do.call(c, x)
        }
    }
    
    bam_field <- names(bam[[1]])
    
    list <- lapply(bam_field, function(y) .unlist(lapply(bam, "[[", y)))
    
    bam_df <- do.call("DataFrame", list)
    names(bam_df) <- bam_field
    
    #return a list that can be called as a data frame
    return(bam_df)
}

# Load the bam file
bam <- readBAM('bmut6.bam')
scanBamHeader('bmut6.bam')

ref = read.fasta('BR59mut_reference.fasta')
illu_result = read.csv('bmut6.csv')


for(i in 1:length(names(ref))){
    ref.name = names(ref)[i]
    ref.seq = ref[[i]]
    
    seq = bam$seq[bam$rname==ref.name]
    cigar = bam$cigar[bam$rname==ref.name]
    flag = bam$flag[bam$rname==ref.name]
    pos = bam$pos[bam$rname==ref.name]
    width = bam$qwidth[bam$rname==ref.name]
    
    #get the seqs by cigar
    seq.clip = mapply(function(x,y){
        cigar.1 = unlist(strsplit(x,"\\d+"))
        cigar.2 = unlist(strsplit(x,"\\D"))
        seq.clip.tmp = c()
        start = 1
        for(j in 1:length(cigar.2)){
            num = as.numeric(cigar.2[j])
            move = cigar.1[j+1]
            if (move %in% c('=','X')){
                seq.clip.tmp = paste0(seq.clip.tmp,toString(y[start:(start+num-1)]))
                start = start + num
            } else if(move == 'D'){
                seq.clip.tmp = paste0(seq.clip.tmp,c2s(rep('-',num)))
            } else {
                start = start +num
            }
        }
        return(seq.clip.tmp)
    },cigar,seq)
    
    #filter to get full-length seqs
    check = sapply(seq.clip,nchar)>length(ref.seq)*0.95
    seq.clip.filter = seq.clip[check]
    pos.filter = pos[check]

    df = data.frame()
    
    for(ii in 1:length(seq.clip.filter)){
        pos.tmp = pos.filter[ii] - 1
        df = rbind(df, c(rep('-',pos.tmp),s2c(seq.clip.filter[ii])))
    }
    
    ilu.r.tmp = illu_result[illu_result[,1] == ref.name,]
    df.tmp = df[,ilu.r.tmp[,2]]
    
    ##fill the gap with consensus
    df.tmp = as.data.frame(apply(df.tmp,2,function(x){x[x=='-']=names(which.max(table(x)));
    return(x)}))
    df.tmp = df.tmp[,apply(df.tmp,2,function(x){names(which.max(table(x)))})!='-']
    
    df.tmp.agg = aggregate(cbind(df.tmp[0],Count=1), by=df.tmp, length)
    df.tmp.agg = df.tmp.agg[order(df.tmp.agg[,ncol(df.tmp.agg)],decreasing = T),]
    
    names(df.tmp.agg)[1:nrow(ilu.r.tmp)] = paste0(ilu.r.tmp[,3],ilu.r.tmp[,2],ilu.r.tmp[,4])
    df.tmp.agg.t = t(df.tmp.agg)
    colnames(df.tmp.agg.t) = paste0('a',1:ncol(df.tmp.agg.t))
    
    df.tmp.agg.t.m = melt(df.tmp.agg.t[1:(nrow(df.tmp.agg.t)-1),])

    ggplot(df.tmp.agg.t.m,aes(x=Var2,y=Var1)) + geom_tile(aes(fill=value,
                            width=0.9, height=1)) +theme_classic() + 
        theme(axis.ticks = element_blank(),
              axis.line = element_blank())  + xlab('Count') +ylab('Position')+
        scale_x_discrete(labels = df.tmp.agg.t[nrow(df.tmp.agg.t),]) +
        ggtitle(ref.name) + 
        scale_fill_manual(name = "Nucleotide",
                          values = c('#5398D9','#F4E3B1','#D96B0C','#A53A3B'))+
        #'#F1F3F2',
        scale_y_discrete(limits = rev(levels(df.tmp.agg.t.m[,1])))
    ggsave(paste0(ref.name,'.jpg'),device = 'jpeg',width = 15,height = 15,dpi = 400)
}





