#Copyright (c) 2017 Haogao Gu. All rights reserved.
library(shiny)
library(Rsamtools)
library(seqinr)
library(reshape2)
library(ggplot2)
options(stringsAsFactors = F)
options(shiny.maxRequestSize=30*1024^2)
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

shinyServer(function(input, output, session) {
    data <- reactive({
        inFile <- input$reffile   
        if (is.null(inFile))
            return(NULL)
        return(read.fasta(inFile$datapath))
    })
    bamfile.r <- reactive({
        file1 = input$bamfile   
        if (is.null(file1))
            return(NULL)
        return(readBAM(file1$datapath))
    })
    ill.r <- reactive({
        file2 = input$snpfile   
        if (is.null(file2))
            return(NULL)
        return(read.csv(file2$datapath))
    })
    refname = reactive({input$var})
    
    
    output$selector <- renderUI({
    selectInput("var", "Choose segment:", as.list(names(data()))) 
    })
    
    # ifsnp <- reactive({
    #     if(!input$usesnp) return(NULL)
    #     return(input$usesnp)
    # })
    # 
    # observe({
    #     updateCheckboxInput(session, "usesnp",value = input$usesnp)
    # })
    

    # observeEvent(input$insertBtn, {
    #         insertUI(
    #             selector = '#placeholder',
    #             #ui = textInput('txt8', 'Insert some text'),
    #             immediate = TRUE,
    #             ui = fileInput('snpfile', 'Choose SNPs File')
    #         )
    # },once=TRUE)
    
    
    output$myImage <- renderImage({
        validate(
            need(refname(), "Sorry, there is no data for your requested analysis."
            )
        )
        ################################################################################        
        ref = data()
        bam = bamfile.r()
        illu_result = ill.r()
        ref.seq = ref[[refname()]]

        seq = bam$seq[bam$rname==refname()]
        cigar = bam$cigar[bam$rname==refname()]
        flag = bam$flag[bam$rname==refname()]
        pos = bam$pos[bam$rname==refname()]
        width = bam$qwidth[bam$rname==refname()]
        
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

        for(jj in 1:length(seq.clip.filter)){
            pos.tmp = pos.filter[jj] - 1
            df = rbind(df, c(rep('-',pos.tmp),s2c(seq.clip.filter[jj])))
        }

        #illumina snps filtering
        if(input$usesnp){
            ilu.r.tmp = illu_result[illu_result[,1] == refname(),]
            df.tmp = df[,ilu.r.tmp[,2]]
        } else {
            df.tmp = df
        }


        ##fill the gap with consensus
        df.tmp = as.data.frame(apply(df.tmp,2,function(x){x[x=='-']=names(which.max(table(x)));
        return(x)}))
        names.tmp = paste0(ilu.r.tmp[,3],ilu.r.tmp[,2],
                           ilu.r.tmp[,4])[apply(df.tmp,2,function(x){names(which.max(table(x)))})!='-']
        df.tmp = df.tmp[,apply(df.tmp,2,function(x){names(which.max(table(x)))})!='-']

        df.tmp.agg = aggregate(cbind(df.tmp[0],Count=1), by=df.tmp, length)
        df.tmp.agg = df.tmp.agg[order(df.tmp.agg[,ncol(df.tmp.agg)],decreasing = T),]
        
        if(input$usesnp){
            names(df.tmp.agg)[1:(ncol(df.tmp.agg)-1)] = names.tmp
        } else {
            df.tmp.agg = df.tmp.agg[,apply(df.tmp.agg,2,function(x){length(table(x))>1})]
            names(df.tmp.agg)[1:(ncol(df.tmp.agg)-1)] = gsub('\\D',
                                                             '',names(df.tmp.agg)[1:(ncol(df.tmp.agg)-1)])
        }

        df.tmp.agg.t = t(df.tmp.agg)
        colnames(df.tmp.agg.t) = paste0('a',1:ncol(df.tmp.agg.t))
        df.tmp.agg.t.m = melt(df.tmp.agg.t[1:(nrow(df.tmp.agg.t)-1),])
        
        ##ploting
        ggplot(df.tmp.agg.t.m,aes(x=Var2,y=factor(Var1))) + geom_tile(aes(fill=value,
                                                                  width=0.9, height=1)) +theme_classic() +
            theme(axis.ticks = element_blank(),
                  axis.line = element_blank())  + xlab('Count') +ylab('Position')+
            scale_x_discrete(labels = df.tmp.agg.t[nrow(df.tmp.agg.t),]) +
            ggtitle(refname()) +
            scale_fill_manual(name = "Nucleotide",
                              values = c('#5398D9','#F4E3B1','#D96B0C','#A53A3B'))+
            #'#F1F3F2',
            scale_y_discrete(limits = rev(levels(df.tmp.agg.t.m[,1])))
        if(input$usesnp){
            ggsave('./images/image1.jpg',device = 'jpeg',width = 15,height = 15,dpi = 300)
        } else {
            ggsave('./images/image1.jpg',device = 'jpeg',width = 60,height = 60,dpi = 150,limitsize = F)
        }

        ################################################################################    

        if(length(dir('images/'))>1){
            filename <-'./images/image1.jpg'
        } else {
            filename <-'./images/tset.png'
        }
        
        # Return a list containing the filename and alt text
        list(src = filename,
             width = '800',
             height = '800',
             alt = 'Result')
        
    }, deleteFile = FALSE)
})