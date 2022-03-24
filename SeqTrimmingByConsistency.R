###Alignment cleaning tool to trim amd mask erroneous regions based on consistency
library(seqinr)


#setting parameters

#trim on both ends
windowSize=25
minGoodSites=18

#mask internally
windowSize_int=25
minGoodprop_int=0

#sequence-wide proportion of informative site threshold
prop_good=0

#trim end and mask internally until a windowsize contain at least 'minGoodSites' goodbases

TrimMaskSeq<-function(nSites,myseq,y,outputfile){  #define function, process one seq each time
    
    nGood=0
    for (site in 1:nSites){
        if (myseq[site]=='a' | myseq[site]=='t' | myseq[site]=='c' | myseq[site]=='g'){nGood=nGood+1}
    }
    #if (nGood>0.1*nSites){
    if (nGood>0){
    #trim left end
    firstGoodBase=nSites

    for (site in 1:(nSites-windowSize-1)){
        nGood=0
        lastBadPos=-1   #potentially need revise
        for (l in 1:windowSize){
            #if (is.na(y[site+l])){nGood=nGood+1}
            if(myseq[site+l]==y[site+l]){nGood=nGood+1}
            else {lastBadPos=l}
        }
        if (nGood>=minGoodSites){
            #firstGoodBase=site+lastBadPos+1
            firstGoodBase=site
            break
        }
    }
    #trim right end
    lastGoodBase=-1
    
    for (site in nSites:(windowSize+1)){
        nGood=0
        lastBadPos=-1
        for (l in 1:windowSize){
            #if (is.na(y[site-l])){nGood=nGood+1}
            if(myseq[site-l]==y[site-l]){nGood=nGood+1}
            else {lastBadPos=l}
        }
        if (nGood>=minGoodSites){
            #lastGoodBase=site-(lastBadPos+1)
            lastGoodBase=site
            break
        }
    }
    #mask internally
    if (firstGoodBase<nSites-windowSize-1 & lastGoodBase>firstGoodBase){
    w=list()
    for (site in (firstGoodBase+windowSize):(lastGoodBase-windowSize_int-1)){
        nGood=0
        gap=0
        #bad=list()
        #containedGap=FALSE
        for (l in 1:windowSize_int){
            #if (is.na(y[site+i])){nGood=nGood+1}
            if (myseq[site+l]==y[site+l]){
            	nGood=nGood+1
            }
            else if (myseq[site+l]=='-'){gap=gap+1}
        }
        if (windowSize_int!=gap){
        	if (nGood/(windowSize_int-gap)<minGoodprop_int){
            	for (i in 1:windowSize_int){w[site+l]=1}
        	}
        }
    }
    
    #clean left ends of bad internal regions (XXXXATGGGC)
    bad_int=list()
    bad_int=which(w %in% c(1))
    if (length(bad_int>0)){
    	#print(attributes(myseq)$name)
    	for (site in bad_int){
    		back_pos=0
    		repeat{
    			w[site-back_pos]=1
  				#examine the sites before, expand 'bad' window until there are three consecutive good site
    			s1=is.null(w[[site-back_pos-1]]) & myseq[site-back_pos-1]==y[site-back_pos-1]
    			s2=is.null(w[[site-back_pos-2]]) & myseq[site-back_pos-2]==y[site-back_pos-2]
    			s3=is.null(w[[site-back_pos-3]]) & myseq[site-back_pos-3]==y[site-back_pos-3]
  				back_pos=back_pos+1
  				if(s1+s2+s3==3 | site-back_pos==0){
    				break
 				}
			}
    	}
    }
    
    
    #adjust the sequence accordingly
    for (site in 1:nSites){
        if (site<firstGoodBase | site>lastGoodBase){myseq[site]='-'} #trim
        else if (site<=length(w) && !is.null(w[[site]])){myseq[site]='-'} #mask
    }
    
    #filter taxa by proportion of good sites
    nGood=0
    for (site in 1:nSites){
        if (myseq[site]=='a' | myseq[site]=='t' | myseq[site]=='c' | myseq[site]=='g'){nGood=nGood+1}
    }
    if (nGood>prop_good*nSites){	
    	write.fasta(myseq,attr(myseq,'name'),file.out=outputfile,open='a')
    }
    }
    }
}



locilist=list.files(pattern='*.nolongbr.aln.fas')
jdata=read.table('jansen_lab_data.txt')

for (i in 1:length(locilist)){ 
	x=read.fasta(locilist[i],seqtype = 'DNA')
	#get consensus seq based on majority rule
	con=consensus(read.alignment(locilist[i],format='fasta'),method='majority',type='DNA')
	alignment_length=length(x[[1]])
	for (j in 1:length(names(x))){
		if (attr(x[[j]],'name') %in% jdata$V1){
			#only trim jansen data
			TrimMaskSeq(alignment_length,x[[j]],con, paste(locilist[i],'.Rtrim.fas',sep=''))
		}else{
			#directly write zhao's data
			write.fasta(x[[j]],attr(x[[j]],'name'),file.out=paste(locilist[i],'.Rtrim.fas',sep=''),open='a')
		}
	}
	#iterate among taxa and write the modified alignemnt to new file
}
