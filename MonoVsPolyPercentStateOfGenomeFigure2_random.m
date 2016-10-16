clear
directory='C:\Data\ChipSeqPaper\Data\RandomlyDownsampled\BamsIntersectedWithEncode\'
groupNames={'Mono'; 'Poly'}
BTLblue=[15 121 191]/255 %0    0.4706    0.6706  is a little nicer

ps=0;

annTypes ={
    'CTCF'
    'PF'
    'TSS'
    'T'
     'E' 
     'WE'
     
    'R'
}

colorscheme=[ 
 0.8549    0.6784    0.5255 %   'CTCF' 
  1.0000    0.4000    0.6980;%PF'
 0.0588    0.4745    0.7490; %   TSS 
 0.5294    0.8078    0.9804; %   'T'
0.2431    0.6863    0.4627; %   'E'
 0.5882    0.9765    0.4824; %   'WE'
  

 1.0000    .3   .3;%   'R'
]
    





sampleNamesVals={};
ctrVals=1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

sampleNames={


'WCE_K562.bam.txt'
};


sampleNamesShort={
'WCE'
};


group=[1 1 2 2];


clear data
clear thisdata
for i=1:length(sampleNames)
    
    fileName=strcat(directory, sampleNames{i});
    [ chr, startPosK27, endPosK27, overlappingReads, overlappingBases, featureLength, annotationType ] = ReadBedtoolsCoverage( fileName );
    data(:,i)=overlappingReads;
    
end



for i=1:size(data,2)
    total=sum(data(:,i));
    for j=1:length(annTypes)
        tableVals(ctrVals,j)=sum(data(strcmp(annotationType, annTypes{j}), i)); 
    end
    sampleNamesVals{ctrVals}=sampleNames{i};
    tablePerc(ctrVals,:)=tableVals(ctrVals,:)/total;    
  thisdata(i,:)=tableVals(ctrVals,:) 
    
    yTickLabels{ctrVals}=sampleNamesShort{i};
    ctrVals=ctrVals+1;    
end

%[p,tbl] = anova2(thisdata,3);
%ps(length(ps)+1)=p(2);

tablePerc(ctrVals,:)=zeros(1, size(tablePerc,2));
yTickLabels{ctrVals}='';
ctrVals=ctrVals+1;  


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

sampleNames={

'mAb-K27ac-1.dupsRemoved.bam.bam.13407899.bam.txt'
'mAb-K27ac-2.dupsRemoved.bam.bam.13407899.bam.txt'
'pAb-K27ac-1.dupsRemoved.bam.bam.13407899.bam.txt'
'pAb-K27ac-2.dupsRemoved.bam.bam.13407899.bam.txt'
};


sampleNamesShort={
'H3K27ac Mono Rep 1'
'H3K27ac Mono Rep 2'
'H3K27ac Poly Rep 1'
'H3K27ac Poly Rep 2'
};


group=[1 1 2 2];


clear data
clear thisdata
for i=1:length(sampleNames)
    
    fileName=strcat(directory, sampleNames{i});
    [ chr, startPosK27, endPosK27, overlappingReads, overlappingBases, featureLength, annotationType ] = ReadBedtoolsCoverage( fileName );
    data(:,i)=overlappingReads;
    
    
end

for j=1:length(annTypes)
      good=strcmp(annotationType, annTypes{j});
      type(j)=sum(endPosK27(good)-startPosK27(good)+1);
end
    


for i=1:size(data,2)
    total=sum(data(:,i));
    for j=1:length(annTypes)
        tableVals(ctrVals,j)=sum(data(strcmp(annotationType, annTypes{j}), i)); 
    end
    sampleNamesVals{ctrVals}=sampleNames{i};
    tablePerc(ctrVals,:)=tableVals(ctrVals,:)/total; 
  thisdata(i,:)=tableVals(ctrVals,:) 
    
    yTickLabels{ctrVals}=sampleNamesShort{i};
    ctrVals=ctrVals+1;    
end


% [p,tbl] = anova2(thisdata,2);
%ps(length(ps)+1)=p(2);


tablePerc(ctrVals,:)=zeros(1, size(tablePerc,2))
yTickLabels{ctrVals}='';
ctrVals=ctrVals+1;  




mono=sum(data(:, group==1),2);
poly=sum(data(:, group==2),2);
figure
loglog(mono(strcmp(annotationType, 'TSS')), poly(strcmp(annotationType, 'TSS')), '.', 'Color', [ 0.0588    0.2045    0.8090])
hold on
loglog(mono(strcmp(annotationType, 'E')), poly(strcmp(annotationType, 'E')), '.', 'Color', [0.1231    0.6863    0.2627] )
legend('TSS', 'Enhancers', 'Location', 'NorthWest')
xlabel('Monoclonal H3K27ac')
ylabel('Polyclonal H3K27ac')
plotBeef(12)


%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

sampleNames={


'1_Exp4-K562-2A-K27me3-mono_H9V71ADXX.1.aligned.duplicates_marked.dupsRemoved.bam.bam.9356733.bam.txt'
'1_Exp4-K562-2B-K27me3-mono_H9V71ADXX.1.aligned.duplicates_marked.dupsRemoved.bam.bam.9356733.bam.txt'
'1_Exp4-K562-2C-K27me3-mono_H9V71ADXX.1.aligned.duplicates_marked.dupsRemoved.bam.bam.9356733.bam.txt'
'1_Exp4-K562-2D-K27me3-mono_H9V71ADXX.1.aligned.duplicates_marked.dupsRemoved.bam.bam.9356733.bam.txt'
'1_Exp4-K562-2E-K27me3-poly_H9V71ADXX.1.aligned.duplicates_marked.dupsRemoved.bam.bam.9356733.bam.txt'
'1_Exp4-K562-2F-K27me3-poly_H9V71ADXX.1.aligned.duplicates_marked.dupsRemoved.bam.bam.9356733.bam.txt'
'1_Exp4-K562-2G-K27me3-poly_H9V71ADXX.1.aligned.duplicates_marked.dupsRemoved.bam.bam.9356733.bam.txt'
};


sampleNamesShort={
'H3K27me3 Mono Rep 1'
'H3K27me3 Mono Rep 2'
'H3K27me3 Mono Rep 3'
'H3K27me3 Mono Rep 4'
'H3K27me3 Poly Rep 1'
'H3K27me3 Poly Rep 2'
'H3K27me3 Poly Rep 3'
};


group=[1 1 2 2];


clear data
clear thisdata
for i=1
    
    fileName=strcat(directory, sampleNames{i});
    [ chr, startPosK27, endPosK27, overlappingReads, overlappingBases, featureLength, annotationType ] = ReadBedtoolsCoverage( fileName );
    data(:,i)=overlappingReads;
    
end



for i=1:size(data,2)
    total=sum(data(:,i));
    for j=1:length(annTypes)
        tableVals(ctrVals,j)=sum(data(strcmp(annotationType, annTypes{j}), i)); 
    end
    sampleNamesVals{ctrVals}=sampleNames{i};
    tablePerc(ctrVals,:)=tableVals(ctrVals,:)/total;    
  thisdata(i,:)=tableVals(ctrVals,:) 
    
    yTickLabels{ctrVals}=sampleNamesShort{i};
    ctrVals=ctrVals+1;    
end

%[p,tbl] = anova2(thisdata,3);
%ps(length(ps)+1)=p(2);

tablePerc(ctrVals,:)=zeros(1, size(tablePerc,2));
yTickLabels{ctrVals}='';
ctrVals=ctrVals+1;  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
sampleNames={

'mAb-K36me3-1.dupsRemoved.bam.bam.16088020.bam.txt'
'mAb-K36me3-2.dupsRemoved.bam.bam.16088020.bam.txt'
'pAb-K36me3-1.dupsRemoved.bam.bam.16088020.bam.txt'
'pAb-K36me3-2.dupsRemoved.bam.bam.16088020.bam.txt'
};


sampleNamesShort={
'H3K36me3 Mono Rep 1'
'H3K36me3 Mono Rep 2'
'H3K36me3 Poly Rep 1'
'H3K36me3 Poly Rep 2'
}

group=[1 1 2 2];


clear data
clear thisdata
for i=1:length(sampleNames)
    
    fileName=strcat(directory, sampleNames{i});
    [ chr, startPosK27, endPosK27, overlappingReads, overlappingBases, featureLength, annotationType ] = ReadBedtoolsCoverage( fileName );
    data(:,i)=overlappingReads;
    
end



for i=1:size(data,2)
    total=sum(data(:,i));
    for j=1:length(annTypes)
        tableVals(ctrVals,j)=sum(data(strcmp(annotationType, annTypes{j}), i)); 
    end
    sampleNamesVals{ctrVals}=sampleNames{i};
    tablePerc(ctrVals,:)=tableVals(ctrVals,:)/total;   
  thisdata(i,:)=tableVals(ctrVals,:) 
    thisdataVals(i,:)=tableVals(ctrVals,:);  
    
    yTickLabels{ctrVals}=sampleNamesShort{i};
    ctrVals=ctrVals+1;    
end

% [p,tbl] = anova2(thisdata,2);
%ps(length(ps)+1)=p(2);

tablePerc(ctrVals,:)=zeros(1, size(tablePerc,2))
ctrVals=ctrVals+1;  


%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

sampleNames={

'K562_K4me1_Ido.dupsRemoved.bam.bam.20565471.bam.txt'
'K562_K4me1_LizExp11.dupsRemoved.bam.bam.20565471.bam.txt'
'mAb-K4me1-1.dupsRemoved.bam.bam.20565471.bam.txt'
'pAb-K4me1-1.dupsRemoved.bam.bam.20565471.bam.txt'
'pAb-K4me1-2.dupsRemoved.bam.bam.20565471.bam.txt'

};

sampleNamesShort={
'H3K4me1 Mono Rep 1'
'H3K4me1 Mono Rep 2'
'H3K4me1 Mono Rep 3'
'H3K4me1 Poly Rep 1'
'H3K4me1 Poly Rep 2'
}


group=[1 1 1 2 2];


clear data
clear thisdata
for i=1:length(sampleNames)
    
    fileName=strcat(directory, sampleNames{i});
    [ chr, startPosK27, endPosK27, overlappingReads, overlappingBases, featureLength, annotationType ] = ReadBedtoolsCoverage( fileName );
    data(:,i)=overlappingReads;
    
end



for i=1:size(data,2)
    total=sum(data(:,i));
    for j=1:length(annTypes)
        tableVals(ctrVals,j)=sum(data(strcmp(annotationType, annTypes{j}), i)); 
    end
    sampleNamesVals{ctrVals}=sampleNames{i};
    tablePerc(ctrVals,:)=tableVals(ctrVals,:)/total;   
  thisdata(i,:)=tableVals(ctrVals,:) 
    yTickLabels{ctrVals}=sampleNamesShort{i};
    ctrVals=ctrVals+1;    
end

% [p,tbl] = anova2(thisdata,2);
%ps(length(ps)+1)=p(2);

tablePerc(ctrVals,:)=zeros(1, size(tablePerc,2));
yTickLabels{ctrVals}='';
ctrVals=ctrVals+1;  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

sampleNames={

'1_Exp4-K562-1A-K4me3-mono_H9V71ADXX.1.aligned.duplicates_marked.dupsRemoved.bam.bam.4955064.bam.txt'
'1_Exp4-K562-1B-K4me3-mono_H9V71ADXX.1.aligned.duplicates_marked.dupsRemoved.bam.bam.4955064.bam.txt'
'1_Exp4-K562-1C-K4me3-mono_H9V71ADXX.1.aligned.duplicates_marked.dupsRemoved.bam.bam.4955064.bam.txt'
'1_Exp4-K562-1D-K4me3-mono_H9V71ADXX.1.aligned.duplicates_marked.dupsRemoved.bam.bam.4955064.bam.txt'
'1_Exp4-K562-1E-K4me3-poly_H9V71ADXX.1.aligned.duplicates_marked.dupsRemoved.bam.bam.4955064.bam.txt'
'1_Exp4-K562-1F-K4me3-poly_H9V71ADXX.1.aligned.duplicates_marked.dupsRemoved.bam.bam.4955064.bam.txt'
'1_Exp4-K562-1G-K4me3-poly_H9V71ADXX.1.aligned.duplicates_marked.dupsRemoved.bam.bam.4955064.bam.txt'
'1_Exp4-K562-1H-K4me3-poly_H9V71ADXX.1.aligned.duplicates_marked.dupsRemoved.bam.bam.4955064.bam.txt'

};


sampleNamesShort={
'H3K4me3 Mono Rep 1'
'H3K4me3 Mono Rep 2'
'H3K4me3 Mono Rep 3'
'H3K4me3 Mono Rep 4'
'H3K4me3 Poly Rep 1'
'H3K4me3 Poly Rep 2'
'H3K4me3 Poly Rep 3'
'H3K4me3 Poly Rep 4'
};

group=[1 1 2 2];


clear data
clear thisdata
for i=1:length(sampleNames)
    
    fileName=strcat(directory, sampleNames{i});
    [ chr, startPosK27, endPosK27, overlappingReads, overlappingBases, featureLength, annotationType ] = ReadBedtoolsCoverage( fileName );
    data(:,i)=overlappingReads;
    
end



for i=1:size(data,2)
    total=sum(data(:,i));
    for j=1:length(annTypes)
        tableVals(ctrVals,j)=sum(data(strcmp(annotationType, annTypes{j}), i)); 
    end
    sampleNamesVals{ctrVals}=sampleNames{i};
    tablePerc(ctrVals,:)=tableVals(ctrVals,:)/total;
  thisdata(i,:)=tableVals(ctrVals,:) 
    yTickLabels{ctrVals}=sampleNamesShort{i};
    ctrVals=ctrVals+1;    
end





%[p,tbl] = anova2(thisdata,4);
%ps(length(ps)+1)=p(2);

tablePerc(ctrVals,:)=zeros(1, size(tablePerc,2))
yTickLabels{ctrVals}='';
ctrVals=ctrVals+1;  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

sampleNames={
'K562_K9me3_Ido.dupsRemoved.bam.bam.18740026.bam.txt'
'K562_K9me3_LizExp11.dupsRemoved.bam.bam.18740026.bam.txt'
'mAb-K9me3-2.dupsRemoved.bam.bam.18740026.bam.txt'
'pAb-K9me3-1.dupsRemoved.bam.bam.18740026.bam.txt'
'pAb-K9me3-2.dupsRemoved.bam.bam.18740026.bam.txt'

};

 


sampleNamesShort={
'H3K9me3 Mono Rep 1'
'H3K9me3 Mono Rep 2'
'H3K9me3 Mono Rep 3'
'H3K9me3 Poly Rep 1'
'H3K9me3 Poly Rep 2'
};

group=[1 1 1 2 2];


clear data
clear thisdata
clear thisdataVals
for i=1:length(sampleNames)
    
    fileName=strcat(directory, sampleNames{i});
    [ chr, startPosK27, endPosK27, overlappingReads, overlappingBases, featureLength, annotationType ] = ReadBedtoolsCoverage( fileName );
    data(:,i)=overlappingReads;
    
end



for i=1:size(data,2)
    total=sum(data(:,i));
    for j=1:length(annTypes)
        tableVals(ctrVals,j)=sum(data(strcmp(annotationType, annTypes{j}), i)); 
    end
    sampleNamesVals{ctrVals}=sampleNames{i};
    tablePerc(ctrVals,:)=tableVals(ctrVals,:)/total;    
  thisdata(i,:)=tableVals(ctrVals,:) 
    
    yTickLabels{ctrVals}=sampleNamesShort{i};
    ctrVals=ctrVals+1;    
end

% [p,tbl] = anova2(thisdata,2);
%ps(length(ps)+1)=p(2);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all
fig=figure();
barh(tablePerc*100, 'stacked')
colormap(colorscheme)

set(gca, 'YTick', 1:ctrVals-1)
set(gca, 'YTickLabel', yTickLabels)
set(gca, 'TickLength', [0 0]);
%set(gca,'View',[0 -90])
axis([0 100.5 -.025  size(tablePerc,1)+1])

legend(annTypes, 'location', 'SouthOutside','Orientation','horizontal')

title('% Reads Mapping to Cannonical Genome Regions')
set(gca,'View',[0 -90])
outFileName='C:\Data\ChipSeqPaper\Figures\ReadsByEncodeTrackRandom.tif'
 fig = gcf;
    fig.PaperUnits = 'inches';
    fig.PaperPosition = [0.6 0.6 6 8];
    fig.PaperPositionMode = 'manual';
    print(outFileName,'-dtiff','-r350')

    

outFileName='C:\Data\ChipSeqPaper\Figures\ReadsByEncodeTrackRandom.pdf'
 fig = gcf;
    fig.PaperUnits = 'inches';
    fig.PaperPosition = [0.6 0.6 6 8];
    fig.PaperPositionMode = 'manual';
    print(outFileName,'-dpdf','-r350')




