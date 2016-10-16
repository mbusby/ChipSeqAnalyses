clear
directory='C:\Data\ChipSeqPaper\Data\DownsampledByInsertSize\BamsIntersectedWithEncode\'
outFileName='C:\Data\ChipSeqPaper\Figures\ReadsByEncodeTrack_insertSize.pdf'
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


group=[1];


clear data
clear thisdata
for i=1:length(sampleNames)
    
    fileName=strcat(directory, sampleNames{i});
    [ ~, ~, ~, overlappingReads, overlappingBases, featureLength, annotationType ] = ReadBedtoolsCoverage( fileName );
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


tablePerc(ctrVals,:)=zeros(1, size(tablePerc,2))
yTickLabels{ctrVals}='';
ctrVals=ctrVals+1;  




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

sampleNames={

'mAb-K27ac-1.dupsRemoved.bam.txt'
'mAb-K27ac-2.dupsRemoved.bam.txt'
'pAb-K27ac-1.dupsRemoved.bam.txt'
'pAb-K27ac-2.dupsRemoved.bam.txt'

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
    [ ~, ~, ~, overlappingReads, overlappingBases, featureLength, annotationType ] = ReadBedtoolsCoverage( fileName );
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


tablePerc(ctrVals,:)=zeros(1, size(tablePerc,2))
yTickLabels{ctrVals}='';
ctrVals=ctrVals+1;  


mono=sum(data(:, group==1),2);
poly=sum(data(:, group==2),2);

fig=figure()
subplot(1,3,1)

TSSMono=mono(strcmp(annotationType, 'TSS'));
TSSPoly=poly(strcmp(annotationType, 'TSS'));
r=corr(TSSMono, TSSPoly)^2;
rT=sprintf('%.2f',r)
rText=strcat('R^2=',rT);

loglog(mono(strcmp(annotationType, 'TSS')), poly(strcmp(annotationType, 'TSS')), '.', 'Color', [ 0.0588    0.2045    0.8090])
hold on
loglog(mono(strcmp(annotationType, 'TSS')), mono(strcmp(annotationType, 'TSS')), 'Color', 'r')
h=legend('TSS','Location', 'NorthWest')
text(70000,2, rText, 'HorizontalAlignment', 'right')

set(h,'FontSize',9);
xlabel('Monoclonal H3K27ac')
ylabel('Polyclonal H3K27ac')
plotBeef(10)

subplot(1,3,2)

loglog(mono, poly, '.', 'Color', [ 0.7 0.7 0.7])
hold on
loglog(mono(strcmp(annotationType, 'TSS')), poly(strcmp(annotationType, 'TSS')), '.', 'Color', [ 0.0588    0.2045    0.8090])

loglog(mono(strcmp(annotationType, 'E')), poly(strcmp(annotationType, 'E')), '.', 'Color', [0.1231    0.6863    0.2627] )
h=legend('All Points', 'TSS', 'Enhancers', 'Location', 'NorthWest')
r=corr(mono, poly)^2;
rT=sprintf('%.2f',r)
rText=strcat('R^2 (All)=',rT);

text(70000,2, rText, 'HorizontalAlignment', 'right')

set(h,'FontSize',9);
xlabel('Monoclonal H3K27ac')
ylabel('Polyclonal H3K27ac')
plotBeef(10)

subplot(1,3,3)


r=corr(mono(strcmp(annotationType, 'E')), poly(strcmp(annotationType, 'E')))^2;
rT=sprintf('%.2f',r)
rText=strcat('R^2=',rT);

loglog(mono(strcmp(annotationType, 'E')), poly(strcmp(annotationType, 'E')), '.', 'Color', [0.1231    0.6863    0.2627])
hold on
loglog(mono(strcmp(annotationType, 'E')), mono(strcmp(annotationType, 'E')),  'Color', 'r')
h=legend('Enhancers','Location', 'NorthWest')
set(h,'FontSize',9);
xlabel('Monoclonal H3K27ac')
ylabel('Polyclonal H3K27ac')
text(70000,2, rText, 'HorizontalAlignment', 'right')
plotBeef(10)
axis([0,  10^5, 0,  10^5])


    fig.PaperUnits = 'inches';
    fig.PaperPosition = [0.6 0.6 10.5 2.9];
    fig.PaperPositionMode = 'manual';
    fig.PaperOrientation = 'landscape';
    %print('C:\Data\ChipSeqPaper\Figures\scatterH3K27acByTSSvsEnhancer_InsertSample','-dpdf','-r350')
    
    fig.PaperUnits = 'inches';
    fig.PaperPosition = [0.6 0.6 10.5 2.9];
    fig.PaperPositionMode = 'manual';
    fig.PaperOrientation = 'landscape';
    %print('C:\Data\ChipSeqPaper\Figures\scatterH3K27acByTSSvsEnhancer_InsertSample.tif','-dtiff','-r350')
    
    
%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

sampleNames={


'1_Exp4-K562-2A-K27me3-mono_H9V71ADXX.1.aligned.duplicates_marked.dupsRemoved.bam.txt'
'1_Exp4-K562-2B-K27me3-mono_H9V71ADXX.1.aligned.duplicates_marked.dupsRemoved.bam.txt'
'1_Exp4-K562-2C-K27me3-mono_H9V71ADXX.1.aligned.duplicates_marked.dupsRemoved.bam.txt'
'1_Exp4-K562-2D-K27me3-mono_H9V71ADXX.1.aligned.duplicates_marked.dupsRemoved.bam.txt'
'1_Exp4-K562-2E-K27me3-poly_H9V71ADXX.1.aligned.duplicates_marked.dupsRemoved.bam.txt'
'1_Exp4-K562-2F-K27me3-poly_H9V71ADXX.1.aligned.duplicates_marked.dupsRemoved.bam.txt'
'1_Exp4-K562-2G-K27me3-poly_H9V71ADXX.1.aligned.duplicates_marked.dupsRemoved.bam.txt'

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
for i=1:length(sampleNames)
    
    fileName=strcat(directory, sampleNames{i});
    [ ~, ~, ~, overlappingReads, overlappingBases, featureLength, annotationType ] = ReadBedtoolsCoverage( fileName );
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

'mAb-K36me3-1.dupsRemoved.bam.16088020.bam.txt'
'mAb-K36me3-2.dupsRemoved.bam.16088020.bam.txt'
'pAb-K36me3-1.dupsRemoved.bam.16088020.bam.txt'
'pAb-K36me3-2.dupsRemoved.bam.16088020.bam.txt'
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
    [ ~, ~, ~, overlappingReads, overlappingBases, featureLength, annotationType ] = ReadBedtoolsCoverage( fileName );
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

'K562_K4me1_Ido.dupsRemoved.bam.txt' 
'K562_K4me1_LizExp11.dupsRemoved.bam.txt' 
'mAb-K4me1-1.dupsRemoved.bam.txt'
'pAb-K4me1-1.dupsRemoved.bam.txt'
'pAb-K4me1-2.dupsRemoved.bam.txt'

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
    [ ~, ~, ~, overlappingReads, overlappingBases, featureLength, annotationType ] = ReadBedtoolsCoverage( fileName );
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

'1_Exp4-K562-1A-K4me3-mono_H9V71ADXX.1.aligned.duplicates_marked.dupsRemoved.bam.txt'
'1_Exp4-K562-1B-K4me3-mono_H9V71ADXX.1.aligned.duplicates_marked.dupsRemoved.bam.txt'
'1_Exp4-K562-1C-K4me3-mono_H9V71ADXX.1.aligned.duplicates_marked.dupsRemoved.bam.txt'
'1_Exp4-K562-1D-K4me3-mono_H9V71ADXX.1.aligned.duplicates_marked.dupsRemoved.bam.txt'
'1_Exp4-K562-1E-K4me3-poly_H9V71ADXX.1.aligned.duplicates_marked.dupsRemoved.bam.txt'
'1_Exp4-K562-1F-K4me3-poly_H9V71ADXX.1.aligned.duplicates_marked.dupsRemoved.bam.txt'
'1_Exp4-K562-1G-K4me3-poly_H9V71ADXX.1.aligned.duplicates_marked.dupsRemoved.bam.txt'
'1_Exp4-K562-1H-K4me3-poly_H9V71ADXX.1.aligned.duplicates_marked.dupsRemoved.bam.txt'

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
    [ ~, ~, ~, overlappingReads, overlappingBases, featureLength, annotationType ] = ReadBedtoolsCoverage( fileName );
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
'K562_K9me3_Ido.dupsRemoved.bam.txt'   
'K562_K9me3_LizExp11.dupsRemoved.bam.txt '
'mAb-K9me3-2.dupsRemoved.bam.txt'
'pAb-K9me3-1.dupsRemoved.bam.txt'
'pAb-K9me3-2.dupsRemoved.bam.txt'

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
    [ ~, ~, ~, overlappingReads, overlappingBases, featureLength, annotationType ] = ReadBedtoolsCoverage( fileName );
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



outFileName='C:\Data\ChipSeqPaper\Figures\ReadsByEncodeTrack_insertSize.tif'

set(gca,'View',[0 -90])
 fig = gcf;
    fig.PaperUnits = 'inches';
    fig.PaperPosition = [0.6 0.6 6 8];
    fig.PaperPositionMode = 'manual';
    print(outFileName,'-dtiff','-r350')

    

outFileName='C:\Data\ChipSeqPaper\Figures\ReadsByEncodeTrack_insertSize.pdf'


 fig = gcf;
    fig.PaperUnits = 'inches';
    fig.PaperPosition = [0.6 0.6 6 8];
    fig.PaperPositionMode = 'manual';
 %   print(outFileName,'-dpdf','-r350')


    
    
    
    
    
    
    


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%HELA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

directory='C:\Data\ChipSeqPaper\Data\IntersectedWithCombinedMarksHeLa\'
groupNames={'Mono'; 'Poly'}


sampleNames={
'1_HeLa_synchronized_K27ac_-_active_motif_mAb_HAL74ADXX.1.aligned.duplicates_marked.dupsRemoved.bam.bam.13500000.bam.txt'
'1_HeLa_synchronized_K27ac_-_active_motif_pAb_HAL74ADXX.1.aligned.duplicates_marked.dupsRemoved.bam.bam.13500000.bam.txt'
'1_HeLa_synchronized_K27ac_-_cell_signaling_HAL74ADXX.1.aligned.duplicates_marked.dupsRemoved.bam.bam.13500000.bam.txt'
};

group=[1 2 3];
groupNames={'Active Mono'; 'Active Poly'; 'CST Mono'}


clear data
for i=1:length(sampleNames)
    
    fileName=strcat(directory, sampleNames{i});
    [ ~, startPos, endPos, overlappingReads, overlappingBases, featureLength, annotationType ] = ReadBedtoolsCoverage( fileName );
    data(:,i)=overlappingReads;
    
end

annTypes={'E'; 'TSS'}



fig=figure()
subplot(1,3,1)
mono=data(:,3);
poly=data(:,2);

rE=corr(mono(strcmp(annotationType, 'E')), poly(strcmp(annotationType, 'E')))^2;
rET=sprintf('%.2f',rE)

rTSS=corr(mono(strcmp(annotationType, 'TSS')), poly(strcmp(annotationType, 'TSS')))^2;
rTSST=sprintf('%.2f',rTSS)
rAll=sprintf('%.2f',corr(mono,poly)^2);

loglog(mono, poly, '.', 'Color', [0.7 0.7 0.7])
hold on
loglog(mono(strcmp(annotationType, 'TSS')), poly(strcmp(annotationType, 'TSS')), '.', 'Color', [ 0.0588    0.2045    0.8090])
loglog(mono(strcmp(annotationType, 'E')), poly(strcmp(annotationType, 'E')), '.', 'Color', [0.1231    0.6863    0.2627] )

xlabel('Monoclonal H3K27ac (CST)')
ylabel('Polyclonal H3K27ac (Active Motif)')
axis([0 10^4 0 10^4])
plotBeef(10)
h=legend(strcat('All Points R^2=', rAll), strcat('TSS R^2=', rTSST), strcat('Enhancers R^2=', rET), 'Location', 'NorthWest')
set(h,'FontSize',9);
%text(7000, 2, strcat('R^2 (all points)=', rAll),'HorizontalAlignment', 'right' )


subplot(1,3,2)
mono=data(:,3);
poly=data(:,1);


rE=corr(mono(strcmp(annotationType, 'E')), poly(strcmp(annotationType, 'E')))^2;
rET=sprintf('%.2f',rE)

rTSS=corr(mono(strcmp(annotationType, 'TSS')), poly(strcmp(annotationType, 'TSS')))^2;
rTSST=sprintf('%.2f',rTSS);
rAll=sprintf('%.2f',corr(mono,poly)^2);

loglog(mono, poly, '.', 'Color', [0.7 0.7 0.7])
hold on
loglog(mono(strcmp(annotationType, 'TSS')), poly(strcmp(annotationType, 'TSS')), '.', 'Color', [ 0.0588    0.2045    0.8090])
loglog(mono(strcmp(annotationType, 'E')), poly(strcmp(annotationType, 'E')), '.', 'Color', [0.1231    0.6863    0.2627] )

xlabel('Monoclonal H3K27ac (CST)')
ylabel('Monoclonal H3K27ac (Active Motif)')
%title('H3K27ac In HeLa Cells',  'FontSize',10)
axis([0 10^4 0 10^4])
plotBeef(10)
h=legend(strcat('All Points R^2=', rAll), strcat('TSS R^2=', rTSST), strcat('Enhancers R^2=', rET), 'Location', 'NorthWest')
set(h,'FontSize',9);
%text(7000, 2, strcat('R^2 (all)=', rAll),'HorizontalAlignment', 'right' )



subplot(1,3,3)
mono=data(:,1);
poly=data(:,2);

rE=corr(mono(strcmp(annotationType, 'E')), poly(strcmp(annotationType, 'E')))^2;
rET=sprintf('%.2f',rE)

rTSS=corr(mono(strcmp(annotationType, 'TSS')), poly(strcmp(annotationType, 'TSS')))^2;
rTSST=sprintf('%.2f',rTSS)
rAll=sprintf('%.2f',corr(mono,poly)^2)


loglog(mono, poly, '.', 'Color', [0.7 0.7 0.7])
hold on
loglog(mono(strcmp(annotationType, 'TSS')), poly(strcmp(annotationType, 'TSS')), '.', 'Color', [ 0.0588    0.2045    0.8090])

loglog(mono(strcmp(annotationType, 'E')), poly(strcmp(annotationType, 'E')), '.', 'Color', [0.1231    0.6863    0.2627] )


xlabel('Monoclonal H3K27ac (Active Motif)')
ylabel('Polyclonal H3K27ac (Active Motif)')
axis([0 10^4 0 10^4])
plotBeef(10)
h=legend(strcat('All R^2=', rAll), strcat('TSS R^2=', rTSST), strcat('Enhancers R^2=', rET), 'Location', 'NorthWest')
set(h,'FontSize',9);
%text(7000, 2, strcat('R^2 (all)=', rAll),'HorizontalAlignment', 'right' )



    fig.PaperUnits = 'inches';
    fig.PaperPosition = [0.6 0.6 10.5 2.9];
    fig.PaperPositionMode = 'manual';
    fig.PaperOrientation = 'landscape';
    print('C:\Data\ChipSeqPaper\Figures\scatterH3K27acByTSSvsEnhancer_HELA.tif','-dtiff','-r350')
    

    fig.PaperUnits = 'inches';
    fig.PaperPosition = [0.6 0.6 10.5 2.9];
    fig.PaperPositionMode = 'manual';
    fig.PaperOrientation = 'landscape';
   print('C:\Data\ChipSeqPaper\Figures\scatterH3K27acByTSSvsEnhancer_HELA.pdf','-dpdf','-r350')
    
 