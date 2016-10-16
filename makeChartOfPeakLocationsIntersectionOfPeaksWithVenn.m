clear
    dirName ='C:\Data\ChipSeqPaper\Data\MergedPeaks'
    outFileName='C:\Data\ChipSeqPaper\Data\MergedPeaks\EncodeVsReads.png'
    
    annTypes ={
    'CTCF'  %1
    'PF'    %2
    'TSS'   %3
    'T'     %4
    'E'     %5
    'WE'    %6   
    'R'     %7
    };

 annTypesLong ={
    'CTCF'  %1
    'Promoter'    %2
    'Tx Start'   %3
    'Transcribed'     %4
    'Enhancer'     %5
    'Weak Enchancer'    %6   
    'Repressed'     %7
    };

    colorscheme=[ 
         0.8549    0.6784    0.5255  %  CTCF 
         1.0000    0.4000    0.6980; %  PF
         0.0588    0.4745    0.7490; %  TSS 
         0.5294    0.8078    0.9804; %  T
         0.2431    0.6863    0.4627; %  E
         0.5882    0.9765    0.4824; %  WE
         1.0000    0.3       0.3;    %  R
    ];

 colorsChart=[ 
        
         0.5294    0.1412    0.2039; %  maroon
          0.7020    0.7020    0.7020; %  T
           0.5294    0.1412    0.2039; %  maroon
          0.7020    0.7020    0.7020; %  T
         
           0.5294    0.1412    0.2039; %  maroon
          0.7020    0.7020    0.7020; %  T
         
           0.5294    0.1412    0.2039; %  maroon
          0.7020    0.7020    0.7020; %  T
         
           0.5294    0.1412    0.2039; %  maroon
          0.7020    0.7020    0.7020; %  T
         
         
    ];

    ctrVals=1;
    
    myBeds = dir(dirName) %yields vector of structures with myBeds.name is the name
    
    

files= dir(dirName);

for i=3:length(files)
    fileNames{i-2}=files(i).name;
end

fileNames={
    'K27acMono_K27acPoly.bed.txt'
    'K27acMono_notK27acPoly.bed.txt'
    'K27acPolyNotMono.bed.txt'

    'K4me1Mono_K4me1Poly.bed.txt'
    'K4me1Mono_notK4me1Poly.bed.txt'
    'K4me1_polyNotMono.bed.txt'
    
    'K4me3Mono_K4me3Poly.bed.txt'
    'K4me3Mono_notK4me3Poly.bed.txt'
    'K4me3_polyNotMono.bed.txt'
  }

for i=1:length(fileNames)
shortNames{i}=strrep(strrep(strrep(strrep(fileNames{i}, '.bam', ''), '.txt', ''), 'out', ''), '_', ' ');
end

    shortNamesSorted=shortNames;
    
shortNames={    'H3K27ac Mono and Poly'
    'H3K27ac Mono Only'
    'H3K27ac Poly Only'
    'H3K4me1 Mono and Poly'
    'H3K4me1 Mono Only'
    'H3K4me1 Poly Only'
    'H3K4me3 Mono and Poly'
    'H3K4me3 Mono Only'
    'H3K4me3 Poly Only'
    }

antibodyNames={    'H3K27ac'
    'H3K27ac'
    'H3K27ac'

    'H3K4me1'
    'H3K4me1'
    'H3K4me1'
    'H3K4me3'
    'H3K4me3'
    'H3K4me3'
    }


chromosomes=[];

for i=1:length(fileNames)
     fileName=strcat(dirName,filesep, fileNames{i});
     fid=fopen(strcat(dirName,filesep, fileNames{i}), 'r');    
  
       
     importedData = textscan( fid, '%s %f %f %s %f %f %s %s %s %s %s  %s %f', 'Headerlines', 0, 'delimiter', '\t' );
     
     chromosomes=unique([chromosomes; unique(importedData{1})]);
              
     fclose(fid);
end    



for i=1:length(fileNames)
     fileName=strcat(dirName,filesep, fileNames{i})
          fid=fopen(fileName, 'r');      
     
     importedData = textscan( fid, '%s %f %f %s %f %f %s %s %s %s %s  %s %f', 'Headerlines', 0, 'delimiter', '\t' );
     
     chromosome=importedData{1};
     startPos=importedData{2};
     endPos=importedData{3};   
     featureType=importedData{7};
     overlappingBases=importedData{13};
    
     total=sum(overlappingBases);
     for j=1:length(annTypes)
         overlapping(ctrVals,j)=sum(overlappingBases(strcmp(featureType, annTypes{j})));
         overlappingPerc(ctrVals,j)=sum(overlappingBases(strcmp(featureType, annTypes{j})))/total;
     end

      yTickLabels{ctrVals}=shortNames{i};

      ctrVals=ctrVals+1;
      
 end    

    figure
    barh(overlappingPerc*100, 'stacked')
     colormap(colorscheme)

    set(gca, 'YTick', 1:ctrVals-1)
    set(gca, 'YTickLabel', yTickLabels)
    set(gca, 'TickLength', [0 0]);
    set(gca,'View',[0 -90])    
     axis([0 100.3 -.025  size(overlappingPerc,1)+1])
   
    %legend(annTypes, 'location', 'SouthOutside','Orientation','horizontal')    

    set(gca,'fontsize',9);
    title('ENCODE Regions of Bases Mapping To Peaks')
 
    display('You may need to turn it right side up using set(gca,View,[0 -90]) again ...')
    xlabel('Percentage of Annotated Genomic Bases')
    legend(annTypes, 'Orientation', 'horizontal', 'Location', 'SouthOutside')

 
 

goodOverlaps(1:3,1)=sum(overlapping(1:3, [3,5]),2)

goodOverlaps(4:6,1)=sum(overlapping(4:6, [5]),2)

goodOverlaps(7:9,1)=sum(overlapping(7:9, [3]),2)

goodOverlaps(:,2)=sum(overlapping,2)-goodOverlaps(:,1);

 allGood=[goodOverlaps(2,:) goodOverlaps(1,:) goodOverlaps(3,:); goodOverlaps(5,:) goodOverlaps(4,:) goodOverlaps(6,:);   ...
     goodOverlaps(8,:) goodOverlaps(7,:) goodOverlaps(9,:);  ]


 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
bothRow=1
monoRow=2
polyRow=3

mono=sum(goodOverlaps(bothRow,:))+sum(goodOverlaps(monoRow,:));
poly=sum(goodOverlaps(polyRow,:))+sum(goodOverlaps(bothRow,:));
both=sum(goodOverlaps(bothRow,:));

figure
barh(1, goodOverlaps(monoRow,1), 'facecolor', 'r')
hold on
barh(2, goodOverlaps(bothRow,1), 'facecolor', [  0.4000    0.3020    0.6000])
barh(3, goodOverlaps(polyRow,1), 'facecolor',[  0.4000    0.7020    1.0000])
barh(1,0-goodOverlaps(monoRow,2), 'facecolor', 'r')
barh(2,0-goodOverlaps(bothRow,2), 'facecolor', [ 0.4000    0.3020    0.6000])
barh(3,0-goodOverlaps(polyRow,2), 'facecolor', [ 0.4000    0.7020    1.0000])

set(gca, 'YTick', 1:3)
    set(gca, 'YTickLabel', { 'Mono Only', 'Mono and Poly', 'Poly only'})
    set(gca, 'TickLength', [0 0]); 
title(antibodyNames{bothRow})
set(gca, 'XTick',[-30000000 0 30000000])
set(gca, 'XTickLabel', {'30,000,000'; '0' ;'30,000,000'})
xlabel('Bases Called In Peaks')
    set(gca,'View',[0 -90])  
axis([-30500000 30500000 -0.5 3.75])
text( -500000, 0.25, 'Other Regions', 'HorizontalAlignment', 'right')
text( 500000, 0.25, 'Expected Regions', 'HorizontalAlignment', 'left')
plot([0 0], [-1 4], '-', 'Color', 'k')
title(antibodyNames{bothRow})
plotBeef(12)


figure



A = [mono poly ]; I = [both ]
venn(A,I,'FaceColor',{'r',[0 .5 1]},'FaceAlpha',{1,0.6},'EdgeColor','black')
set(gca,'xtick',[])  
set(gca,'ytick',[])  
%text(-3000, 2500, 'Mono', 'HorizontalAlignment', 'left')
%text(3000, 2500, 'Poly', 'HorizontalAlignment', 'right')
title(antibodyNames{bothRow})
plotBeef(14)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

bothRow=4
monoRow=5
polyRow=6

mono=sum(goodOverlaps(bothRow,:))+sum(goodOverlaps(monoRow,:));
poly=sum(goodOverlaps(polyRow,:))+sum(goodOverlaps(bothRow,:));
both=sum(goodOverlaps(bothRow,:));

figure
barh(1, goodOverlaps(monoRow,1), 'facecolor', 'r')
hold on
barh(2, goodOverlaps(bothRow,1), 'facecolor', [  0.4000    0.3020    0.6000])
barh(3, goodOverlaps(polyRow,1), 'facecolor',[  0.4000    0.7020    1.0000])
barh(1,0-goodOverlaps(monoRow,2), 'facecolor', 'r')
barh(2,0-goodOverlaps(bothRow,2), 'facecolor', [ 0.4000    0.3020    0.6000])
barh(3,0-goodOverlaps(polyRow,2), 'facecolor', [ 0.4000    0.7020    1.0000])

set(gca, 'YTick', 1:3)
    set(gca, 'YTickLabel', { 'Mono Only', 'Mono and Poly', 'Poly only'})
    set(gca, 'TickLength', [0 0]); 
title(antibodyNames{bothRow})
set(gca, 'XTick',[-20000000 0 20000000])
set(gca, 'XTickLabel', {'20,000,000'; '0' ;'20,000,000'})
xlabel('Bases Called In Peaks')
    set(gca,'View',[0 -90])  
axis([-20000000 20000000  -0.5 3.75])
text( -500000, 0.25, 'Other Regions', 'HorizontalAlignment', 'right')
text( 500000, 0.25, 'Expected Regions', 'HorizontalAlignment', 'left')
plot([0 0], [-1 4], '-', 'Color', 'k')
title(antibodyNames{bothRow})
plotBeef(12)


figure



A = [mono poly ]; I = [both ]
venn(A,I,'FaceColor',{'r',[0 .5 1]},'FaceAlpha',{1,0.6},'EdgeColor','black')
set(gca,'xtick',[])  
set(gca,'ytick',[])  
%text(-3000, 2500, 'Mono', 'HorizontalAlignment', 'left')
%text(3000, 2500, 'Poly', 'HorizontalAlignment', 'right')
title(antibodyNames{bothRow})
plotBeef(14)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

bothRow=7
monoRow=8
polyRow=9

mono=sum(goodOverlaps(bothRow,:))+sum(goodOverlaps(monoRow,:));
poly=sum(goodOverlaps(polyRow,:))+sum(goodOverlaps(bothRow,:));
both=sum(goodOverlaps(bothRow,:));

figure
barh(1, goodOverlaps(monoRow,1), 'facecolor', 'r')
hold on
barh(2, goodOverlaps(bothRow,1), 'facecolor', [  0.4000    0.3020    0.6000])
barh(3, goodOverlaps(polyRow,1), 'facecolor',[  0.4000    0.7020    1.0000])
barh(1,0-goodOverlaps(monoRow,2), 'facecolor', 'r')
barh(2,0-goodOverlaps(bothRow,2), 'facecolor', [ 0.4000    0.3020    0.6000])
barh(3,0-goodOverlaps(polyRow,2), 'facecolor', [ 0.4000    0.7020    1.0000])

set(gca, 'YTick', 1:3)
    set(gca, 'YTickLabel', { 'Mono Only', 'Mono and Poly', 'Poly only'})
    set(gca, 'TickLength', [0 0]); 
title(antibodyNames{bothRow})
set(gca, 'XTick',[-22000000 0 22000000])
set(gca, 'XTickLabel', {'22,000,000'; '0' ;'22,000,000'})
xlabel('Bases Called In Peaks')
    set(gca,'View',[0 -90])  
axis([-22000000 22000000  -0.5 3.75])
text( -500000, 0.25, 'Other Regions', 'HorizontalAlignment', 'right')
text( 500000, 0.25, 'Expected Regions', 'HorizontalAlignment', 'left')
plot([0 0], [-1 4], '-', 'Color', 'k')
plotBeef(12)


figure



A = [mono poly ]; I = [both ]
venn(A,I,'FaceColor',{'r',[0 .5 1]},'FaceAlpha',{1,0.6},'EdgeColor','black')
set(gca,'xtick',[])  
set(gca,'ytick',[])  
%text(-3000, 2500, 'Mono', 'HorizontalAlignment', 'left')
%text(3000, 2500, 'Poly', 'HorizontalAlignment', 'right')
title(antibodyNames{bothRow})
plotBeef(14)