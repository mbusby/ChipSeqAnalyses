clear
close all
directory='C:\Data\ChipSeqPaper\Data\RandomlyDownsampled\BamsIntersectedWith2000Bp\'
groupNames={'Mono'; 'Poly'}
BTLblue=[15 121 191]/255
plotCtr=1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure

sampleNames={

'mAb-K27ac-1.dupsRemoved.bam.bam.13407899.bam.txt'
'mAb-K27ac-2.dupsRemoved.bam.bam.13407899.bam.txt'
'pAb-K27ac-1.dupsRemoved.bam.bam.13407899.bam.txt'
'pAb-K27ac-2.dupsRemoved.bam.bam.13407899.bam.txt'

};

group=[1 1 2 2];


clear data
for i=1:length(sampleNames)
    
    fileName=strcat(directory, sampleNames{i});
    [ chr, startPos, endPos, overlappingReads, overlappingBases, featureLength ] = ReadBedtoolsCoverage( fileName);
    data(:,i)=overlappingReads;
    
end


clear fig
fig=figure()
plotCtr=1;
subplot(1,3,plotCtr)
loglog(data(:,1), data(:,2),'.', 'Color', [0.6 0.6 0.6], 'MarkerSize', 7)
r2=corr(data(:,1), data(:,2))^2;
text(1.5, 60000, strcat('R^2=',num2str(r2,2)), 'FontSize', 10) 
xlabel('Monoclonal Rep 1', 'FontSize',12)
ylabel('Monoclonal Rep 2', 'FontSize',12)
axis([0,  10^5, 0,  10^5])

set(gca,'XTick',[1 10 10^2 10^3 10^4 10^5])
xt = get(gca, 'XTick');
set(gca, 'FontSize', 9)
set(gca,'YTick',[1 10 10^2 10^3 10^4 10^5])
xt = get(gca, 'YTick');
set(gca, 'FontSize', 9)

plotCtr=plotCtr+1;

firstPoly=find(group==2, 1);
subplot(1,3,plotCtr)
loglog(data(:,1), data(:,firstPoly), '.', 'Color', BTLblue, 'MarkerSize', 7);
r2=corr(data(:,1), data(:,firstPoly))^2;
text(1.5, 60000, strcat('R^2=',num2str(r2,2)), 'FontSize', 10) 
title('H3K27ac','FontSize', 14)
xlabel('Monoclonal Rep 1', 'FontSize',12)
ylabel('Polyclonal Rep 1', 'FontSize',12)
axis([0,  10^5, 0,  10^5])
plotCtr=plotCtr+1;
set(gca,'XTick',[1 10 10^2 10^3 10^4 10^5])
xt = get(gca, 'XTick');
set(gca, 'FontSize', 9)
set(gca,'YTick',[1 10 10^2 10^3 10^4 10^5])
xt = get(gca, 'YTick');
set(gca, 'FontSize', 9)

subplot(1,3,plotCtr)
loglog(data(:,firstPoly), data(:,firstPoly+1),'.', 'Color', [0.6 0.6 0.6], 'MarkerSize', 7)
r2=corr(data(:,firstPoly), data(:,firstPoly+1))^2;
text(1.5, 60000, strcat('R^2=',num2str(r2,2)), 'FontSize', 10) 
xlabel('Polyclonal Rep 1', 'FontSize', 12)
ylabel('Polyclonal Rep 2', 'FontSize', 12)
axis([0,  10^5, 0,  10^5])
set(gca,'XTick',[1 10 10^2 10^3 10^4 10^5])
xt = get(gca, 'XTick');
set(gca, 'FontSize', 9)
set(gca,'YTick',[1 10 10^2 10^3 10^4 10^5])
xt = get(gca, 'YTick');
set(gca, 'FontSize', 9)
plotCtr=plotCtr+1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




    fig.PaperUnits = 'inches';
    fig.PaperPosition = [0.6 0.6 10.5 2.9];
    fig.PaperPositionMode = 'manual';
    fig.PaperOrientation = 'landscape';
   print('C:\Data\ChipSeqPaper\Figures\random_scatterH3K27ac','-dpng','-r350')
    


sampleNames={

'1_Exp4-K562-2A-K27me3-mono_H9V71ADXX.1.aligned.duplicates_marked.dupsRemoved.bam.bam.9356733.bam.txt'
'1_Exp4-K562-2B-K27me3-mono_H9V71ADXX.1.aligned.duplicates_marked.dupsRemoved.bam.bam.9356733.bam.txt'
'1_Exp4-K562-2C-K27me3-mono_H9V71ADXX.1.aligned.duplicates_marked.dupsRemoved.bam.bam.9356733.bam.txt'
'1_Exp4-K562-2D-K27me3-mono_H9V71ADXX.1.aligned.duplicates_marked.dupsRemoved.bam.bam.9356733.bam.txt'
'1_Exp4-K562-2E-K27me3-poly_H9V71ADXX.1.aligned.duplicates_marked.dupsRemoved.bam.bam.9356733.bam.txt'
'1_Exp4-K562-2F-K27me3-poly_H9V71ADXX.1.aligned.duplicates_marked.dupsRemoved.bam.bam.9356733.bam.txt'
'1_Exp4-K562-2G-K27me3-poly_H9V71ADXX.1.aligned.duplicates_marked.dupsRemoved.bam.bam.9356733.bam.txt'
'1_Exp4-K562-2H-K27me3-poly_H9V71ADXX.1.aligned.duplicates_marked.dupsRemoved.bam.bam.9356733.bam.txt'

};

group=[1 1 1 1 2 2 2 2];


clear data
for i=1:length(sampleNames)
    
    fileName=strcat(directory, sampleNames{i});
    [ chr, startPos, endPos, overlappingReads, overlappingBases, featureLength ] = ReadBedtoolsCoverage( fileName );
    data(:,i)=overlappingReads;
    
end

fig=figure()
plotCtr=1;

subplot(1,3,plotCtr)
loglog(data(:,1), data(:,2),'.', 'Color', [0.6 0.6 0.6], 'MarkerSize', 7)
r2=corr(data(:,1), data(:,2))^2;
text(1.5, 60000, strcat('R^2=',num2str(r2,2)), 'FontSize', 10) 
xlabel('Monoclonal Rep 1', 'FontSize', 12)
ylabel('Monoclonal Rep 2', 'FontSize', 12)
axis([0,  10^5, 0,  10^5])
set(gca,'XTick',[1 10 10^2 10^3 10^4 10^5])
set(gca,'YTick',[1 10 10^2 10^3 10^4 10^5])
xt = get(gca, 'YTick');
set(gca, 'FontSize', 9)
plotCtr=plotCtr+1;

firstPoly=find(group==2, 1);
subplot(1,3,plotCtr)
loglog(data(:,1), data(:,firstPoly), '.', 'Color', BTLblue, 'MarkerSize', 7);
r2=corr(data(:,1), data(:,firstPoly))^2;
text(1.5, 60000, strcat('R^2=',num2str(r2,2)), 'FontSize', 10) 
title('H3K27me3', 'FontSize', 14)
xlabel('Monoclonal Rep 1', 'FontSize', 12)
ylabel('Polyclonal Rep 1', 'FontSize', 12)
axis([0,  10^5, 0,  10^5])
set(gca,'XTick',[1 10 10^2 10^3 10^4 10^5])
set(gca,'YTick',[1 10 10^2 10^3 10^4 10^5])
xt = get(gca, 'YTick');
set(gca, 'FontSize', 9)
plotCtr=plotCtr+1;

subplot(1,3,plotCtr)
loglog(data(:,firstPoly), data(:,firstPoly+1),'.', 'Color', [0.6 0.6 0.6], 'MarkerSize', 7)
r2=corr(data(:,firstPoly), data(:,firstPoly+1))^2;
text(1.5, 60000, strcat('R^2=',num2str(r2,2)), 'FontSize', 10) 
xlabel('Polyclonal Rep 1', 'FontSize', 12)
ylabel('Polyclonal Rep 2', 'FontSize', 12)
%plotBeef(16)
axis([0,  10^5, 0,  10^5])
set(gca,'XTick',[1 10 10^2 10^3 10^4 10^5])
plotCtr=plotCtr+1;
set(gca,'YTick',[1 10 10^2 10^3 10^4 10^5])
xt = get(gca, 'YTick');
set(gca, 'FontSize', 9)



    fig.PaperUnits = 'inches';
    fig.PaperPosition = [0.6 0.6 10.5 2.9];
    fig.PaperPositionMode = 'manual';
    fig.PaperOrientation = 'landscape';
   print('C:\Data\ChipSeqPaper\Figures\random_scatterH3K27me3','-dpng','-r350')
    

    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


groupNames={'Mono'; 'Poly'}

sampleNames={
'K562_K4me1_Ido.dupsRemoved.bam.bam.20565471.bam.txt'  
'K562_K4me1_LizExp11.dupsRemoved.bam.bam.20565471.bam.txt' 
'mAb-K4me1-1.dupsRemoved.bam.bam.20565471.bam.txt'
'pAb-K4me1-1.dupsRemoved.bam.bam.20565471.bam.txt'
'pAb-K4me1-2.dupsRemoved.bam.bam.20565471.bam.txt'
};



group=[1 1 1  2 2];


clear data
for i=1:length(sampleNames)
    
    fileName=strcat(directory, sampleNames{i});
    [ chr, startPos, endPos, overlappingReads, overlappingBases, featureLength ] = ReadBedtoolsCoverage( fileName );
    data(:,i)=overlappingReads;
    
end

fig=figure()
plotCtr=1;

subplot(1,3,plotCtr)
loglog(data(:,1), data(:,2),'.', 'Color', [0.6 0.6 0.6], 'MarkerSize', 7)
hold on
r2=corr(data(:,1), data(:,2))^2;
text(1.5, 60000, strcat('R^2=',num2str(r2,2)), 'FontSize', 10) 
xlabel('Monoclonal Rep 1', 'FontSize', 12)
ylabel('Monoclonal Rep 2', 'FontSize', 12)
axis([0,  10^5, 0,  10^5])
set(gca,'XTick',[1 10 10^2 10^3 10^4 10^5])
plotCtr=plotCtr+1;
set(gca,'YTick',[1 10 10^2 10^3 10^4 10^5])
xt = get(gca, 'YTick');
set(gca, 'FontSize', 9)

firstPoly=find(group==2, 1);
subplot(1,3,plotCtr)
loglog(data(:,1), data(:,firstPoly), '.', 'Color', BTLblue, 'MarkerSize', 7);
r2=corr(data(:,1), data(:,firstPoly))^2;
text(1.5, 60000, strcat('R^2=',num2str(r2,2)), 'FontSize', 10) 
title('H3K4me1', 'FontSize', 14)
xlabel('Monoclonal Rep 1', 'FontSize', 12)
ylabel('Polyclonal Rep 1', 'FontSize', 12)
axis([0,  10^5, 0,  10^5])
set(gca,'XTick',[1 10 10^2 10^3 10^4 10^5])
plotCtr=plotCtr+1;
set(gca,'YTick',[1 10 10^2 10^3 10^4 10^5])
xt = get(gca, 'YTick');
set(gca, 'FontSize', 9)

subplot(1,3,plotCtr)
loglog(data(:,firstPoly), data(:,firstPoly+1),'.', 'Color', [0.6 0.6 0.6], 'MarkerSize', 7)
r2=corr(data(:,firstPoly), data(:,firstPoly+1))^2;
text(1.5, 60000, strcat('R^2=',num2str(r2,2)), 'FontSize', 10) 
xlabel('Polyclonal Rep 1', 'FontSize', 12)
ylabel('Polyclonal Rep 2', 'FontSize', 12)
%plotBeef(16)
axis([0,  10^5, 0,  10^5])
set(gca,'XTick',[1 10 10^2 10^3 10^4 10^5])
plotCtr=plotCtr+1;
set(gca,'YTick',[1 10 10^2 10^3 10^4 10^5])
xt = get(gca, 'YTick');
set(gca, 'FontSize', 9)

    fig.PaperUnits = 'inches';
    fig.PaperPosition = [0.6 0.6 10.5 2.9];
    fig.PaperPositionMode = 'manual';
    fig.PaperOrientation = 'landscape';
   print('C:\Data\ChipSeqPaper\Figures\random_scatterH3K4me1.png','-dpng','-r350')





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




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

group=[1 1 1 1 2 2 2 2];


clear data
for i=1:length(sampleNames)
    
    fileName=strcat(directory, sampleNames{i});
    [ chr, startPos, endPos, overlappingReads, overlappingBases, featureLength ] = ReadBedtoolsCoverage( fileName );
    data(:,i)=overlappingReads;
    
end

fig=figure()
plotCtr=1;

subplot(1,3,plotCtr)
loglog(data(:,1), data(:,2),'.', 'Color', [0.6 0.6 0.6], 'MarkerSize', 7)
r2=corr(data(:,1), data(:,2))^2;
text(1.5, 60000, strcat('R^2=',num2str(r2,2)), 'FontSize', 10) 
xlabel('Monoclonal Rep 1', 'FontSize', 12)
ylabel('Monoclonal Rep 2', 'FontSize', 12)
axis([0,  10^5, 0,  10^5])
set(gca,'XTick',[1 10 10^2 10^3 10^4 10^5])
plotCtr=plotCtr+1;
set(gca,'YTick',[1 10 10^2 10^3 10^4 10^5])
xt = get(gca, 'YTick');
set(gca, 'FontSize', 9)

firstPoly=find(group==2, 1);
subplot(1,3,plotCtr)
loglog(data(:,1), data(:,firstPoly), '.', 'Color', BTLblue, 'MarkerSize', 7);
r2=corr(data(:,1), data(:,firstPoly))^2;
text(1.5, 60000, strcat('R^2=',num2str(r2,2)), 'FontSize', 10) 
title('H3K4me3', 'FontSize', 14)
xlabel('Monoclonal Rep 1', 'FontSize', 12)
ylabel('Polyclonal Rep 1', 'FontSize', 12)
axis([0,  10^5, 0,  10^5])
set(gca,'XTick',[1 10 10^2 10^3 10^4 10^5])
plotCtr=plotCtr+1;
set(gca,'YTick',[1 10 10^2 10^3 10^4 10^5])
xt = get(gca, 'YTick');
set(gca, 'FontSize', 9)

subplot(1,3,plotCtr)
loglog(data(:,firstPoly), data(:,firstPoly+1),'.', 'Color', [0.6 0.6 0.6], 'MarkerSize', 7)
r2=corr(data(:,firstPoly), data(:,firstPoly+1))^2;
text(1.5, 60000, strcat('R^2=',num2str(r2,2)), 'FontSize', 10) 
xlabel('Polyclonal Rep 1', 'FontSize', 12)
ylabel('Polyclonal Rep 2', 'FontSize', 12)
%plotBeef(16)
axis([0,  10^5, 0,  10^5])
set(gca,'XTick',[1 10 10^2 10^3 10^4 10^5])
plotCtr=plotCtr+1;
set(gca,'YTick',[1 10 10^2 10^3 10^4 10^5])
xt = get(gca, 'YTick');
set(gca, 'FontSize', 9)


    fig.PaperUnits = 'inches';
    fig.PaperPosition = [0.6 0.6 10.5 2.9];
    fig.PaperPositionMode = 'manual';
    fig.PaperOrientation = 'landscape';
   print('C:\Data\ChipSeqPaper\Figures\random_scatterH3K4me3.png','-dpng','-r350')





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


groupNames={'Mono'; 'Poly'}

sampleNames={
'K562_K9me3_Ido.dupsRemoved.bam.bam.18740026.bam.txt'
'K562_K9me3_LizExp11.dupsRemoved.bam.bam.18740026.bam.txt'
'mAb-K9me3-2.dupsRemoved.bam.bam.18740026.bam.txt'
'pAb-K9me3-1.dupsRemoved.bam.bam.18740026.bam.txt'
'pAb-K9me3-2.dupsRemoved.bam.bam.18740026.bam.txt'
};



group=[1 1 1 2 2];


clear data
for i=1:length(sampleNames)
    
    fileName=strcat(directory, sampleNames{i});
    [ chr, startPos, endPos, overlappingReads, overlappingBases, featureLength ] = ReadBedtoolsCoverage( fileName );
    data(:,i)=overlappingReads;
    
end

fig=figure()
plotCtr=1;

subplot(1,3,plotCtr)
loglog(data(:,1), data(:,2),'.', 'Color', [0.6 0.6 0.6], 'MarkerSize', 7)
r2=corr(data(:,1), data(:,2))^2;
text(1.5, 60000, strcat('R^2=',num2str(r2,2)), 'FontSize', 10)
xlabel('Monoclonal Rep 1', 'FontSize', 12)
ylabel('Monoclonal Rep 2', 'FontSize', 12)
axis([0,  10^5, 0,  10^5])
set(gca,'XTick',[1 10 10^2 10^3 10^4 10^5])
plotCtr=plotCtr+1;
set(gca,'YTick',[1 10 10^2 10^3 10^4 10^5])
xt = get(gca, 'YTick');
set(gca, 'FontSize', 9)

firstPoly=find(group==2, 1);
subplot(1,3,plotCtr)
loglog(data(:,1), data(:,firstPoly), '.', 'Color', BTLblue, 'MarkerSize', 7);
r2=corr(data(:,1), data(:,firstPoly))^2;
text(1.5, 60000, strcat('R^2=',num2str(r2,2)), 'FontSize', 10) 
title('H3K9me3', 'FontSize', 14)
xlabel('Monoclonal Rep 1', 'FontSize', 12)
ylabel('Polyclonal Rep 1', 'FontSize', 12)
axis([0,  10^5, 0,  10^5])
set(gca,'XTick',[1 10 10^2 10^3 10^4 10^5])
plotCtr=plotCtr+1;
set(gca,'YTick',[1 10 10^2 10^3 10^4 10^5])
xt = get(gca, 'YTick');
set(gca, 'FontSize', 9)

subplot(1,3,plotCtr)
loglog(data(:,firstPoly), data(:,firstPoly+1),'.', 'Color', [0.6 0.6 0.6], 'MarkerSize', 7)
r2=corr(data(:,firstPoly), data(:,firstPoly+1))^2;
text(1.5, 60000, strcat('R^2=',num2str(r2,2)), 'FontSize', 10) 
xlabel('Polyclonal Rep 1', 'FontSize', 12)
ylabel('Polyclonal Rep 2', 'FontSize', 12)
%plotBeef(16)
axis([0,  10^5, 0,  10^5])
set(gca,'XTick',[1 10 10^2 10^3 10^4 10^5])
plotCtr=plotCtr+1;
set(gca,'YTick',[1 10 10^2 10^3 10^4 10^5])
xt = get(gca, 'YTick');
set(gca, 'FontSize', 9)


    fig.PaperUnits = 'inches';
    fig.PaperPosition = [0.6 0.6 10.5 2.9];
    fig.PaperPositionMode = 'manual';
    fig.PaperOrientation = 'landscape';
    print('C:\Data\ChipSeqPaper\Figures\random_scatterH3K9me3.png','-dpng','-r350')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

directory='C:\Temp\MonoVsPolyPaper\IntersectedWithCombinedMarksHeLa\'
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
    [ chr, startPos, endPos, overlappingReads, overlappingBases, featureLength, annotationType ] = ReadBedtoolsCoverage( fileName );
    data(:,i)=overlappingReads;
    
end

%getPCA( data, group, groupNames )
%clusterSamples( data, sampleNames )
%getAllLogLogPlots( data,sampleNames )


[h pTest]=ttest2(data(:, 1:2)' , data(:, 2:4)');
sum(pTest<0.05)/length(pTest)

[h pNull]=ttest2(data(:, [1 3 ] )' , data(:, [2 4 ])');
sum(pNull<0.05)/length(pNull)


annTypes=unique(annotationType);

g=unique(group);


sampleNamesShort={
'K27ac Active Motif Monoclonal'
'K27ac Active Motif Polyclonal'
'K27ac CST Monoclonal'
};

annTypes={'E'; 'TSS'}

   for i=1
    for j=3
        for k=1:length(annTypes)
            for m=1:length(annTypes)
                if k~=m
                    figure 
                    good=strcmp(annotationType, annTypes{m});
                    loglog(data(good, i), data(good,j), '.')                
                    hold on
                    good=strcmp(annotationType, annTypes{k});
                    loglog(data(good, i), data(good,j), '.')
                    xlabel(sampleNamesShort{i})
                    ylabel(sampleNamesShort{j})
                    at{1}=annTypes{m};
                    at{2}=annTypes{k};
                    legend(at)
                end
            end
        end       
        
    end
 end

 


%%%%5%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





