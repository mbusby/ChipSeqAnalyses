clear
    dirName ='C:\Data\ChipSeqPaper\Data\DownsampledByInsertSize\PeaksIntersectedWithEncodeAtBaseLevel\'
    outFileName='C:\Data\ChipSeqPaper\Figures\InsertSize_EncodeVsPeaks.png'
    
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

    ctrVals=1;
    
    myBeds = dir(dirName) %yields vector of structures with myBeds.name is the name
    
    fileNamesSorted={
'mAb-K27ac-1.normalized.bed.txt'
'mAb-K27ac-2.normalized.bed.txt'
'pAb-K27ac-1.normalized.bed.txt'
'pAb-K27ac-2.normalized.bed.txt'
'spacer'
'1_Exp4-K562-2A-K27me3-mono_H9V71ADXX.1.normalized.bed.txt'
'1_Exp4-K562-2B-K27me3-mono_H9V71ADXX.1.normalized.bed.txt'
'1_Exp4-K562-2C-K27me3-mono_H9V71ADXX.1.normalized.bed.txt'
'1_Exp4-K562-2D-K27me3-mono_H9V71ADXX.1.normalized.bed.txt'
'1_Exp4-K562-2E-K27me3-poly_H9V71ADXX.1.normalized.bed.txt'
'1_Exp4-K562-2F-K27me3-poly_H9V71ADXX.1.normalized.bed.txt'
'1_Exp4-K562-2G-K27me3-poly_H9V71ADXX.1.normalized.bed.txt'
'spacer'
'K562_K4me1_Ido.normalized.bed.txt'
'K562_K4me1_LizExp11.normalized.bed.txt'
'mAb-K4me1-1.normalized.bed.txt'
'pAb-K4me1-1.normalized.bed.txt'
'pAb-K4me1-2.normalized.bed.txt'
'spacer'
'1_Exp4-K562-1A-K4me3-mono_H9V71ADXX.1.normalized.bed.txt'
'1_Exp4-K562-1B-K4me3-mono_H9V71ADXX.1.normalized.bed.txt'
'1_Exp4-K562-1C-K4me3-mono_H9V71ADXX.1.normalized.bed.txt'
'1_Exp4-K562-1D-K4me3-mono_H9V71ADXX.1.normalized.bed.txt'
'1_Exp4-K562-1E-K4me3-poly_H9V71ADXX.1.normalized.bed.txt'
'1_Exp4-K562-1F-K4me3-poly_H9V71ADXX.1.normalized.bed.txt'
'1_Exp4-K562-1G-K4me3-poly_H9V71ADXX.1.normalized.bed.txt'
'1_Exp4-K562-1H-K4me3-poly_H9V71ADXX.1.normalized.bed.txt'
'spacer'
'K562_K9me3_Ido.normalized.bed.txt'
'K562_K9me3_LizExp11.normalized.bed.txt'
'mAb-K9me3-2.normalized.bed.txt'
'pAb-K9me3-1.normalized.bed.txt'
'pAb-K9me3-2.normalized.bed.txt'


};

shortNames={'H3K27ac Mono Rep 1'
'H3K27ac Mono Rep 2'
'H3K27ac Poly Rep 1'
'H3K27ac Poly Rep 1'
''
'H3K27me3 Mono Rep 1'
'H3K27me3 Mono Rep 2'
'H3K27me3 Mono Rep 3'
'H3K27me3 Mono Rep 4'
'H3K27me3 Poly Rep 1'
'H3K27me3 Poly Rep 2'
'H3K27me3 Poly Rep 3'
''
'H3K4me1 Mono Rep 1'
'H3K4me1 Mono Rep 2'
'H3K4me1 Mono Rep 3'
'H3K4me1 Poly Rep 1'
'H3K4me1 Poly Rep 2'
''
'H3K4me3 Mono Rep 1'
'H3K4me3 Mono Rep 2'
'H3K4me3 Mono Rep 3'
'H3K4me3 Mono Rep 4'
'H3K4me3 Poly Rep 1'
'H3K4me3 Poly Rep 2'
'H3K4me3 Poly Rep 3'
'H3K4me3 Poly Rep 4'
''
'H3K9me3 Mono Rep 1'
'H3K9me3 Mono Rep 2'
'H3K9me3 Mono Rep 3'
'H3K9me3 Poly Rep 1'
'H3K9me3 Poly Rep 2'


};
    
    
    
    for i=1:length(fileNamesSorted)
        
         if strcmp(fileNamesSorted{i}, 'spacer')==0

             fileName=strcat(dirName,filesep, fileNamesSorted{i});
             [ chr, startPos, endPos, featureType, nBases ] = ReadBedtoolsIntersectOfBedsWithReads( fileName );
             total=sum(nBases);
             for j=1:length(annTypes)
                peakTotals(ctrVals,j)=sum(nBases(strcmp(featureType, annTypes{j})));
                tablePerc(ctrVals,j)=sum(nBases(strcmp(featureType, annTypes{j})))/total;
                thisData(i,j)=sum(nBases(strcmp(featureType, annTypes{j})))/total;
              end

             yTickLabels{ctrVals}=shortNamesSorted{i};

             ctrVals=ctrVals+1;
         else
            tablePerc(ctrVals,:)=zeros(1, size(tablePerc,2)) ;
                       
            ctrVals=ctrVals+1;
                         
         end
                 
    end    

    figure
    barh(tablePerc*100, 'stacked')
     colormap(colorscheme)

    set(gca, 'YTick', 1:ctrVals-1)
    set(gca, 'YTickLabel', yTickLabels)
    set(gca, 'TickLength', [0 0]);
    set(gca,'View',[0 -90])    
     axis([0 100.3 -.025  size(tablePerc,1)+1])
   
    legend(annTypes, 'location', 'SouthOutside','Orientation','horizontal')    

    set(gca,'fontsize',9);
    title('% Peaks Mapping To ENCODE Genome Regions')
 
    display('You may need to turn it right side up using set(gca,View,[0 -90]) again ...')
    
    outFileName='C:\Data\ChipSeqPaper\Figures\InsertSize_EncodeVsPeaks450.tif'
    fig = gcf;
    fig.PaperUnits = 'inches';
    fig.PaperPosition = [0.6 0.6 5.5 10];
    fig.PaperPositionMode = 'manual';
    print(outFileName,'-dtiff','-r650')
    
    fclose all;


