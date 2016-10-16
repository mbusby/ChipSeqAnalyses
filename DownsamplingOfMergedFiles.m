clear
    dirName ='C:\Temp\MonoVsPolyPaper\DownsampledMergedPeaks'
     dirName ='C:\Data\ChipSeqPaper\Data\MergedPeaks\Downsampling'
    
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
    

files= dir(dirName);

for i=3:length(files)
    fileNames{i-2}=files(i).name;
end


for i=1:length(fileNames)
shortNames{i}=strrep(strrep(strrep(strrep(fileNames{i}, '.bam', ''), '.txt', ''), 'out', ''), '_', ' ');
end


for i=1:length(fileNames)
shortNames{i}=shortNames{i}(1:20);
end

shortNamesSorted=shortNames;
    
for i=1:length(fileNames)
    f=fileNames{i};
    stops=strfind(f, '.');
    antibody{i}=f(1:stops(1)-1);
    readDepth(i)=str2num(f(stops(1)+1:stops(2)-1));    
end

antibodies=unique(antibody);
    
%# Column Headers:					1	2	3	4	5	6	7	8
%#PeakID	chr	start	end	strand	Normalized Tag Count	region size	findPeaks Score	Total Tags (normalized to Control Experiment)	Control Tags	Fold Change vs Control	p-value vs Control	Clonal Fold Change
%6-14476	6	27094730	27116323	+	7225.1	21593	24	4627.3	149.5	30.95	0.00E+00	0.95


chromosomes=[];


for i=1:length(fileNames)
     fileName=strcat(dirName,filesep, fileNames{i});
     fid=fopen(strcat(dirName,filesep, fileNames{i}), 'r');    
  
       
     importedData = textscan( fid, '%s %s %f %f %s %f %f %f %f %f %f %f %f', 'Headerlines', 39, 'delimiter', '\t' );
     
     chromosome=importedData{2};
     startPos=importedData{3};
     endPos=importedData{4};
             
     readsInPeaks(i)=sum(endPos-startPos+1);    
     
        
     fclose(fid)
end    



antibodies={  'K27acMono'
    'K27acPoly'
    'K4me1Mono'
    'K4me1Poly'
    'K4me3Mono'
    'K4me3Poly'
}



antibodiesShort={  'H3K27ac'
    'H3K27ac'
    'H3K4me1'
    'H3K4me1'
    'H3K4me3'
    'H3K4me3'
}

goodTop=[51613011/2; 51613011/2; 71172433/2; 71172433/2 ; 36823102/2 ;36823102/2]

ctr=1
fig=figure()
for i=1:2:length(antibodies)
    good=strcmp(antibody, antibodies{i})
    data=sortrows([[0; readDepth(good)'] [0; readsInPeaks(good)']])
   subplot(1, 3,ctr)
    plot(data(:,1), data(:,2), 'r')
    
    hold on
    title(antibodiesShort{i})
    
    good=strcmp(antibody, antibodies{i+1});
    data=sortrows([[0; readDepth(good)'] [0; readsInPeaks(good)']]);
    plot(data(:,1), data(:,2), 'b')
    legend('Mono', 'Poly', 'location', 'SouthEast')
    xlabel('Read Pairs')
    ylabel('Bases Designated as Peaks')
    
    axis([0  goodTop(i)+goodTop(i)*.02 0 max(data(:,2))+(max(data(:,2))*0.05)])
    plotBeef(12)
    ctr=ctr+1
    
    
end




    fig.PaperUnits = 'inches';
    fig.PaperPosition = [0.6 0.6 10.5 2.9];
    fig.PaperPositionMode = 'manual';
    fig.PaperOrientation = 'landscape';
   print('C:\Data\ChipSeqPaper\Figures\saturation\saturationOfActivatingMarks.pdf','-dpdf','-r350')
    




antibodies={ 'K27me3Mono'
    'K27me3Poly'
    'K9me3Mono'
    'K9me3Poly'
}



antibodiesShort={  'H3K27me3'
    'H3K27me3'
    'H3K9me3'
    'H3K9me3'
}

goodTop=[55855959/2; 55855959/2; 56712933/2; 56712933/2 ;]


fig =figure()

ctr=1
for i=1:2:length(antibodies)
    
 
    good=strcmp(antibody, antibodies{i})
    data=sortrows([[0; readDepth(good)'] [0; readsInPeaks(good)']])
   subplot(1, 3,ctr)
    plot(data(:,1), data(:,2), 'r')
    
    hold on
    title(antibodiesShort{i})
    
    good=strcmp(antibody, antibodies{i+1});
    data=sortrows([[0; readDepth(good)'] [0; readsInPeaks(good)']]);
    plot(data(:,1), data(:,2), 'b')
    legend('Mono', 'Poly', 'location', 'NorthEast')
    xlabel('Read Pairs')
    ylabel('Bases Designated as Peaks')
    
    axis([0  goodTop(i)+goodTop(i)*.02 0 max(data(:,2))+(max(data(:,2))*0.05)])
    plotBeef(12)
    ctr=ctr+1
    
    
end



    fig.PaperUnits = 'inches';
    fig.PaperPosition = [0.6 0.6 10.5 2.9];
    fig.PaperPositionMode = 'manual';
    fig.PaperOrientation = 'landscape';
   print('C:\Data\ChipSeqPaper\Figures\saturation\saturationOfRepressiveMarks.pdf','-dpdf','-r350')
    









