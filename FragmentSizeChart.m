clear

dirName='C:\Data\ChipSeqPaper\Data\WCE_ByInsertSize\'
outFileName='C:\Data\ChipSeqPaper\Figures\WCE_ByInsertSize.pdf'


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
    
    tablePerc=0
    
    
    for i=1:length(myBeds)
        if myBeds(i).isdir==0
            s=myBeds(i).name;
            p=strfind(s,'.');

            if isempty(p)==0
                s=s(1:p(1)-1);
            end    

           p=strfind(s,'_');

           if isempty(p)==0
               sortStr(i)=str2num(s(max(p)+1:length(s))); %Assumes antibody is last string
           else
               sortStr(i)=1;
           end
         
            shortNames{i}=strrep(s, '_', ' ');
            
            %pull out the insertSize, unique them, then build a new list.
            %This should work if the insertSize aren't there in the name -
            %it will just sort by the full name then
        else 
            sortStr(i)=0;
            shortNames{i}='NA';
        end
        
    end
    
      insertSize=unique(sort(sortStr(sortStr~=0)));
      
      ctr=1;
      for i=1:length(insertSize)
          ab=insertSize(i);
          for j=1:length(sortStr)
             if myBeds(j).isdir==0 
                 if sortStr(j)==ab
                     fileNamesSorted{ctr}=myBeds(j).name;
                     shortNamesSorted{ctr}=shortNames(j);
                     ctr=ctr+1;
                 end  
             end
          end
          
      end
      
        
    
    
    for i=1:length(fileNamesSorted)
        
         if strcmp(fileNamesSorted{i}(1), '.')==0
               
              i
              fileNamesSorted{i}
             fileName=strcat(dirName,filesep, fileNamesSorted{i})
             [ ~, ~, ~, ~, overlappingBases, ~,  featureType ] = ReadBedtoolsCoverageLowMem( fileName );
             total=sum(overlappingBases)
             
             for j=1:length(annTypes)
                peakTotals(ctrVals,j)=sum(overlappingBases(strcmp(featureType, annTypes{j})));
                tablePerc(ctrVals,j)=sum(overlappingBases(strcmp(featureType, annTypes{j})))/total;
                thisData(i,j)=sum(overlappingBases(strcmp(featureType, annTypes{j})))/total;
             end
              
             clear  overlappingBases;
             clear featureType;

             yTickLabels(ctrVals)=shortNamesSorted{i};
             xTickVals(ctrVals)=insertSize(i);

             ctrVals=ctrVals+1;
         else
            tablePerc(ctrVals,:)=zeros(1, size(tablePerc,2)) ;
                       
            ctrVals=ctrVals+1;
            
             
         end
                 
    end    

    figure
    barh(tablePerc*100, 'stacked')
     colormap(colorscheme)

    set(gca, 'YTick', 1:ctrVals-2)
    set(gca, 'YTickLabel', yTickLabels)
    set(gca, 'TickLength', [0 0]);
    set(gca,'View',[0 -90])
    axis([0 100 0  size(tablePerc,1)])
   
    legend(annTypes, 'location', 'SouthOutside','Orientation','horizontal')    

    set(gca,'fontsize',7);
    title('% Peaks Mapping To ENCODE Genome Regions')
 
figure
subplot(2,1,1)
plot(xTickVals(1:631),  tablePerc(1:631,3)*100, 'w.')
hold on
plot(xTickVals(1:631),  tablePerc(1:631,3)*100, 'r.')
r2= corr(xTickVals(1:631)',  tablePerc(1:631,3))^2
r2text=['R^2=' num2str(r2,2)]
legend(r2text, 'location', 'NorthWest')
legend('boxoff')
ylabel('% TSS (open)')
title('Reads Aligning to Annotated Open and Closed Chromatin')  
subplot(2,1,2)
plot(xTickVals(1:631),  tablePerc(1:631,7)*100, 'w.')
hold on
plot(xTickVals(1:631),  tablePerc(1:631,7)*100, 'b.')
%axis([0 700 0 85])
r2= corr(xTickVals(1:631)',  tablePerc(1:631,7))^2
r2text=['R^2=0.70']
legend(r2text, 'location', 'NorthWest')
legend('boxoff')
ylabel('% Repressed (closed)')
xlabel('Insert Size')  



for i=1:7
    annTypes{i}
    corr(xTickVals(1:631)',  tablePerc(1:631,i))^2
end


mean(tablePerc( (insertSize<=120)', 3))/mean(tablePerc( (insertSize<=700 & insertSize>650 )', 3))


mean(tablePerc( (insertSize<=700 & insertSize>650 )', 7))/mean(tablePerc( (insertSize<=120)', 7))

    fig = gcf;
    fig.PaperUnits = 'inches';
    fig.PaperPosition = [0.6 0.6 5.5 10];
    fig.PaperPositionMode = 'manual';
    print(outFileName,'-dpdf','-r350')
    
    fclose all;