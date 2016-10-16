

 blu=[0.0588    0.4745    0.7490]; %  TSS 
     
gre=[ 1.0000    0.3       0.3; ]; %  E

monoNames={'M1';
    'M2';
    'M3';
    'M4';
    'M5';}

polyNames={'P1';
    'P2';
    'P3';
    'P4';
    'P5';}

    
%{
         % Downsampled normalized by insert size
K27ac	29320	mAb-K27ac-1.dupsRemoved.bam.normalized.calls		
K27ac	28413	mAb-K27ac-2.dupsRemoved.bam.normalized.calls		
K27ac	33784	pAb-K27ac-1.dupsRemoved.bam.normalized.calls		
K27ac	33575	pAb-K27ac-2.dupsRemoved.bam.normalized.calls		
K27me3	10229	1_Exp4-K562-2A-K27me3-mono_H9V71ADXX.1.aligned.duplicates_marked.dupsRemoved.bam.normalized.calls		
K27me3	10148	1_Exp4-K562-2B-K27me3-mono_H9V71ADXX.1.aligned.duplicates_marked.dupsRemoved.bam.normalized.calls		
K27me3	10002	1_Exp4-K562-2C-K27me3-mono_H9V71ADXX.1.aligned.duplicates_marked.dupsRemoved.bam.normalized.calls		
K27me3	10196	1_Exp4-K562-2D-K27me3-mono_H9V71ADXX.1.aligned.duplicates_marked.dupsRemoved.bam.normalized.calls		
K27me3	13984	1_Exp4-K562-2E-K27me3-poly_H9V71ADXX.1.aligned.duplicates_marked.dupsRemoved.bam.normalized.calls		
K27me3	12937	1_Exp4-K562-2F-K27me3-poly_H9V71ADXX.1.aligned.duplicates_marked.dupsRemoved.bam.normalized.calls		
K27me3	18721	1_Exp4-K562-2G-K27me3-poly_H9V71ADXX.1.aligned.duplicates_marked.dupsRemoved.bam.normalized.calls		
K4me1	52118	K562_K4me1_Ido.dupsRemoved.bam.normalized.calls		
K4me1	62227	K562_K4me1_LizExp11.dupsRemoved.bam.normalized.calls		
K4me1	54468	mAb-K4me1-1.dupsRemoved.bam.normalized.calls		
K4me1	59427	pAb-K4me1-1.dupsRemoved.bam.normalized.calls		
K4me1	58730	pAb-K4me1-2.dupsRemoved.bam.normalized.calls		
k4me3	15002	1_Exp4-K562-1A-K4me3-mono_H9V71ADXX.1.aligned.duplicates_marked.dupsRemoved.bam.normalized.calls		
k4me3	15167	1_Exp4-K562-1B-K4me3-mono_H9V71ADXX.1.aligned.duplicates_marked.dupsRemoved.bam.normalized.calls		
k4me3	15185	1_Exp4-K562-1C-K4me3-mono_H9V71ADXX.1.aligned.duplicates_marked.dupsRemoved.bam.normalized.calls		
k4me3	15053	1_Exp4-K562-1D-K4me3-mono_H9V71ADXX.1.aligned.duplicates_marked.dupsRemoved.bam.normalized.calls		
k4me3	15904	1_Exp4-K562-1E-K4me3-poly_H9V71ADXX.1.aligned.duplicates_marked.dupsRemoved.bam.normalized.calls		
k4me3	15982	1_Exp4-K562-1F-K4me3-poly_H9V71ADXX.1.aligned.duplicates_marked.dupsRemoved.bam.normalized.calls		
k4me3	15956	1_Exp4-K562-1G-K4me3-poly_H9V71ADXX.1.aligned.duplicates_marked.dupsRemoved.bam.normalized.calls		
k4me3	16160	1_Exp4-K562-1H-K4me3-poly_H9V71ADXX.1.aligned.duplicates_marked.dupsRemoved.bam.normalized.calls		
K9me3	6083	K562_K9me3_Ido.dupsRemoved.bam.normalized.calls		
K9me3	5326	K562_K9me3_LizExp11.dupsRemoved.bam.normalized.calls		
K9me3	5533	mAb-K9me3-2.dupsRemoved.bam.normalized.calls		
K9me3	9399	pAb-K9me3-1.dupsRemoved.bam.normalized.calls		
K9me3	7411	pAb-K9me3-2.dupsRemoved.bam.normalized.calls		
%}
    
mono=[29320 28413];
poly=[33784 33575];
mMono=mean(mono);
mPoly=mean(poly);
(mPoly-mMono)/mPoly;

[h pK27ac]=ttest2(mono, poly)

figure
subplot(3,2,1)
bar(1:length(mono), mono, 'FaceColor', blu)
hold on
bar(length(mono)+1:(length(poly)+length(mono)), poly, 'FaceColor', gre)
tickLabels={monoNames{1:length(mono)}  polyNames{1:length(poly)} }
set(gca, 'xTick', 1:(length(poly)+length(mono)) )
set(gca, 'xTickLabel', tickLabels)
set(gca, 'XTickLabelRotation', 0)
ylabel('Peaks Detected')
title('H3K27ac')



mono=[10229 10148 10002 10196 ];
poly= [13984 12937 18721];
mMono=mean(mono);
mPoly=mean(poly);
(mPoly-mMono)/mPoly;
[h pK27me3]=ttest2([10229 10148 10002 10196 ], [13984 12937 18721])

subplot(3,2,2)
bar(1:length(mono), mono, 'FaceColor', blu)
hold on
bar(length(mono)+1:(length(poly)+length(mono)), poly, 'FaceColor', gre)
tickLabels={monoNames{1:length(mono)}  polyNames{1:length(poly)} }
set(gca, 'xTick', 1:(length(poly)+length(mono)) )
set(gca, 'xTickLabel', tickLabels)
set(gca, 'XTickLabelRotation', 0)
ylabel('Peaks Detected')
title('H3K27me3')



%No difference

mono=[52118 62227 54468	];
poly= [ 59427 58730];
mMono=mean(mono);
mPoly=mean(poly);
(mPoly-mMono)/mPoly;

[h pK4me1]=ttest2(mono, poly)	

subplot(3,2,3)
bar(1:length(mono), mono, 'FaceColor', blu)
hold on
bar(length(mono)+1:(length(poly)+length(mono)), poly, 'FaceColor', gre)
tickLabels={monoNames{1:length(mono)}  polyNames{1:length(poly)} }
set(gca, 'xTick', 1:(length(poly)+length(mono)) )
set(gca, 'xTickLabel', tickLabels)
set(gca, 'XTickLabelRotation', 0)
ylabel('Peaks Detected')
title('H3K4me1')




mono=[15002 15167	15185 15053	];
poly= [ 15904 15982 15956 16160];
mMono=mean(mono);
mPoly=mean(poly);
(mPoly-mMono)/mPoly;

[h pK4me3]=ttest2([15002 15167	15185 15053], [ 15904 15982 15956 16160])
%Poly

subplot(3,2,4)
bar(1:length(mono), mono, 'FaceColor', blu)
hold on
bar(length(mono)+1:(length(poly)+length(mono)), poly, 'FaceColor', gre)
tickLabels={monoNames{1:length(mono)}  polyNames{1:length(poly)} }
set(gca, 'xTick', 1:(length(poly)+length(mono)) )
set(gca, 'xTickLabel', tickLabels)
set(gca, 'XTickLabelRotation', 0)
ylabel('Peaks Detected')
title('H3K4me3')

mono=[6083 5326  5533]
poly=[ 9399 7411]
subplot(3,2,5)
bar(1:length(mono), mono, 'FaceColor', blu)
hold on
bar(length(mono)+1:(length(poly)+length(mono)), poly, 'FaceColor', gre)
tickLabels={monoNames{1:length(mono)}  polyNames{1:length(poly)} }
set(gca, 'xTick', 1:(length(poly)+length(mono)) )
set(gca, 'xTickLabel', tickLabels)
ylabel('Peaks Detected')
set(gca, 'XTickLabelRotation', 0)
title('H3K9me3')


[h pK9me3]=ttest2([6083 5326  5533], [ 9399 7411])
%poly

peaks=[29320
28413
33784
33575
10229
10148
10002
10196
13984
12937
18721
52118
62227
54468
59427
58730
15002
15167
15185
15053
15904
15982
15956
16160
6083
5326
5533
9399
7411
]





mAb-K27ac-1.dupsRemoved.bam.normalized			0.5589
mAb-K27ac-2.dupsRemoved.bam.normalized			0.5149
pAb-K27ac-1.dupsRemoved.bam.normalized			0.5496
pAb-K27ac-2.dupsRemoved.bam.normalized			0.5442
1_Exp4-K562-2A-K27me3-mono_H9V71ADXX.1.aligned.duplicates_marked.dupsRemoved.bam.normalized			0.0112
1_Exp4-K562-2B-K27me3-mono_H9V71ADXX.1.aligned.duplicates_marked.dupsRemoved.bam.normalized			0.0112
1_Exp4-K562-2C-K27me3-mono_H9V71ADXX.1.aligned.duplicates_marked.dupsRemoved.bam.normalized			0.0112
1_Exp4-K562-2D-K27me3-mono_H9V71ADXX.1.aligned.duplicates_marked.dupsRemoved.bam.normalized			0.0113
1_Exp4-K562-2E-K27me3-poly_H9V71ADXX.1.aligned.duplicates_marked.dupsRemoved.bam.normalized			0.016
1_Exp4-K562-2F-K27me3-poly_H9V71ADXX.1.aligned.duplicates_marked.dupsRemoved.bam.normalized			0.0147
1_Exp4-K562-2G-K27me3-poly_H9V71ADXX.1.aligned.duplicates_marked.dupsRemoved.bam.normalized			0.0231
K562_K4me1_Ido.dupsRemoved.bam.normalized			0.175
K562_K4me1_LizExp11.dupsRemoved.bam.normalized			0.258
mAb-K4me1-1.dupsRemoved.bam.normalized			0.1902
pAb-K4me1-1.dupsRemoved.bam.normalized			0.2512
pAb-K4me1-2.dupsRemoved.bam.normalized			0.247
1_Exp4-K562-1A-K4me3-mono_H9V71ADXX.1.aligned.duplicates_marked.dupsRemoved.bam.normalized			0.239
1_Exp4-K562-1B-K4me3-mono_H9V71ADXX.1.aligned.duplicates_marked.dupsRemoved.bam.normalized			0.2689
1_Exp4-K562-1C-K4me3-mono_H9V71ADXX.1.aligned.duplicates_marked.dupsRemoved.bam.normalized			0.267
1_Exp4-K562-1D-K4me3-mono_H9V71ADXX.1.aligned.duplicates_marked.dupsRemoved.bam.normalized			0.2564
1_Exp4-K562-1E-K4me3-poly_H9V71ADXX.1.aligned.duplicates_marked.dupsRemoved.bam.normalized			0.3301
1_Exp4-K562-1F-K4me3-poly_H9V71ADXX.1.aligned.duplicates_marked.dupsRemoved.bam.normalized			0.3189
1_Exp4-K562-1G-K4me3-poly_H9V71ADXX.1.aligned.duplicates_marked.dupsRemoved.bam.normalized			0.3107
1_Exp4-K562-1H-K4me3-poly_H9V71ADXX.1.aligned.duplicates_marked.dupsRemoved.bam.normalized			0.2992
K562_K9me3_Ido.dupsRemoved.bam.normalized			0.0539
K562_K9me3_LizExp11.dupsRemoved.bam.normalized			0.0503
mAb-K9me3-2.dupsRemoved.bam.normalized			0.0411
pAb-K9me3-1.dupsRemoved.bam.normalized			0.0504
pAb-K9me3-2.dupsRemoved.bam.normalized			0.0458





mono=[0.5589 0.5149];
poly=[0.5496 0.5442];
mMono=mean(mono);
mPoly=mean(poly);
(mPoly-mMono)/mPoly;

[h pK27ac]=ttest2(mono, poly)

figure
subplot(3,2,1)
bar(1:length(mono), 100*mono, 'FaceColor', blu)
hold on
bar(length(mono)+1:(length(poly)+length(mono)), 100*poly, 'FaceColor', gre)
tickLabels={monoNames{1:length(mono)}  polyNames{1:length(poly)} }
set(gca, 'xTick', 1:(length(poly)+length(mono)) )
set(gca, 'xTickLabel', tickLabels)
ylabel('%Reads in Peaks')
set(gca, 'XTickLabelRotation', 0)
title('H327ac')


mono=[0.0112 0.0112 0.0112 0.0113];
poly=[0.016 0.0147 0.0231];
mMono=mean(mono);
mPoly=mean(poly);
(mPoly-mMono)/mPoly;

[h pK27me3]=ttest2(mono, poly)

subplot(3,2,2)
bar(1:length(mono), 100*mono, 'FaceColor', blu)
hold on
bar(length(mono)+1:(length(poly)+length(mono)), 100*poly, 'FaceColor', gre)
tickLabels={monoNames{1:length(mono)}  polyNames{1:length(poly)} }
set(gca, 'xTick', 1:(length(poly)+length(mono)) )
set(gca, 'xTickLabel', tickLabels)
ylabel('%Reads in Peaks')
set(gca, 'XTickLabelRotation', 0)
title('H3K27me3')



mono=[0.175 0.258 0.1902];
poly=[0.2512 0.247];
mMono=mean(mono);
mPoly=mean(poly);
(mPoly-mMono)/mPoly;

[h pK4me1]=ttest2(mono, poly)

subplot(3,2,3)
bar(1:length(mono), 100*mono, 'FaceColor', blu)
hold on
bar(length(mono)+1:(length(poly)+length(mono)), 100*poly, 'FaceColor', gre)
tickLabels={monoNames{1:length(mono)}  polyNames{1:length(poly)} }
set(gca, 'xTick', 1:(length(poly)+length(mono)) )
set(gca, 'xTickLabel', tickLabels)
ylabel('%Reads in Peaks')
set(gca, 'XTickLabelRotation', 0)
title('H3K4me1')


mono=[0.239 0.2689 0.267 0.2564];
poly=[0.3301 0.3189 0.3107 0.2992];
mMono=mean(mono);
mPoly=mean(poly);
(mPoly-mMono)/mPoly;

[h pk4me3]=ttest2(mono, poly)
subplot(3,2,4)
bar(1:length(mono), 100*mono, 'FaceColor', blu)
hold on
bar(length(mono)+1:(length(poly)+length(mono)), 100*poly, 'FaceColor', gre)
tickLabels={monoNames{1:length(mono)}  polyNames{1:length(poly)} }
set(gca, 'xTick', 1:(length(poly)+length(mono)) )
set(gca, 'xTickLabel', tickLabels)
ylabel('%Reads in Peaks')
set(gca, 'XTickLabelRotation', 0)
title('H3K4me3')



mono=[0.0539 0.0503 0.0411];
poly=[0.0504 0.0458];
mMono=mean(mono);
mPoly=mean(poly);
(mPoly-mMono)/mPoly;

[h pK9me3]=ttest2(mono, poly)

subplot(3,2,5)
bar(1:length(mono), 100*mono, 'FaceColor', blu)
hold on
bar(length(mono)+1:(length(poly)+length(mono)), 100*poly, 'FaceColor', gre)
tickLabels={monoNames{1:length(mono)}  polyNames{1:length(poly)} }
set(gca, 'xTick', 1:(length(poly)+length(mono)) )
set(gca, 'xTickLabel', tickLabels)
ylabel('%Reads in Peaks')
set(gca, 'XTickLabelRotation', 0)
title('H3K9me3')

[pK27ac; pK27me3; pK4me1; pK4me3; pK9me3]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55

%{
'H3K27ac	Mono	Rep1'	0.875827462
'H3K27ac	Mono	Rep2'	0.88654666
'H3K27ac	Poly	Rep	0.803149702
'H3K27ac	Poly	Rep	0.802304752
[]			
'H3K27me3	Mono	Rep	0.982650746
'H3K27me3	Mono	Rep	0.982597444
'H3K27me3	Mono	Rep	0.98383682
'H3K27me3	Mono	Rep	0.982263075
'H3K27me3	Poly	Rep	0.983874525
'H3K27me3	Poly	Rep	0.983876945
'H3K27me3	Poly	Rep	0.984716163
[]			
'H3K4me1	Mono	Rep	0.302589707
'H3K4me1	Mono	Rep	0.294970667
'H3K4me1	Mono	Rep	0.335833689
'H3K4me1	Poly	Rep	0.320974501
'H3K4me1	Poly	Rep	0.325164315
[]			
'H3K4me3	Mono	Rep	0.930006784
'H3K4me3	Mono	Rep	0.926957566
'H3K4me3	Mono	Rep	0.926552589
'H3K4me3	Mono	Rep	0.92796546
'H3K4me3	Poly	Rep	0.915612392
'H3K4me3	Poly	Rep	0.914946881
'H3K4me3	Poly	Rep	0.915670143
'H3K4me3	Poly	Rep	0.915049168
[]			
'H3K9me3	Mono	Rep	0.740185475
'H3K9me3	Mono	Rep	0.720775397
'H3K9me3	Mono	Rep	0.630380819
'H3K9me3	Poly	Rep	0.690185041
'H3K9me3	Poly	Rep	0.675887342
%}
0.875827462
0.88654666
0.803149702
0.802304752

0.982650746
0.982597444
0.98383682
0.982263075
0.983874525
0.983876945
0.984716163

0.302589707
0.294970667
0.335833689
0.320974501
0.325164315

0.930006784
0.926957566
0.926552589
0.92796546
0.915612392
0.914946881
0.915670143
0.915049168

0.740185475
0.720775397
0.630380819
0.690185041
0.675887342




mono=[0.875827462 0.88654666];
poly=[0.803149702 0.802304752];
mMono=mean(mono);
mPoly=mean(poly);
(mPoly-mMono)/mPoly;

[h pK27ac]=ttest2(mono, poly)

figure
subplot(3,2,1)
bar(1:length(mono), 100*mono, 'FaceColor', blu)
hold on
bar(length(mono)+1:(length(poly)+length(mono)), 100*poly, 'FaceColor', gre)
tickLabels={monoNames{1:length(mono)}  polyNames{1:length(poly)} }
set(gca, 'xTick', 1:(length(poly)+length(mono)) )
set(gca, 'xTickLabel', tickLabels)
ylabel('%Reads in Peaks')
set(gca, 'XTickLabelRotation', 0)
title('H327ac')


mono=[0.982650746 0.982597444 0.98383682 0.982263075 ];
poly=[0.983874525 0.983876945 0.984716163];
mMono=mean(mono);
mPoly=mean(poly);
(mPoly-mMono)/mPoly;

[h pK27me3]=ttest2(mono, poly)

subplot(3,2,2)
bar(1:length(mono), 100*mono, 'FaceColor', blu)
hold on
bar(length(mono)+1:(length(poly)+length(mono)), 100*poly, 'FaceColor', gre)
tickLabels={monoNames{1:length(mono)}  polyNames{1:length(poly)} }
set(gca, 'xTick', 1:(length(poly)+length(mono)) )
set(gca, 'xTickLabel', tickLabels)
ylabel('%Reads in Peaks')
set(gca, 'XTickLabelRotation', 0)
title('H3K27me3')



mono=[0.302589707 0.294970667 0.335833689 ];
poly=[0.320974501 0.325164315];
mMono=mean(mono);
mPoly=mean(poly);
(mPoly-mMono)/mPoly;

[h pK4me1]=ttest2(mono, poly)

subplot(3,2,3)
bar(1:length(mono), 100*mono, 'FaceColor', blu)
hold on
bar(length(mono)+1:(length(poly)+length(mono)), 100*poly, 'FaceColor', gre)
tickLabels={monoNames{1:length(mono)}  polyNames{1:length(poly)} }
set(gca, 'xTick', 1:(length(poly)+length(mono)) )
set(gca, 'xTickLabel', tickLabels)
ylabel('%Reads in Peaks')
set(gca, 'XTickLabelRotation', 0)
axis([0 6 0 100])
title('H3K4me1')


mono=[0.930006784 0.926957566 0.926552589 0.92796546];
poly=[0.915612392 0.914946881 0.915670143 0.915049168];
mMono=mean(mono);
mPoly=mean(poly);
(mPoly-mMono)/mPoly;

[h pk4me3]=ttest2(mono, poly)
subplot(3,2,4)
bar(1:length(mono), 100*mono, 'FaceColor', blu)
hold on
bar(length(mono)+1:(length(poly)+length(mono)), 100*poly, 'FaceColor', gre)
tickLabels={monoNames{1:length(mono)}  polyNames{1:length(poly)} }
set(gca, 'xTick', 1:(length(poly)+length(mono)) )
set(gca, 'xTickLabel', tickLabels)
ylabel('%Reads in Peaks')
set(gca, 'XTickLabelRotation', 0)
title('H3K4me3')



mono=[0.740185475 0.720775397 0.630380819 ];
poly=[0.690185041 0.675887342];
mMono=mean(mono);
mPoly=mean(poly);
(mPoly-mMono)/mPoly;

[h pK9me3]=ttest2(mono, poly)

subplot(3,2,5)
bar(1:length(mono), 100*mono, 'FaceColor', blu)
hold on
bar(length(mono)+1:(length(poly)+length(mono)), 100*poly, 'FaceColor', gre)
tickLabels={monoNames{1:length(mono)}  polyNames{1:length(poly)} }
set(gca, 'xTick', 1:(length(poly)+length(mono)) )
set(gca, 'xTickLabel', tickLabels)
ylabel('%Reads in Peaks')
set(gca, 'XTickLabelRotation', 0)
title('H3K9me3')

[pK27ac; pK27me3; pK4me1; pK4me3; pK9me3]



