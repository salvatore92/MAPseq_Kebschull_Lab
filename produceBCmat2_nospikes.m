function produceBCmat2_nospikes(bcn,prefix,brain)
%process MAPseq data to produce barcode matrix

%cd thresholds;
%bcn=212:229;
%bcn=1:length(bcn);

%filter out SSI with no barcodes
filesize=[];
for i=1:length(bcn)
    file=dir(['pre/',prefix,'_pre',int2str(bcn(i)),'_counts.txt']);
    filesize(i)=file.bytes;
end

bcnfilt=bcn(filesize>0);
data=struct();
for i=1:length(bcnfilt)
data(i).counts=dlmread(['pre/',prefix,'_pre',int2str(bcnfilt(i)),'_counts.txt']);
data(i).reads=int8(char(textread(['pre/',prefix,'_pre',int2str(bcnfilt(i)),'_seq.txt'],'%s')));
end
 
%save(['data_',brain,'.mat','data','-v7.3');
 
 
 
 
%% Finish error correction by reading in bowtie alignments and finding connected graph components
 
%load alignments
positions1=[];
for i=1:length(bcnfilt)
positions1(i).x=dlmread(['pre/','bowtiepre',int2str(bcnfilt(i)),'_2u_1.txt']);
positions1(i).y=dlmread(['pre/','bowtiepre',int2str(bcnfilt(i)),'_2u_3.txt']);
clustermatrix1(i).C=sparse(positions1(i).x,positions1(i).y,1); %make a sparse matrix using the bowtie columns 1 and 3 as x and y coordinates for nonzero matrix entries
end
 %save('clustermatrix1.mat','clustermatrix1','-v7.3');
% 
 %load clustermatrix1
 
%find connected components
graph=[];
for i=1:length(bcnfilt)
i
    [graph(i).S,graph(i).G]=graphconncomp(clustermatrix1(i).C,'Directed',false); %find the connected graph components
end
% save('graph1.mat','graph');
% 
% load graph1
 
%collapse barcodes to most abundant member of the connected graph component
 
for i=1:length(bcnfilt)
x=1:graph(i).S;
[tf,loc]=ismember(x,graph(i).G,'R2012a');
collapsedreads=data(i).reads(loc,:);
collapsedcounts=accumarray(graph(i).G',data(i).counts);%'
[corrected(i).counts2u,ix]=sort(collapsedcounts,'descend');
corrected(i).reads2u=collapsedreads(ix,:);
data(i).BCseqff=corrected(i).reads2u;
data(i).BCcountsff=corrected(i).counts2u;
end
%save('data.mat','data','-v7.3');
 
 
 
 %% load data.m and match the barcodes from the different areas against each other to produce a barcode matrix
% this is code modified to in situ MAPseq, so it will run without an
% injection site.


%% find overlap of barcodes
%% set up reference set of barcodes
%load data
%load data.mat

%collect all suquences detected in the target sites
refbarcodes_tmp=[];
for i=1:length(data)
refbarcodes_tmp=[refbarcodes_tmp;data(i).BCseqff];  
end
refbarcodes=unique(refbarcodes_tmp,'rows'); %unique barocodes to get reference set



%% construct barcodematrix by matching barcodes in target sites to reference barcodes

barcodematrix=zeros(size(refbarcodes,1),length(data));%initiate the barcode matrix

for i=1:length(data)
    %pull out reads and counts into new variables for ease of use
    BCreads=data(i).BCseqff;
    BCmolcounts=data(i).BCcountsff;
    [ind,loc]=ismember(BCreads,refbarcodes,'rows');
    barcodematrix(loc(loc~=0),i)=BCmolcounts(ind);   
end

prebarcodematrix=barcodematrix;
prerefbarcodes=refbarcodes;


count=sum(prebarcodematrix,2);

outFile = sprintf('%s_BCglobal_quickout.txt', prefix);

    % 1) Sort by count (descending)
    [count_sorted, idxref] = sort(count, 'descend');
    ref_sorted = refbarcodes(idxref, :);

    % 2) Convert int8 ASCII to char
    seqChar = char(ref_sorted);

    % 3) Write output file
   out = compose("%d\t%s", count_sorted, string(seqChar));
    writelines(out, outFile);

for i=1:length(data)
    [~,idx]=ismember(data(i).BCseqff,ref_sorted,'rows');
    graph(i).localtopreglobal=idx(graph(i).G);
    filename1 = sprintf('bowtie%d_2u_1.txt', bcn(i));
    localidx = (1:length(graph(i).localtopreglobal))';
    writematrix(localidx, filename1);
    filename3 = sprintf('bowtie%d_2u_3.txt', bcn(i));
    writematrix(graph(i).localtopreglobal, filename3);
end


save(['prebarcodematrix',prefix,'_',brain,'.mat'],'prebarcodematrix','prerefbarcodes','ref_sorted','-v7.3')

						   