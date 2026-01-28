function produceBCmat3(bcn,prefix,spikeseq,brain)
%to provide globally collapsed refbarcodes for bowtie

%%create global data fist
dataglobal=struct();
dataglobal.counts=dlmread(['thresholds/',prefix,'_global_counts.txt']);
dataglobal.reads=int8(char(textread(['thresholds/',prefix,'_global_seq.txt'],'%s')));
for i=1:size(dataglobal.reads,1)
    isspikes(i)=all(dataglobal.reads(i,25:32)==int8(spikeseq));
end

positionsglobal=[];
positionsglobal.x=dlmread(['thresholds/','bowtieglobal_2u_1.txt']);
positionsglobal.y=dlmread(['thresholds/','bowtieglobal_2u_3.txt']);
clustermatrixglobal.C=sparse(positionsglobal.x,positionsglobal.y,1); %make a sparse matrix using the bowtie columns 1 and 3 as x and y coordinates for nonzero matrix entries

%find connected components
[globalS, globalG]=graphconncomp(clustermatrixglobal.C,'Directed',false); %find the connected graph components


uniq_globalG=unique(globalG,'stable');
uniq_spikes_globalG=unique(globalG(isspikes));
uniq_nonspikes_globalG=uniq_globalG(~ismember(uniq_globalG,uniq_spikes_globalG));
tf_spikes_globalG=ismember(globalG,uniq_spikes_globalG);
spikes_globalG=globalG(tf_spikes_globalG);
tf_nonspikes_globalG=ismember(globalG,uniq_nonspikes_globalG);
nonspikes_globalG=globalG(tf_nonspikes_globalG);


%collapse barcodes to most abundant member of the connected graph component
%in the global
dataglobal_spikes_reads=dataglobal.reads(tf_spikes_globalG,:);
spikes_x=1:length(uniq_spikes_globalG);
[~,loc_spikes_globalG]=ismember(uniq_spikes_globalG(spikes_x),spikes_globalG,'R2012a');
collapsedspikes_globalreads=dataglobal_spikes_reads(loc_spikes_globalG,:);

dataglobal_nonspikes_reads=dataglobal.reads(tf_nonspikes_globalG,:);
nonspikes_x=1:length(uniq_nonspikes_globalG);
[~,loc_nonspikes_globalG]=ismember(uniq_nonspikes_globalG(nonspikes_x),nonspikes_globalG,'R2012a');
collapsednonspikes_globalreads=dataglobal_nonspikes_reads(loc_nonspikes_globalG,:);


%%dealing with each sample
%filter out SSI with no barcodes
filesize=[];
for i=1:length(bcn)
    file=dir(['thresholds/',prefix,'_',int2str(bcn(i)),'_counts.txt']);
    filesize(i)=file.bytes;
end

bcnfilt=bcn(filesize>0);
data=struct();
spikes=struct();

for i=1:length(bcnfilt)
data(i).counts=dlmread(['thresholds/',prefix,'_',int2str(bcnfilt(i)),'_counts.txt']);
data(i).reads=int8(char(textread(['thresholds/',prefix,'_',int2str(bcnfilt(i)),'_seq.txt'],'%s')));
end

%load alignments
positions=[];
graph=[];
spikes=[];
correctedspikes=[];
correctednonspikes=[];
for i=1:length(bcnfilt)
positions(i).x=dlmread(['thresholds/','bowtie',int2str(bcnfilt(i)),'_2u_1.txt']);
positions(i).y=dlmread(['thresholds/','bowtie',int2str(bcnfilt(i)),'_2u_3.txt']);
graph(i).G=globalG(positions(i).y);
tf_spikesG=ismember(graph(i).G,spikes_globalG);
tf_nonspikesG=ismember(graph(i).G,nonspikes_globalG);
graph(i).spikesG=graph(i).G(tf_spikesG);
graph(i).nonspikesG=graph(i).G(tf_nonspikesG);

[uniqspikesG,~,re_spikesG] = unique(graph(i).spikesG); %spikes_globalG remapped from 1...k
collapsed_spikescounts=accumarray(re_spikesG,data(i).counts(tf_spikesG));
collapsed_spikesreads=collapsedspikes_globalreads(ismember(uniq_spikes_globalG,uniqspikesG),:);

[uniqnonspikesG,~,re_nonspikesG] = unique(graph(i).nonspikesG); %spikes_globalG remapped from 1...k
collapsed_nonspikescounts=accumarray(re_nonspikesG,data(i).counts(tf_nonspikesG));
collapsed_nonspikesreads=collapsednonspikes_globalreads(ismember(uniq_nonspikes_globalG,uniqnonspikesG),:);

[correctedspikes(i).counts2u,ixspikes]=sort(collapsed_spikescounts,'descend');
correctedspikes(i).reads2u=collapsed_spikesreads(ixspikes,:);

[correctednonspikes(i).counts2u,ixnonspikes]=sort(collapsed_nonspikescounts,'descend');
correctednonspikes(i).reads2u=collapsed_nonspikesreads(ixnonspikes,:);
end


 
%% remove reads containing homopolymers
minrunlength=7; % as 0.25^7*23=0.0014 or less than 1% of barcodes will have this by chance?
for i=1:length(bcnfilt)
    a_spikes=findhomopolymers(correctedspikes(i).reads2u,minrunlength);
    correctedspikes(i).freads=correctedspikes(i).reads2u(~a_spikes,:);
    correctedspikes(i).fcounts=correctedspikes(i).counts2u(~a_spikes,:);
    a_nonspikes=findhomopolymers(correctednonspikes(i).reads2u,minrunlength);
    correctednonspikes(i).freads=correctednonspikes(i).reads2u(~a_nonspikes,:);
    correctednonspikes(i).fcounts=correctednonspikes(i).counts2u(~a_nonspikes,:);
    data(i).BCseqff=correctednonspikes(i).freads;
    data(i).BCcountsff=correctednonspikes(i).fcounts;
    spikes(i).reads2u=correctedspikes(i).freads;
    spikes(i).counts2u=correctedspikes(i).fcounts;
end

save(['spikes',prefix,'_',brain,'.mat'],'spikes')%save spikes


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


%set thresholds of how many molecule counts a barcode has to have to be
%considered real. by default we leave this at 0
threshold=0; %lowest molecule count considered trustworthy
barcodematrix=barcodematrix(sum(barcodematrix>threshold,2)~=0,:);
refbarcodes=refbarcodes(sum(barcodematrix>threshold,2)~=0,:);

save(['barcodematrix',prefix,'_',brain,'.mat'],'barcodematrix','refbarcodes','-v7.3')


						   