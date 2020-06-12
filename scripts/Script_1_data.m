%% Processes and saves data

%% DNA structural variables
% Data and scripts for computation of structural variables were
% obtained from https://github.com/JanZrimec/DNA_structural_variables

% url{1} = 'https://github.com/JanZrimec/DNA_structural_variables/blob/master/NN_structural_properties.mat';
% url{2} = 'https://github.com/JanZrimec/DNA_structural_variables/blob/master/NNN_structural_properties.mat';
% url{3} = 'https://github.com/JanZrimec/DNA_structural_variables/blob/master/List_structural_variables.csv';
% for i =1:length(url)
%     tmp = split(url{i},'/');
%     urlwrite(url{i},tmp{end})
% end

% Install test tool DNA_structural_variables
cd('..')
system('git clone https://github.com/JanZrimec/DNA_structural_variables.git')
cd('Script_0_data')

%% OriT Positive set
% Data was obtained from Zrimec & Lapanje 2018
% Stored as a Fasta file containing the 64 sequences
% at link http://dnatools.eu/MOB/data_orit_64.fasta

clear all

url{1} = 'https://github.com/JanZrimec/Plasmid_MOB_prediction_oriT/blob/master/data_orit_204.fasta';
url{2} = 'https://github.com/JanZrimec/Plasmid_MOB_prediction_oriT/blob/master/data_orit_64.fasta';
url{3} = 'https://github.com/JanZrimec/Plasmid_MOB_prediction_oriT/blob/master/orit_64_metadata.mat';
for i =1:length(url)
    tmp = split(url{i},'/');
    urlwrite(url{i},tmp{end})
end


fasta64 = fastaread('data_orit_64.fasta');

p64_pos = cell([64,1]);
p64_pos_y = cell([64,1]);
for i=1:64
   p64_pos{i,1} = fasta64(i).Sequence; 
   tmp = split(fasta64(i).Header,'_');
   p64_pos_y{i,1} = tmp{2};
end

%% OriT Negative set
% negative examples - 2 sets - mixed coding noncoding in random from GC

% Load data
% For creation of the negative sets a metadata file from Zrimec & Lapanje 
% 2018 is used, which contains whole plasmid sequences and cutouts for
% obtaining negatives
load('orit_64_metadata.mat')

% reset random number generator
rng default

% nic is in the middle pos 1000/1001
for i=1:64
    if isempty(MET{i,16})
        p64_neg{i,1} = MET{i,11};  
        p64_base(i) = basecount(MET{i,7}.Sequence);
    else
        p64_neg{i,1} = MET{i,16};
        p64_base(i) = basecount(MET{i,13}.Sequence);
    end
end

% generate genomic negatives from surrounding +/- 200 to 800 bp around nic 
ind_rand=rand(64,2);
for i=1:64
    if ~strcmp(p64_neg{i,1},p64_neg{6,1})
        itmp(i,1) = round(ind_rand(i,1)*600)+200; % region 200-800 from nic is set to center
        itmp(i,2) = round(ind_rand(i,2)*600)+200;
        p64_neg{i,2} = p64_neg{i,1}(1000-itmp(i,1)-139:1000-itmp(i,1)+90);
        p64_neg{i,3} = p64_neg{i,1}(1000+itmp(i,2)-139:1000+itmp(i,2)+90);
    end
end

% generate random sequences according to nucleotide frequencies - 4 per
% element - skip elements with too small size
for i=1:64
    if p64_base(i).A+p64_base(i).C+p64_base(i).G+p64_base(i).T > 2000
        p64_rand{i,1} = randseq(230,'FromStructure',p64_base(i));
        p64_rand{i,2} = randseq(230,'FromStructure',p64_base(i));
        p64_rand{i,3} = randseq(230,'FromStructure',p64_base(i));
        p64_rand{i,4} = randseq(230,'FromStructure',p64_base(i));
    end
end

% filter negatives and make final sets
% delete missing elements
i = 1;
while i <= length(p64_neg)
    if isempty(p64_neg{i,2})
        p64_neg(i,:) = [];
    end
    i = i+1;
end
i=1;
while i <= length(p64_rand)
    if isempty(p64_rand{i,1})
        p64_rand(i,:) = [];
    end
    i = i+1;
end

% negative genomic: take randomly either first or second sequences and then 
% fill missing by taking additional 5 from the rest randomly
ind_rand2 = round(rand(59,1));
ind_rand3 = round(rand(5,1)*59);

for i=1:59
    p64_neg_sel{i} = p64_neg{i,2+ind_rand2(i,1)};
    obratni{i} = p64_neg{i,2+1-ind_rand2(i,1)};
end
for i=1:5
    p64_neg_sel{59+i}=obratni{ind_rand3(i,1)};
end
p64_neg_sel = p64_neg_sel';

% negative random: take first and random from second
ind_rand4 = round(rand(6,1)*59);

p64_rand_sel = p64_rand(:,1);
for i=1:6
    p64_rand_sel{58+i}=p64_rand{ind_rand4(i,1),2};
end

%% OriT analysis of distances between sequences in the datasets
for i=1:64
    for j=1:64
        D_pos(i,j) = seqpdist([p64_pos{i};p64_pos{j}], 'Method', 'p-distance');
        D_neg(i,j) = seqpdist([p64_neg_sel{i};p64_neg_sel{j}], 'Method', 'p-distance');
        D_rand(i,j) = seqpdist([p64_rand_sel{i};p64_rand_sel{j}], 'Method', 'p-distance');
    end
end

% positive sequences have two close clusters - these are:
% 1 bp difference was observed between RSF1010 (p31) and R1162 (p32) and an 
% 11 bp difference was observed between pSymA (p38) and pSymB (p39), next
% elements were over 33bp apart 
figure(1)
imagesc(D_pos)
colorbar
title('p-distance positive')

p64_pos_pdist_sort = sort(nonzeros(triu(D_pos)));
% how many bps apart
disp(p64_pos_pdist_sort(1:10)*230)

for i=1:2
    [X,Y] = find(D_pos==p64_pos_pdist_sort(i));
    disp([X,Y])
    remove(i) = X(1);
end

% all negative sequences are far appart
figure(2)
imagesc(D_neg)
colorbar
title('p-distance genomic negatives')

figure(3)
imagesc(D_rand)  
colorbar
title('p-distance random negatives')

disp(min(nonzeros(triu(D_neg))))
disp(min(nonzeros(triu(D_rand))))

%% oriT processing and save
% merge into one dataset and remove two elements
disp(remove)
p186 = [p64_pos;p64_neg_sel;p64_rand_sel];
p62_pos_y = p64_pos_y;
p186(remove)=[];
p186(remove+64)=[];
p186(remove+64*2)=[];
p62_pos_y(remove)=[];

% create grouping variables for positive and whole datasets
chars = unique(p62_pos_y);
p62_y = zeros([62,1]);
for i =1:4
	p62_y(ismember(p62_pos_y,chars{i})) = i-1;
end

p186_y = [zeros([62,1]);ones([62,1]);2*ones([62,1])];

% distance matrix
for i=1:186
    for j=1:186
        D_p186(i,j) = seqpdist([p186{i};p186{j}], 'Method', 'p-distance');
    end
end
figure(4)
imagesc(D_p186)  
colorbar
title('p-distance p186')

disp(min(nonzeros(triu(D_p186))))

% save datasets p64_pos and p186 to hdf5
save('Dataset_orit.mat','p186','p186_y','p62_y','-v7.3')

%% Promoter regions - Positive negative
% Positive and negative sets from Gusmao et al. 2011 via personal
% correspondence

clear all
%cd('/Users/zrimec/Box Sync/Projects-active/Clank4A')
rng default

% import Gusmao 100bp sequences x 812 x 3
Gusmaocolipos = read_files_Gusmao('Gusmao_coli_pos.csv');
Gusmaocolimix1 = read_files_Gusmao('Gusmao_coli_mix1.csv');
Gusmaocolictrl = read_files_Gusmao('Gusmao_coli_ctrl.csv');

for i=1:812
    gusmaodata{i,1} = upper(Gusmaocolipos{i});
    gusmaodata{i+812,1} = upper(Gusmaocolimix1{i});
    gusmaodata{i+812*2,1} = upper(Gusmaocolictrl{i});
end

% dilute data to 200 per sample
r1 = randsample(812,200);
r2 = randsample(812,200);
r3 = randsample(812,200);

for i=1:200
    gus600{i,1} = gusmaodata{r1(i)};
    gus600{i+200,1} = gusmaodata{r2(i)+812};
    gus600{i+400,1} = gusmaodata{r3(i)+1624};
end

% create grouping variable
gus600_y = [zeros(200,1);ones(200,1);2*ones(200,1)];

%% Promoter regions - Sigma groups

% data was obtained from RegulonDB (Release: 9.1 Date: 04-07-2016) 
RegDB_all = read_files_Regulon('RegulonDB-PromoterSet.txt');
rng default

[Sigma_names,~,J] = unique(RegDB_all(:,5));
occ = histc(J, 1:numel(Sigma_names));
C = arrayfun(@num2str,occ,'UniformOutput',false);
disp([Sigma_names,C])

% choose 6 sigmas
RegDB_sigma = cell(0,8);
for i = 3:8
    RegDB_sigma = [RegDB_sigma;RegDB_all(ismember(RegDB_all(:,5),Sigma_names{i}),:)];
end

% remove rows with empty sequence fields
RegDB_sigma(ismember(RegDB_sigma(:,6),''),:) = [];

% dilute all groups to size of smallest group
[Sigma_names,~,J] = unique(RegDB_sigma(:,5));
occ = histc(J, 1:numel(Sigma_names));
C = arrayfun(@num2str,occ,'UniformOutput',false);
disp([Sigma_names,C])

RegDB_sigma_600 = cell(0,8);
for i=1:6
    dilute_index = randsample(occ(i),min(occ));
    sliced = RegDB_sigma(ismember(RegDB_sigma(:,5),Sigma_names{i}),:);
    RegDB_sigma_600 = [RegDB_sigma_600;sliced(dilute_index,:)];
end
    
sig600 = cellfun(@upper,RegDB_sigma_600(:,6),'UniformOutput',false);

% create grouping variable
sig600_y = [];
for i = 1:6
    sig600_y = [sig600_y; (i-1)*ones(min(occ),1)];
end

%% Promoter regions distance analysis and save

% distance matrix
for i=1:length(gus600)
    for j=i+1:length(gus600)
        D_gus600(i,j) = seqpdist([gus600{i};gus600{j}],'Method','p-distance');
    end
end

figure(5)
imagesc(D_gus600)  
colorbar
title('p-distance promoters pos/neg dataset')

for i=1:length(sig600)
    for j=i+1:length(sig600)
        D_sig600(i,j) = seqpdist([sig600{i};sig600{j}],'Method','p-distance');
    end
end

figure(6)
imagesc(D_sig600)  
colorbar
title('p-distance promoters Sigma dataset')

disp(min(nonzeros(triu(D_gus600))))
disp(min(nonzeros(triu(D_sig600))))

% remove close elements from sigma set
sig600_pdist_sort = sort(nonzeros(triu(D_sig600)));
% how many bps apart
disp(sig600_pdist_sort(1:10)*length(sig600{1}))

[X,Y] = find(D_sig600==sig600_pdist_sort(1));
disp([X,Y])
tmp = sort(unique(X));
remove = tmp(1:2);

sig600(remove)=[];
sig600_y(remove)=[];

% recheck distances sigma
clear D_sig600
for i=1:length(sig600)
    for j=i+1:length(sig600)
        D_sig600(i,j) = seqpdist([sig600{i};sig600{j}],'Method','p-distance');
    end
end

disp(min(nonzeros(triu(D_sig600))))

% save datasets p64_pos and p186 to hdf5
save('Dataset_prom.mat','gus600','gus600_y','sig600','sig600_y','-v7.3')

%% functions
function tlines = read_files_Gusmao(fname)

fid=fopen(fname);
tline = fgetl(fid);
tlines = cell(0,1);
while ischar(tline)
    tlines{end+1,1} = upper(erase(tline,','));
    tline = fgetl(fid);
end
fclose(fid);
end

function tlines = read_files_Regulon(fname)

fid=fopen(fname);
tline = fgetl(fid);
tlines = cell(0,8);
while ischar(tline)
    tmp = strsplit(tline,'\t','CollapseDelimiters', false);
    tlines(end+1,1:length(tmp)) = tmp;
    tline = fgetl(fid);
end
fclose(fid);

end
