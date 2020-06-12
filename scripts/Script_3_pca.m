%% PCA analysis
% performs analysis of the PCAs and plotting

clear all
%cd '/Users/zrimec/Box Sync/Projects-active/Clank4A/Script_2_pca'

%% load results from script_1_smer
% Zx are PCAs, Ux are loadings, cumvarX are cumulative variances
load('../data/NN3.fasta.mat','Z3','U3','cumvar3')
load('../data/NN5.fasta.mat','Z5','U5','cumvar5')
load('../data/NN7.fasta.mat','Z7','U7','cumvar7')
load('../data/NN9.fasta.mat','Z9','U9','cumvar9')
load('DNA_structural_variables/struct_data_16.mat');

%% make table with properties
table_smer = zeros(4,4);
table_smer(1:4,1) = [3,5,7,9]; % smer
table_smer(1:4,2) = [4^3,4^5,4^7,4^9]; % permuations
table_smer(1:4,3) = [length(U3),length(U5),length(U7),length(U9)]; % permuations
tresh = 0.99;
table_smer(1,4) = find(cumvar3>tresh,1);
table_smer(2,4) = find(cumvar5>tresh,1);
table_smer(3,4) = find(cumvar7>tresh,1);
table_smer(4,4) = find(cumvar9>tresh,1);

%% tables of variance per pca
% for each smer size, mean, std
% separate
maxvar = max(table_smer(:,4));
table_var = zeros(6,maxvar);
table_var(1,:) = [cumvar3(1,1),diff([cumvar3(1:maxvar-1);[cumvar3(2:maxvar)]])];
table_var(2,:) = [cumvar5(1,1),diff([cumvar5(1:maxvar-1);[cumvar5(2:maxvar)]])];
table_var(3,:) = [cumvar7(1,1),diff([cumvar7(1:maxvar-1);[cumvar7(2:maxvar)]])];
table_var(4,:) = [cumvar9(1,1),diff([cumvar9(1:maxvar-1);[cumvar9(2:maxvar)]])];
table_var(5,:) = mean(table_var(1:4,:));
table_var(6,:) = std(table_var(1:4,:));

% cumulative
table_cumvar = zeros(6,maxvar);
table_cumvar(1,:) = cumvar3(1,1:maxvar);
table_cumvar(2,:) = cumvar5(1,1:maxvar);
table_cumvar(3,:) = cumvar7(1,1:maxvar);
table_cumvar(4,:) = cumvar9(1,1:maxvar);
table_cumvar(5,:) = mean(table_cumvar(1:4,:));
table_cumvar(6,:) = std(table_cumvar(1:4,:));

% on average how many variables needed to pass different tresholds
treshs = [0.6,0.8,0.9,0.99];
treshs(2,:) =zeros(1,4);
for i=1:4
    treshs(2,i) = find(table_cumvar(5,:)>treshs(1,i),1);
end

%% plot of pcas for each smer size
% shows clusters of structures are formed - support for clustering
close all
figure(1)
subplot(1,4,1)
biplot(U3(:,1:3),'Score',Z3(:,1:3))
title('3-s-mer')
subplot(1,4,2)
biplot(U5(:,1:3),'Score',Z5(:,1:3))
title('5-s-mer')
% for higher use random subsampling to 1000 points
rng(111)
rnd7 = randsample(length(Z7),1000);
subplot(1,4,3)
biplot(U7(:,1:3),'Score',Z7(rnd7,1:3))
title('7-s-mer')
rnd9 = randsample(length(Z9),1000);
subplot(1,4,4)
biplot(U9(:,1:3),'Score',Z9(rnd9,1:3))
title('9-s-mer')

%% most important variables from pca loadings

% number of variables to display
numvar = 10; 
% cutoff to process the loadings
cutoff = 0.2; 
% cut loadings below cutoff for better analysis
U3_cutoff=(U3>cutoff)|(U3<-cutoff);
U5_cutoff=(U5>cutoff)|(U5<-cutoff);
U7_cutoff=(U7>cutoff)|(U7<-cutoff);
U9_cutoff=(U9>cutoff)|(U9<-cutoff);

% sum loadings acorss different smer sizes to increase power of results
Usum = [U3_cutoff(:,1:numvar);zeros(7,10)] + ...
       [U5_cutoff(:,1:numvar);zeros(2,10)] + ...
       U7_cutoff(:,1:numvar) + U9_cutoff(:,1:numvar);

% plot
Usum_normed = Usum./max(Usum);
figure(20)
imagesc(Usum_normed)
colorbar
xlabel('PC')
ylabel('Structural variable')
title('Relative importance of structural variables')

%% Tables of most important structural variables per pca
% import list of structures
load('DNA_structural_variables/struct_data_16.mat');
 
% select using either cutoff 0.75 or 1 the groups of properties per 10 PCAs
select1 = Usum_normed >= 0.75;
select2 = Usum_normed == 1;
struct_select1 = cell(1,9);
struct_select2 = cell(1,9);
for i=1:10
    struct_select1 = [struct_select1;[repmat({i},(sum(select1(:,i))),1),struct(select1(:,i),:)]];
    struct_select2 = [struct_select2;[repmat({i},(sum(select2(:,i))),1),struct(select2(:,i),:)]];
end

%% most important structural variables overall
% for each smer size take optimal number of pcas and retrieve equal number
% of variables
[I_sort_U3sq(:,2),I_sort_U3sq(:,1)]=sort(sum(power(U3(:,1:14),2),2),'descend');
[I_sort_U5sq(:,2),I_sort_U5sq(:,1)]=sort(sum(power(U5(:,1:17),2),2),'descend');
[I_sort_U7sq(:,2),I_sort_U7sq(:,1)]=sort(sum(power(U7(:,1:17),2),2),'descend');
[I_sort_U9sq(:,2),I_sort_U9sq(:,1)]=sort(sum(power(U9(:,1:18),2),2),'descend');

struct_sort = [I_sort_U3sq(1:20,:),I_sort_U5sq(1:20,:),I_sort_U7sq(1:20,:),I_sort_U9sq(1:20,:)];

% resort based on all 4 - make like an index into which is accumulated,
% then sort
scount = zeros(64,1);
for i = 1:57
    scount(I_sort_U3sq(i,1)) = scount(I_sort_U3sq(i,1))+I_sort_U3sq(i,2);
    scount(I_sort_U5sq(i,1)) = scount(I_sort_U3sq(i,1))+I_sort_U3sq(i,2);
    scount(I_sort_U7sq(i,1)) = scount(I_sort_U3sq(i,1))+I_sort_U3sq(i,2);
    scount(I_sort_U9sq(i,1)) = scount(I_sort_U3sq(i,1))+I_sort_U3sq(i,2);
end
scount = scount/4;

[scount_sort(:,2),scount_sort(:,1)] = sort(scount,'descend');

% lists of top 20 variables per smer size and across all smers
struct_sort = [struct_sort,scount_sort(1:20,:)];

% table of top 20 variables
scount_vars = struct(scount_sort(1:20,1),:);
