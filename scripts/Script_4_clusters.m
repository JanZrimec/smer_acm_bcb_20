%% clustering analysis
% description: 
% ka=[4,8,16,32,64,128,256] -> 2-8 bits of information + raw data

%% make clusters kmeans

% set clustering parameters and open parallel
opts=statset('Display','final','MaxIter',1000,'UseParallel',true);
parpool

% load data
load('../data/NN3.fasta.mat')
load('../data/NN5.fasta.mat')
load('../data/NN7.fasta.mat')
load('../data/NN9.fasta.mat')
cd '../scripts'

% run kmeans
ka=[4,8,16,32];
for i=1:4
    [idx3{i},C3{i},sumd3{i},D3{i}] = ...
    kmeans(Z3,ka(i),'Options',opts,'Replicates',10); 
    [a3{i},b3{i}]=hist(idx3{i},unique(idx3{i}));
end
C3{5}=Z;
idx3{5}=(1:64)';

ka=[4,8,16,32,64,128,256,1024];
for i=1:7
    [idx5{i},C5{i},sumd5{i},D5{i}] = ...
    kmeans(Z5,ka(i),'Options',opts,'Replicates',10);
    [a5{i},b5{i}]=hist(idx5{i},unique(idx5{i}));
end
C5{8}=Z;
idx5{8}=(1:1024)';

ka=[4,8,16,32,64,128,256,16384];
for i=1:7
    [idx7{i},C7{i},sumd7{i},D7{i}] = ...
    kmeans(Z7,ka(i),'Options',opts,'Replicates',10);
    [a7{i},b7{i}]=hist(idx7{i},unique(idx7{i}));
end
C7{8}=Z;
idx7{8}=(1:16384)';

ka=[4,8,16,32,64,128,256,262144];
for i=1:7
    [idx9{i},C9{i},sumd9{i},D9{i}] = ...
    kmeans(Z9,ka(i),'Options',opts,'Replicates',10); 
    [a9{i},b9{i}]=hist(idx9{i},unique(idx9{i}));
end
C9{8}=Z;
idx9{8}=(1:262144)';    

% save clustering results
save('clusters_bits_9_11_16.mat',idx3,C3,sumd3,D3,a3,b3,...
                                 idx5,C5,sumd5,D5,a5,b5,...
                                 idx7,C7,sumd7,D7,a7,b7,...
                                 idx9,C9,sumd9,D9,a9,b9,'-v7.3')

                             
%% cluster analysis
% elbow
% silhouette

load('../data/clusters_bits_9_11_16.mat')
load('../data/NN3.fasta.mat')
load('../data/NN5.fasta.mat')
load('../data/NN7.fasta.mat')
load('../data/NN9.fasta.mat')
cd 'NP_MANOVA_bootstrap'

% run 
for i = 1:5
    [F3(i), Sw3(i), St3(i)] = Anova_F_ratio(Z3,idx3{i});
end

for i = 1:8
    disp(i)
    [F5(i), Sw5(i), St5(i)] = Anova_F_ratio(Z5,idx5{i});
    [F7(i), Sw7(i), St7(i)] = Anova_F_ratio(Z7,idx7{i});
%     [F9A(i), Sw9A(i), St9A(i)] = Anova_F_ratio(Z9A,idx9{i});
%     [F9(i), Sw9(i), St9(i)] = Anova_F_ratio(Z9,idx9{i});
end

% Z9 contrain to computable options - random sampling (10%)
rng(113)
rind = randi(length(Z9),[round(length(Z9)/10),1]); %10-fold dilution
for i = 1:8
    disp(i)
    [F9(i), Sw9(i), St9(i)] = Anova_F_ratio(Z9(rind,:),idx9{i}(rind));
end

%% silhouette
% same z9 filtering

figure(20) % dump plots here

for i = 1:4
    [sil3{i},h3{i}] = silhouette(Z3,idx3{i});
    sil3_m(i) = mean(sil3{i});
end
for i =1:7
    disp(i)
    [sil5{i},h5{i}]=silhouette(Z5,idx5{i});
    [sil7{i},h7{i}]=silhouette(Z7,idx7{i});
    [sil9{i},h9{i}]=silhouette(Z9(rind),idx9{i}(rind));
    sil5_m(i) = mean(sil5{i});
    sil7_m(i) = mean(sil7{i});
    sil9_m(i) = mean(sil9{i});
end


%% plot
% §
% https://stats.stackexchange.com/questions/11175/elbow-criteria-to-determine-number-of-cluster

close all
figure(11)
hold on
subplot(2,4,1)
semilogx([4,8,16,32],(St3(1:4)-Sw3(1:4))./St3(1:4),'r')
title('Elbow method')
xlabel('Cluster size k')
ylabel('Variance explained')
subplot(2,4,5)
semilogx([4,8,16,32],sil3_m(1:4),'r')
title('Silhouette method')
xlabel('Cluster size k')
ylabel('Mean Silhouette value')

subplot(2,4,2)
semilogx([4,8,16,32,64,128,256],(St5(1:7)-Sw5(1:7))./St5(1:7),'r')
title('Elbow method')
xlabel('Cluster size k')
ylabel('Variance explained')
subplot(2,4,6)
semilogx([4,8,16,32,64,128,256],sil5_m,'r')
title('Silhouette method')
xlabel('Cluster size k')
ylabel('Mean Silhouette value')

subplot(2,4,3)
semilogx([4,8,16,32,64,128,256],(St7(1:7)-Sw7(1:7))./St7(1:7),'r')
title('Elbow method')
xlabel('Clusters (bits)')
ylabel('Variance explained')
subplot(2,4,7)
semilogx([4,8,16,32,64,128,256],sil7_m,'r')
title('Silhouette method')
xlabel('Cluster size k')
ylabel('Mean Silhouette value')

subplot(2,4,4)
semilogx([4,8,16,32,64,128,256],(St9(1:7)-Sw9(1:7))./St9(1:7),'r')
title('Elbow method')
xlxlabel('Cluster size k')
ylabel('Variance explained')
subplot(2,4,8)
semilogx([4,8,16,32,64,128,256],sil9_m,'r')
title('Silhouette method')
xlabel('Cluster size k')
ylabel('Mean Silhouette value')

% export data to table
table_cl = zeros(8,7);
table_cl(1,1:4) = (St3(1:4)-Sw3(1:4))./St3(1:4);
table_cl(2,1:7) = (St5(1:7)-Sw5(1:7))./St5(1:7);
table_cl(3,1:7) = (St7(1:7)-Sw7(1:7))./St7(1:7);
table_cl(4,1:7) = (St9(1:7)-Sw9(1:7))./St9(1:7);
table_cl(5,1:4) = sil3_m(1:4);
table_cl(6,1:7) = sil5_m;
table_cl(7,1:7) = sil7_m;
table_cl(8,1:7) = sil9_m;


%% determine optimal k
% evalcluster

klists = {[4,8,16,32],[4,8,16,32,64,128,256],...
            [4,8,16,32,64,128,256],[4,8,16,32,64,128,256]};
        
methods = {'CalinskiHarabasz','DaviesBouldin','gap','silhouette'}; % ,'gap','silhouette'
      
myfunc = @(X,K)(kmeans(X, K, 'emptyaction','singleton',...
    'replicate',3));

results = cell(length(fnames),length(methods));
optimal_k = zeros(length(fnames),length(methods));

for ii=1:length(methods)
    disp(methods(ii))
    results{1,ii} = evalclusters(Z3,myfunc,methods(ii),'klist',klists{1});
    results{2,ii} = evalclusters(Z5,myfunc,methods(ii),'klist',klists{2});
    results{3,ii} = evalclusters(Z7,myfunc,methods(ii),'klist',klists{3});
    results{4,ii} = evalclusters(Z9(rind),myfunc,methods(ii),'klist',klists{4});
    for i=1:4
        optimal_k(i,ii) = results{i,ii}.OptimalK;
    end
end
