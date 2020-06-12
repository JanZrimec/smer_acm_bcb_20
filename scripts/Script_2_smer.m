%% Script to generate s-mers
clear all
close all

% generate input sequence permuations
initialize_permn();
[NN3, fname3] = make_seq_perms_fasta(3);
[NN5, fname5] = make_seq_perms_fasta(5);
[NN7, fname7] = make_seq_perms_fasta(7);
[NN9, fname9] = make_seq_perms_fasta(9);

% get structures
path_structures = 'DNA_structural_variables/';
path_script1 = '../';
cd(path_structures)
X3 = get_structures_par([path_script1,fname3], 3);
X5 = get_structures_par([path_script1,fname5], 5);
X7 = get_structures_par([path_script1,fname7], 7);
X9 = get_structures_par([path_script1,fname9], 9);
cd(path_script1)

% get PCAs and save
[Z3,cumvar3,X3_norm,U3,S3,V3] = get_pca(X3);
[Z5,cumvar5,X5_norm,U5,S5,V5] = get_pca(X5);
[Z7,cumvar7,X7_norm,U7,S7,V7] = get_pca(X7);
[Z9,cumvar9,X9_norm,U9,S9,V9] = get_pca(X9);

save([fname3,'.mat'],'Z3','cumvar3','X3_norm','U3','S3','V3','-v7.3')
save([fname5,'.mat'],'Z5','cumvar5','X5_norm','U5','S5','V5','-v7.3')
save([fname7,'.mat'],'Z7','cumvar7','X7_norm','U7','S7','V7','-v7.3')
save([fname9,'.mat'],'Z9','cumvar9','X9_norm','U9','S9','V9','-v7.3')

% functions
function [] = initialize_permn()
% Installs permn if needed from saved .zip file
fname_parfor = 'permn.m';
if ~isfile(fname_parfor)
   disp('%%% Installing permn tool %%%')
   system('unzip Additional_tools/permn.zip')
end
end

function [NNx, fname] = make_seq_perms_fasta(n)
% returns sequence permutations and saves fasta
v = [1,2,3,4];
NNx = permn(v,n);
for i=1:4^n
    NNxd(i).Sequence = int2nt(NNx(i,:));
    NNxd(i).Header = sprintf('%d',i);
end 
fname = sprintf('NN%d.fasta',n);
fastawrite(fname,NNxd) 
end

function [Z,cumvar,X_norm,U,S,V] = get_pca(X)
% returns PCA variables via SVD function
[m,~] = size(X); %m examples, n features
Xmn = repmat(mean(X),m,1);
Xsd = repmat(std(X),m,1);
X_norm = (X - Xmn)./(Xsd);
Sigma = (1/m)*X_norm'*X_norm;
[U,S,V] = svd(Sigma);
%Ureduce = U(1,1:k);
Z = X_norm*U;
cumvar = cumsum(var(Z))/sum(var(Z));
end