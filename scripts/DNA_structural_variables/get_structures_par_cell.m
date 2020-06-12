function [out] = get_structures_par_cell(fname,window)
% function that calculates structural variables from input sequences 
% fname ... input fasta file
% window ... size of sliding window for predictions (3 to 20 bps)
% out ... output sturctures in cell array format 

% process input
seqs = struct2cell(fastaread(fname));
seqs(1,:) = [];
seqs = seqs';
tmps = length(seqs);
x = length(seqs{1});
initialize_ParforProgress();
pp = ParforProgress; %https://se.mathworks.com/matlabcentral/fileexchange/48705-parforprogress-class

% assert statements and define additional parameters
assert(x>2,"Input must be at least 3-mers.")
assert(window>2&window<=20,"Allowed window sizes 3 to 20.")
assert(window<=x,"Window cannot be larger than k.")
if x < 5
Weka = false; R = false; Perl = false; melt = false;
elseif x < 7
Weka = true; R = true; Perl = false; melt = false;
else
Weka = true; R = true; Perl = true; melt = true;
end

% initialize Weka
initialize_Weka()
neighbor = 6; % best model
if x < 30
    neighbor = 0;
end
W = floor(window/5)*5;

% run Matlab and Weka
d1 = load('NN_structural_properties.mat');
d2 = load('NNN_structural_properties.mat');
out_Mat = cell(tmps,1);
parfor i = 1:tmps
    tmp = NN_sliced(seqs{i},window,Weka,neighbor,W,melt,i,x,d1.NNsp,d2.NNNsp);
    tmp
    out_Mat{i} = tmp(1,3:end);
    iteration_number = step(pp, i);
    disp(fprintf('Finished iteration %d of %d\n',iteration_number,i)); 
end

% run R and Perl
% {'Hel';'MGW';'Pro';'Rol'}
out_R = cell(tmps,4); out_Perl = cell(tmps,1);
if R
    out_R = get_R_variables(fname);
end
if Perl
    out_Perl = get_Perl_variables(fname);
end

% process results
out = cell(tmps,64);
for i=1:tmps
    out(i,:) = [out_Mat{i}(1:57),out_R(i,:),out_Perl(i,:),out_Mat{i}(58:59)];
end

end