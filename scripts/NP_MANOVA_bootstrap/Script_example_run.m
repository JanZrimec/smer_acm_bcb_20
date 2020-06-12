% data
in{1} = 'ACGTCGGCAGCGTCGCTCGAACCGTGCCTGCCGAGGCACG';
in{2} = 'ACGTCGGCAGCGTCGCTCGTTCCGTGCCTGCCGAGGCACG';
in{3} = 'ACGTCGGCAGCGTCGCTCGATCCGTGCCTGCCGAGGCACG';
in{4} = 'GAGTGTTTTTGGGGCACCCCCCTGTGTCCTCCCCGGCTCT';
in{5} = 'GAGTGTTTTTGGGGCACCCGGCTGTGTCCTCCCCGGCTCT';
in{6} = 'GAGTGTTTTTGGGGCACCCGCCTGTGTCCTCCCCGGCTCT';

% bootstraps
num_seqs = length(in);
num_boots = 10;
seq_len = length(in{1});
rng(111)
boots = make_bootstraps(in,num_seqs,num_boots,seq_len);

% anova
[F, Sw, St] = Anova_F_ratio_pairwise(char(in),[0,0,0,1,1,1]);
F_boot = zeros([num_boots,2,2]);
Sw_boot = zeros([num_boots,1]);
St_boot = zeros([num_boots,1]);
for i=1:num_boots
    [F_boot(i,:,:),Sw_boot(i),St_boot(i)] = Anova_F_ratio_pairwise(char(boots{i}),[0,0,1,1]);
end

% norm_to_boots
[F_norm,F_boot_norm] = norm_to_boots(F(1,2),F_boot(:,1,2));

% get_R2
R2 = get_R2(Sw,St);

% pvalues_boots
p = pvalues_boots(F(1,2),F_boot(:,1,2),num_boots);
