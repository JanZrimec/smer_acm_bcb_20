function boots = make_bootstraps(in,num_seqs,num_boots,seq_len)
% make bootstraps as specified on matlab.com
% num_seqs ... number of sequences
% num_boots ... number of bootstrap replicates
% seq_len ... length of seqeunces

boots = cell(num_boots,1);
for n = 1:num_boots
    reorder_index = randsample(seq_len,seq_len,true);
    for i = num_seqs:-1:1 %reverse order to preallocate memory;
        bootseq{i} = strrep(in{i}(reorder_index),'-','');
    end
    boots{n} = bootseq;
end

end