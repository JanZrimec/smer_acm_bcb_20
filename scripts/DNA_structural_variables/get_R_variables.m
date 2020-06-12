function out = get_R_variables(fname)
% necessary R packages are installed within r script
fname_r_vars = 'R_variables.csv';
path_to_rscript = '/Library/Frameworks/R.framework/Resources';
system(['cp ',fname,' R_tmp.fa'])
system([path_to_rscript,'/Rscript get_R_variables.R'])
system('rm R_tmp.fa*')

% read data
data_r = importdata(fname_r_vars);
tmp = char(data_r.colheaders);
data_r_colnames = unique(cellstr(tmp(:,2:4)));
data_r_slice = cell(length(data_r_colnames),1);
out = cell(size(data_r.data,1),4);
for i = 1:length(data_r_colnames)
    idx = strcmp(cellstr(tmp(:,2:4)),data_r_colnames{i});
    data_r_slice = data_r.data(:,idx);
    out(:,i) = mat2cell(data_r_slice,ones(1,size(data_r_slice,1)));
end
system(['rm ',fname_r_vars])
end