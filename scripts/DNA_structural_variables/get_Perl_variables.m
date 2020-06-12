function out = get_Perl_variables(fname)
% the necessary modules are installed if needed, except BioPerl
url_perl = 'http://dna.bu.edu/orchid/orchid2_server_predictions.tar.gz';
fname_perl = 'orchid2_server_predictions.tar.gz';
fname_perl2 = 'print_orchid2_pattern_from_fasta.pl';
path_perl = 'orchid2_server_predictions';
fname_perl_vars = 'Perl_variables.txt';
if ~isfile(fname_perl2)
    disp('%%% Installing Perl modules, make sure BioPerl is installed %%%')
    websave(fname_perl,url_perl)
    system(['tar xopf ',fname_perl])
    system(['mv ',path_perl,'/',fname_perl2,' ',fname_perl2])
    system(['mv ',path_perl,'/PredictCleavagePatterns2.pm PredictCleavagePatterns2.pm'])
    system(['rm -rf ',path_perl])
    system(['rm ',fname_perl])
end

% run
system(['cp ',fname,' perl_tmp.fa'])
system(['perl print_orchid2_pattern_from_fasta.pl -f perl_tmp.fa > ',fname_perl_vars])
system('rm perl_tmp.fa')

% read data
data_perl = importdata(fname_perl_vars);
out = data_perl.data;
system(['rm ',fname_perl_vars])
end