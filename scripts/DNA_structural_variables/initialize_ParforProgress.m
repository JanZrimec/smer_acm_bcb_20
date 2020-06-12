function [] = initialize_ParforProgress()
% Installs ParforProgress if needed
fname_parfor = 'ParforProgress.m';
if ~isfile(fname_parfor)
   disp('%%% Installing ParforProgress tool %%%')
   system('git clone https://github.com/dgolden1/ParforProgress.git')
   system('mv ParforProgress/*.m .')
   system('rm -rf ParforProgress')
end
end