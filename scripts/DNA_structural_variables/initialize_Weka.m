function [] = initialize_Weka()
% installs weka if needed
fname_weka = 'weka.jar';
path_weka = 'weka-3-7-9';
fname_weka_dl = 'weka-3-7-9.zip';
url_weka = 'https://sourceforge.net/projects/weka/files/weka-3-7/3.7.9/weka-3-7-9.zip/download';
if ~isfile(fname_weka)
    disp('%%% Installing Weka %%%')
    websave(fname_weka_dl,url_weka)
    system(['unzip ',fname_weka_dl])
    system(['mv ',path_weka,'/',fname_weka,' weka.jar'])
    system(['rm -rf ',path_weka])
    system(['rm ',fname_weka_dl])
end
end