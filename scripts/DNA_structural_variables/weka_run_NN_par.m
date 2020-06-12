function [out,name2]=weka_run_NN_par(klas,fname,nbr,W)

name = sprintf('count_W%d_sum_T1_w0%d_rand.arff',W,nbr);
zero = zeros(size(klas));
para = [W,nbr,3,1]; %W(width), wide(neighboring), 3, T(treshold)
name2 = PBDNN_make_arff_euk_nerand(fname,klas,zero,para(1),para(2),para(3),para(4));    
tmp = [fname,'.test'];
command = ['java -cp weka.jar weka.classifiers.trees.M5P -t ',name,' -T ',name2,' -classifications "weka.classifiers.evaluation.output.prediction.CSV -file ',tmp,'"'];
[~] = unix(command);
tmp2 = dlmread(tmp,',',1,0);
out = tmp2(:,3);

end