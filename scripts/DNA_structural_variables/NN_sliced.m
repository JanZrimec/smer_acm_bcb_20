function NNx2 = NN_sliced(in,win,weka,nbr,W,melt,ii,x,NNsp,NNNsp)
% sliced function for parallel implementation NN_structures_par
    
% load('NN_structural_properties.mat')
% load('NNN_structural_properties.mat')
i = 1;
NNx{i,1} = in;
NNx{i,2} = nt2int(in);
NNx2 = cell(i,61);

for j=1:53 % size NN structural properties

    for k=1:x-1  % size of nn
        % transform nn int to matrix NNsp - 4x(n1-1)+n2
        NNx{i,2+j}(k) = NNsp{2+j,(6 + (4*(NNx{i,2}(k)-1)+NNx{i,2}(k+1)))};
    end

    %sum if specified in matrix
    if strcmp(NNsp{2+j,5},'Y')
        for kk=1:x-win+1
            NNx2{i,2+j}(1,kk) = sum(NNx{i,2+j}(kk:kk+win-2));
        end
    elseif strcmp(NNsp{2+j,5},'M')
        for kk=1:x-win+1
            NNx2{i,2+j}(1,kk) = mean(NNx{i,2+j}(kk:kk+win-2));
        end
    elseif strcmp(NNsp{2+j,5},'I')
        for kk=1:x-win+1
            NNx2{i,2+j}(1,kk) = win/sum(1./NNx{i,2+j}(kk:kk+win-2)); %width/sum(1./a(i:i+(width-1)));
        end
    else                                                % elseif strcmp(NNsp{2+j,5},'N')
        for kk=1:x-win+1
            NNx2{i,2+j}(1,kk) = mean(NNx{i,2+j}(kk:kk+win-2));  % NNx{i,2+j}(x) = [];  % matrix index is out of range for deletion
        end
    end

end

for j=1:4 %size NNN structural properties - shift by NN!

    for k=1:x-2 % size of nnn
        l=1;
        while (strcmp(NNNsp{l,1}(1:3),NNx{i,1}(k:k+2)) || strcmp(NNNsp{l,1}(5:7),NNx{i,1}(k:k+2))) == 0 %loop breaks when NNN is found
            l=l+1;
        end
        NNx{i,2+53+j}(k) = NNNsp{l,1+j};
    end
    for kk=1:x-win+1
        NNx2{i,2+53+j}(1,kk) = mean(NNx{i,2+53+j}(kk:kk+win-3)); %mean dnaze data
    end

end

if weka==true
    fname= sprintf('tmp_id%d',ii);
    [temptidd,fname2] = weka_run_NN_par(NNx{i,1},fname,nbr,W); % first position
    NNx2{i,60} = temptidd';
    [~] = unix(['rm ',fname,'.test ',fname2]);
end

if melt==true
    temp = oligoprop(NNx{i,1});
    NNx2{i,61} = temp.Tm(1);
end

end