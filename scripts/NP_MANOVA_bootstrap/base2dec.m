function out = base2dec(in)
%support function

[n,m] = size(in);
out = zeros(n,1);
for i = 1:m
    out = out + in(:,i)*power(4,m-i);
end

end