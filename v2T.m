function T = v2T(v,n)
%convert column vector to a Topelitz matrix, using offset indices

v = v(:)';
L = length(v);
T = zeros(n+L-1,n);
for col = 1 : n
    T((1:L)+col-1,col) = v;
end

end

