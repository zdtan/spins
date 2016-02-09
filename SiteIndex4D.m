function [ site, i, j, k ] = SiteIndex4D( n, K, L, M )
%SiteIndex4D Returns site index as (site, i, j, k).

if (n<1 || n>8*K*L*M)
    disp(['Error in index n, n = ',num2str(n)]);
end

k = floor((n-1)/(8*K*L)) + 1;
n = mod(n-1,8*K*L) + 1;

j = floor((n-1)/(8*K)) + 1;
n = mod(n-1,8*K) + 1;

i = floor((n-1)/8) + 1;
site = mod(n-1,8) + 1;

end