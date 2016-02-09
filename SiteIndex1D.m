function [ siteIndex ] = SiteIndex1D( site, i, j, k, K, L, M )
%SiteIndex1D Returns site index as number between 1 and 8*K*L*M inclusive.

if (site<1 || site>8)
    disp(['Error in index site, site = ',num2str(site)]);
end
if (i<1 || i>K)
    disp(['Error in index i, i = ',num2str(i)]);
end
if (j<1 || j>L)
    disp(['Error in index j, j = ',num2str(j)]);
end
if (k<1 || k>M)
    disp(['Error in index k, k = ',num2str(k)]);
end

siteIndex = site + 8*(i-1) + 8*K*(j-1) + 8*K*L*(k-1);

end

