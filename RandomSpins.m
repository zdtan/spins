function [ spinMat ] = RandomSpin( N )
%RandomSpin - Generate N random spins of unit length
%       using solid angle method.

rvals = 2*rand(N,1)-1;
elevation = asin(rvals);
azimuth = 2*pi*rand(N,1);
[x,y,z] = sph2cart(azimuth,elevation,1);	% unit vectors
spinMat = horzcat(x,y,z);

end
