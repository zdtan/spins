function [ spinMat ] = GenerateRS( sigma, K, L, M );
%GenerateRS - Generate random spin configuration.
%   sigma gives length of triangular spins (kagome spins of unit length).

%% system size
N = 8*K*L*M; % 8 sites per unit cell

%% create random spins
spinMat = RandomSpins(N);
% scale length of triangular lattice spins
spinMat(4:4:N,:) = sigma*spinMat(4:4:N,:);

end