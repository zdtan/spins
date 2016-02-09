function [ T, Asites, energy, totalM ] = ...
	AutocorrelationFunction( T, spinSol, spinMat, Jmat, alpha, sigma, K, L, M, AFM )
%AutocorrelationFunction Calculates autocorrelation function at each time-step.
%   The switch 'AFM' is true (1) or false (0) for
%     antiferromagnetic or ferromagnetic Jmat respectively.
%	spinSol is a 3D array.

%% system size
N = 8*K*L*M; % 8 sites per unit cell

%% number of time steps
nTimeSteps = size(T,1);	% should be equal to Tmesh

%% sanity check on size of spins
if (size(spinSol,1) ~= N) || (size(spinSol,2) ~= 3) || (size(spinSol,3) ~= nTimeSteps)
    disp('== ERROR: Inconsistency in size of spin solution! ==');
    return
end
if size(spinMat,1) ~= N
    disp('== ERROR: Inconsistency in size of spin matrix! ==');
    return
end

%% autocorrelation functions
Afull = zeros(nTimeSteps,N);	% for each site
Asites = zeros(nTimeSteps,8);	% averaged over unit cells
%~ Akag = zeros(nTimeSteps,1);	% averaged over all kagome sites
%~ Atri = zeros(nTimeSteps,1);	% averaged over all triangular sites

for t = 1:nTimeSteps
	Afull(t,:) = dot(spinMat,spinSol(:,:,t),2);	%need to transpose?
end
% break up based on site
Asites1 = Afull(:,1:8:N);
Asites2 = Afull(:,2:8:N);
Asites3 = Afull(:,3:8:N);
Asites4 = Afull(:,4:8:N);
Asites5 = Afull(:,5:8:N);
Asites6 = Afull(:,6:8:N);
Asites7 = Afull(:,7:8:N);
Asites8 = Afull(:,8:8:N);
% average over sites
Asites = horzcat(sum(Asites1,2),sum(Asites2,2), ...
	sum(Asites3,2),sum(Asites4,2),sum(Asites5,2), ...
	sum(Asites6,2),sum(Asites7,2),sum(Asites8,2))/(N/8);

%~ Akag = sum( horzcat(Asites(:,1),Asites(:,2),Asites(:,3), ...
	%~ Asites(:,5),Asites(:,6),Asites(:,7)), 2 )/6;
%~ Atri = sum( horzcat(Asites(:,4),Asites(:,8)), 2 )/2;

%% also check total energy conserved
energy = zeros(nTimeSteps,1);	% without total magnetisation term
%Jmat = AdjacencyMatrix( alpha, K, L, M, AFM );	% adjacency matrix
for t = 1:nTimeSteps
	LFmat = -Jmat*spinSol(:,:,t);
	for n = 1:N
		energy(t) = -0.5*sum(dot(spinSol(:,:,t),LFmat,2));
	end
	energy(t) = energy(t)/N - (-3/4 - (alpha*sigma)^2/2);
	% per spin, without total magnetisation term
end

%% also check total spin conserved
totalM = zeros(nTimeSteps,3);
for t = 1:nTimeSteps
	totalM(t,:) = sum(spinSol(:,:,t),1);
end

end
