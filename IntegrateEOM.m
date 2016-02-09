function [ T, spinSol ] = IntegrateEOM( spinMat, Jmat, ...
	alpha, K, L, M, AFM, Tmax, Tmesh )
%IntegrateEOM Integrates equations of motion for initial configuration spinMat.
%   The switch 'AFM' is true (1) or false (0) for
%     antiferromagnetic or ferromagnetic Jmat respectively.

%% system size
N = 8*K*L*M; % 8 sites per unit cell

%% sanity check on size of spins
if size(spinMat,1) ~= N
    disp('== ERROR: Inconsistency in size of spin matrix! ==');
    return
end

%%%%% reshape into column vectors %%%%%
%% this converts into column vector
% first component, then second component, then 3rd.
%reshape(mat,3*N,1);
%% this converts back into matrix
%reshape(col,N,3);

%% define odefun 'eom'
function dy = eom(t,y)
	dy = zeros(3*N,1);    % a column vector
	%~ h = -Jmat*reshape(y,N,3);
	%~ dy = reshape(cross(reshape(y,N,3),reshape(h,N,3)),...
			%~ 3*N,1);	% a column vector
	dy = reshape(cross(reshape(y,N,3),-Jmat*reshape(y,N,3)),...
			3*N,1);
end

%% call ode45 solver
options = odeset('RelTol',1e-15);
tspan = linspace(0,Tmax,Tmesh);	% keep Tmesh=1000 points consistently
[T,Y] = ode45(@eom,tspan,reshape(spinMat,3*N,1),options);

%% number of time steps
nTimeSteps = size(T,1);
spinSol = zeros(N,3,nTimeSteps);

%% converts solution back into matrix for each time-step
for t = 1:nTimeSteps
	spinSol(:,:,t) = reshape(Y(t,:),N,3);
end

end
