function [ spinMat, energy, acceptance ] = MonteCarlo( iSpinMat, Jmat, ...
    alpha, beta, delta, sigma, K, L, M, numSweep, AFM, constrainM )
%MonteCarlo - Monte Carlo simulation of spins on extended kagome lattice.
%   alpha = J2/J1, beta = |J1/kT|, delta is "strength" of spin disturbance
%   for each Metropolis step. 
%       Note: for ferromagnetic J1, J2 set AFM = false. beta > 0 always.
%   Number of unit cells in kagome plane is K-by-L, 
%   and number of layers of unit cells is M.
%   iSpinMat is Nx3 matrix of initial spins.
%   numSweep is number of Metropolis steps per spin on average.
%% energy calculated as energy per spin - (-3/4 - (alpha*sigma)^2/2) [diff with GS energy]
% (old) set constrainM to true (total magnetisation near zero) or false
% (current) constrainM is coefficient for term added to Hamiltonian - zero adds no term
% acceptance ratio returned is averaged over the 2nd half of the run.
%%%%%%%%%%%%%%%%%%%%%%%

%% system size
N = 8*K*L*M; % 8 sites per unit cell

%%% create directory for saving
%[s,mess,messid] = mkdir(ID);
%if s == 1
%    % successful in creating directory
%    if strcmp(messid,'MATLAB:MKDIR:DirectoryExists')
%        %~ disp(['== WARNING: directory already exists! ==']);
%        %~ return
%    end
%else
%    disp(['== ABORT: failed to create directory! ==']);
%    return
%end

%% initialise result
spinMat = iSpinMat;
totalM = sum(spinMat,1);

%% system energy per spin
energy = zeros(1,numSweep+1);

% calculate energy of system at the very start
energy(1) = 0.5*sum(dot(spinMat,Jmat*spinMat,2))/N ...
	- (-3/4 - (alpha*sigma)^2/2) ...
	+ constrainM*sum(abs(totalM).^2)/N;
%~ if constrainM
	%~ energy(1) = energy(1) + sum(abs(totalM).^2)/N;
%~ end

%% spinVec modification in Monte Carlo move
% perturb spin (length S) by adding spin (length delta.*S),
%  and scaling back final answer to length S
function spin2 = disturbSpin(spin1,delta,spinX) % spinX is unit-length spin
	spinLength = norm(spin1,2);
	tempSpin = spin1 + delta.*spinLength.*spinX;
	spin2 = tempSpin.*spinLength./norm(tempSpin,2);
end

%% begin sweeps
disp(['SET alpha = ',num2str(alpha), ...
	', beta = ',num2str(beta),' , sigma = ',num2str(sigma),...
	' , delta = ',num2str(delta),...
	' for ',num2str(numSweep),' sweeps']);
acceptance = 0;
for t = 1:numSweep
    %disp(['...Monte Carlo run ',num2str(t),' of ',num2str(numSweep),'...']);
    % tic;
    accept = 0; % reset spin-flip acceptance count for each sweep
    energy(t+1) = energy(t); % initialise energy from previous sweep
    
    %% Re-calculating energy is costly, though once a sweep is probably ok...
    %~ % calculate local field at start of sweep
    %~ LFmat = -Jmat*spinMat;
    %~ % calculate energy of system at start of each sweep
    %~ energy(t+1) = 0.5*sum(dot(spinMat,-LFmat,2))/N ...
	%~ - (-3/4 - (alpha*sigma)^2/2);
    %~ if constrainM
	%~ totalM = sum(spinMat,1);
	%~ energy(t+1) = energy(t+1) + sum(abs(totalM).^2)/N;
    %~ end
    
    % generate random spins for Monte Carlo moves in this sweep
    spinRand = RandomSpins(N);
    % generate random numbers to compare against weight
    rprob = rand(N,1);
    
    %% choose a site once on average for Metropolis algorithm
    for s = 1:N
        n = randi(N);
	   spin1 = spinMat(n,:);
        spin2 = disturbSpin(spin1,delta,spinRand(s,:)); % modified spin candidate
	   localField = -Jmat(n,:)*spinMat;
	   	%localField = -transpose(Jmat(:,n))*spinMat;	% use columns of sparse matrix instead?
        % Metropolis algorithm
	   totalM2 = totalM - spin1 + spin2;
	   weight = dot((spin2-spin1),localField,2) ...
		  + constrainM*(sum(abs(totalM).^2)-sum(abs(totalM2).^2));
	%~ if constrainM
		%~ totalM2 = totalM - spin1 + spin2;
		%~ weight = weight + sum(abs(totalM).^2) - sum(abs(totalM2).^2);
	%~ end
	   if (weight < 0) && (rprob(s) > exp(beta*weight))
		% move rejected, no flipping, no change in energy
	   else
		% move accepted
		accept = accept + 1;
		spinMat(n,:) = spin2;
		
		%% Skip this expensive step!!!
		%~ % update local field
		%~ LFmat = -Jmat*spinMat;
		
		% update energy and totalM
		energy(t+1) = energy(t+1) - weight/N;
		totalM = totalM2;
		%~ if constrainM
			%~ totalM = totalM2;
		%~ end
	end
	% continue until sites visited once on average
    end
    if t > (numSweep/2)
	acceptance = acceptance + accept;	% update overall acceptance
    end
    % continue numSweep times
    %% display acceptance percentage for current run
    %acceptPercent = accept/N;
    %disp(['Acceptance for this run = ', num2str(100*acceptPercent), '%']);
    % disp(['Time taken = ',num2str(toc),' secs']);
end
acceptance = acceptance/(N*numSweep/2);
disp(['Overall acceptance (for 2nd half of run) = ', num2str(100*acceptance), '%']);

end
