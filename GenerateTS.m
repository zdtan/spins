function [ spinMat, energy, acceptance ] = GenerateTS( ID, ...
	alpha, beta, sigma, K, L, M, AFM, constrainM, Jmat, sinit, numSweep, delta0 )
%GenerateTS - Generate thermal spin state configuration using Monte Carlo,
%       saving files to subfolder ID. Save spins after each annealing step.
%   alpha is ratio of J2 to J1.
%   sigma gives length of triangular spins (kagome spins of unit length).
%	beta is the target (inverse) temperature.
% sinit contains information about random number generator before the runs.
% numSweep = #sweeps per spin at each temperature

%% create directory for saving
[s,mess,messid] = mkdir(ID);
if s == 1
    % successful in creating directory
    if strcmp(messid,'MATLAB:MKDIR:DirectoryExists')
        disp(['== WARNING: directory already exists! ==']);
        %return % comment out to let code go on
    end
else
    disp(['== ABORT: failed to create directory! ==']);
    return
end

%% system size
N = 8*K*L*M; % 8 sites per unit cell

%% create random spins
spinMat = GenerateRS(sigma,K,L,M);

%% Monte Carlo simulated annealing
% schedule based on lowering temperature by e each time (except the last time)
ticID = tic;
nStepsTot = ceil(log(beta)) + 1;
acceptance = zeros(1,nStepsTot);	% save and plot a graph later
% annealing steps
nSteps = nStepsTot;
while nSteps > 1	%skip to final step if only 1 step left
	% tempBeta starts from 1 and gets bigger with each step
	tempBeta = exp(nStepsTot-nSteps);	% multiply by e each time
	delta = delta0/sqrt(tempBeta);	% scale as sqrt(T)
	%delta = delta0/(tempBeta^(2/3));	% scale as T^{2/3}? this increases acceptance ratio too much!
	%delta = deltaMin + (nSteps/nStepsTot)*(deltaMax-deltaMin);	% delta gets smaller with each step
	[ spinMat, energy, acc ] = MonteCarlo( spinMat, Jmat, ...
		alpha, tempBeta, delta, sigma, K, L, M, numSweep, AFM, constrainM );
	acceptance(nStepsTot-nSteps+1) = acc;
	nSteps = nSteps - 1;
	% save spin configuration at each tempBeta (and magnetisation)
	totalM = sum(spinMat,1);
	spinfile = [ID,'/spins-a',num2str(alpha),'-b',num2str(tempBeta),...
		'-s',num2str(sigma), ...
		'-',num2str(K),'x',num2str(L),'x',num2str(M),'-m',num2str(norm(totalM,2))];
	if ~AFM
		spinfile = [spinfile,'_FM'];
	end
	save([spinfile,'.mat'],'spinMat','totalM', ...
		'alpha','beta','sigma','K','L','M','constrainM','-mat');
	%% save and plot energy trajectory
	energyfile = [ID,'/energy','-a',num2str(alpha),'-b',num2str(tempBeta), ...
    		'-d',num2str(delta),'-s',num2str(sigma), ...
    		'-',num2str(K),'x',num2str(L),'x',num2str(M), ...
    		'_',num2str(numSweep),'sweeps_',num2str(acc),'acceptance'];
	if ~AFM
    		energyfile = [energyfile,'_FM'];
	end
	save([energyfile,'.mat'],'energy','-mat');
	% linear plot
	h = figure;
	plot(energy);
	axis([1,numSweep+1,-Inf,Inf])
	set(gca,'FontSize',20);
	xlabel('sweeps per spin','FontSize',30)
	ylabel('energy','FontSize',30)
	print(h,'-dpdf',[energyfile,'.pdf']);
	close;
	% log plot
	h = figure;
	plot(log10(energy));
	axis([1,numSweep+1,-Inf,Inf])
	set(gca,'FontSize',20);
	xlabel('sweeps per spin','FontSize',30)
	ylabel('log_{10} (energy)','FontSize',30)
	print(h,'-dpdf',[energyfile,'-LOG','.pdf']);
	close;
end
% final step at target temperature
if nStepsTot <= 1
	delta = delta0;
else
	delta = delta0/sqrt(beta);	% scale as sqrt(T)
	%numSweep = nSweep*max(nStepsTot,10);	% don't cap at 10x?
end
[ spinMat, energy, acc ] = MonteCarlo( spinMat, Jmat, ...
	alpha, beta, delta, sigma, K, L, M, numSweep, AFM, constrainM );
acceptance(nStepsTot)= acc;

timeMC = toc(ticID);
disp(['Time taken on Monte Carlo = ',num2str(timeMC./3600),' hrs']);

%% save final spin configuration + acceptance ratio + rng settings
spinfile = [ID,'/spins-a',num2str(alpha),'-b',num2str(beta),...
    '-s',num2str(sigma), ...
    '-',num2str(K),'x',num2str(L),'x',num2str(M)];
if ~AFM
    spinfile = [spinfile,'_FM'];
end
save([spinfile,'.mat'],'spinMat','acceptance','sinit',...
	'alpha','beta','sigma','K','L','M','constrainM','-mat');
%~ % save initial random number generator settings
%~ save([ID,'/rngsettings.mat'],'sinit','-mat');

%% save and plot energy trajectory
energyfile = [ID,'/energy','-a',num2str(alpha),'-b',num2str(beta), ...
    		'-d',num2str(delta),'-s',num2str(sigma), ...
    		'-',num2str(K),'x',num2str(L),'x',num2str(M), ...
    		'_',num2str(numSweep),'sweeps_',num2str(acc),'acceptance'];
if ~AFM
	energyfile = [energyfile,'_FM'];
end
save([energyfile,'.mat'],'energy','-mat');
% linear plot
h = figure;
plot(energy);
axis([1,numSweep+1,-Inf,Inf])
set(gca,'FontSize',20);
xlabel('sweeps per spin','FontSize',30)
ylabel('energy','FontSize',30)
print(h,'-dpdf',[energyfile,'.pdf']);
close;
% log plot
h = figure;
plot(log10(energy));
axis([1,numSweep+1,-Inf,Inf])
set(gca,'FontSize',20);
xlabel('sweeps per spin','FontSize',30)
ylabel('log_{10} (energy)','FontSize',30)
print(h,'-dpdf',[energyfile,'-LOG','.pdf']);
close;
    
if nStepsTot <= 1
	% don't plot graph otherwise get error!
else
	h = figure;
	plot(acceptance,'--o');
	axis([1,nStepsTot,0,1])
	set(gca,'FontSize',20);
	set(gca,'XTick',[1:nStepsTot])
	ylabel('acceptance','FontSize',30)
	print(h,'-dpdf',[spinfile,'-acceptance.pdf']);
	close;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% call MeasureSpins to check Lagrange multipliers and angular deviations
[ multipliers, angles ] = MeasureSpins( spinMat, Jmat, alpha, K, L, M, AFM );

%% visualise and save data
% break up Lagrange multipliers based on site
site1 = multipliers(1:8:N);
site2 = multipliers(2:8:N);
site3 = multipliers(3:8:N);
site4 = multipliers(4:8:N);
site5 = multipliers(5:8:N);
site6 = multipliers(6:8:N);
site7 = multipliers(7:8:N);
site8 = multipliers(8:8:N);

% filename to use
measurefile = [ID,'/measure-a',num2str(alpha),'-b',num2str(beta),...
    '-s',num2str(sigma), ...
    '-',num2str(K),'x',num2str(L),'x',num2str(M)];
if ~AFM
    measurefile = [measurefile,'_FM'];
end
save([measurefile,'.mat'],'site1','site2','site3','site4',...
    'site5','site6','site7','site8','multipliers','angles','-mat');

h = figure;
plot(horzcat(site1,site2,site3,site5,site6,site7));
axis([1,3*N/4,-Inf,Inf])
set(gca,'FontSize',20);
xlabel('kagome site','FontSize',30)
ylabel('multiplier','FontSize',30)
print(h,'-dpdf',[measurefile,'-kagome.pdf']);
close;

h = figure;
plot(horzcat(site4,site8));
axis([1,N/4,-Inf,Inf])
set(gca,'FontSize',20);
xlabel('triangular site','FontSize',30)
ylabel('multiplier','FontSize',30)
print(h,'-dpdf',[measurefile,'-triangular.pdf']);
close;

h = figure;
plot(angles);
axis([1,N,-Inf,Inf])
set(gca,'FontSize',20);
xlabel('site','FontSize',30)
ylabel('angular deviation (degrees)','FontSize',30)
maxAngle = max(angles);
print(h,'-dpdf',[measurefile,'-angles(',num2str(maxAngle),').pdf']);
close;

disp(['Minimum multipliers = ',num2str(min(multipliers)),...
    ' and Maximum deviation = ',num2str(maxAngle)]);

end
