%% taken from 20140826-spinwaves

function [ spinMat, energy ] = GenerateGS( ID, ...
				alpha, Jmat, beta, sigma, K, L, M, AFM, numSweep )
%GenerateGS - Generate ground state configuration from starting configuration
%		obtained from simulated annealing (to inverse temperature beta) using Monte Carlo.
%	Uses steepest descent, loading and saving files to subfolder ID.
%   alpha is ratio of J2 to J1.
%   sigma gives length of triangular spins (kagome spins of unit length).
% should have 10^5 sweeps for 8x8x8 system I think...

%% create directory for saving
[s,mess,messid] = mkdir(ID);
if s == 1
    % successful in creating directory
    if strcmp(messid,'MATLAB:MKDIR:DirectoryExists')
        %disp(['== ABORT: directory already exists! ==']);
        %return
    else
        disp(['== ABORT: directory does not exist! ==']);
        return
    end
else
    disp(['== ABORT: failed to create directory! ==']);
    return
end

%% system size
N = 8*K*L*M; % 8 sites per unit cell

%% load thermal state
MCspinfile = [ID,'/spins-a',num2str(alpha),'-b',num2str(beta),...
    '-s',num2str(sigma), ...
    '-',num2str(K),'x',num2str(L),'x',num2str(M)];
if ~AFM
    MCspinfile = [MCspinfile,'_FM'];
end
load([MCspinfile,'.mat'],'-mat','spinMat');

%% steepest descent
ticID = tic;
% default parameters
delta = 1;
[ spinMat, energy ] = SteepestDescent( spinMat, Jmat, ...
	alpha, delta, sigma, K, L, M, numSweep, AFM );

timeSD = toc(ticID);
disp(['Time taken on steepest descent = ',num2str(timeSD./3600),' hrs']);

%% save and plot energy trajectory for SD
energyfile = [ID,'/energy-a',num2str(alpha), ...
    '-d',num2str(delta),'-s',num2str(sigma), ...
    '-',num2str(K),'x',num2str(L),'x',num2str(M), ...
    '_',num2str(numSweep),'sweeps_SD'];
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

%% save final spin configuration
spinfile = [ID,'/spins-a',num2str(alpha),'-SD', ...
    '-s',num2str(sigma), ...
    '-',num2str(K),'x',num2str(L),'x',num2str(M)];
if ~AFM
    spinfile = [spinfile,'_FM'];
end
save([spinfile,'.mat'],'spinMat',...
	'alpha','sigma','K','L','M','-mat');

end
