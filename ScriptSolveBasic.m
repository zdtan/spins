%% To use after MC completes to do MD
% Script plots solutions without averaging over different spin states

% cut-off in plots
Amin = 0.1;

%% automatic initialisation
AFM = true;
sigma = 1;
ID = ['a',num2str(alpha),'-b',num2str(beta),'-s',num2str(sigma),...
    '-',num2str(K),'x',num2str(L),'x',num2str(M),'_',UID];
Tmesh = 10^3;	% default to 1000
%% need to be careful with Tmax
%Tmax = min(beta^2,10^5);	% impose fixed integration time
%Tmax = max(beta,0.01*beta^2);
% for exp1, 1*beta
% for exp1.5, 0.1*beta^(1.5)
% for exp2, 0.01*beta^2
Tmax = max(100,2*beta);

% calculate adjacency matrix just once
Jmat = AdjacencyMatrix( alpha, K, L, M, AFM );

%% load final spin configuration
spinfile = [ID,'/spins-a',num2str(alpha),'-b',num2str(beta),...
    '-s',num2str(sigma), ...
    '-',num2str(K),'x',num2str(L),'x',num2str(M)];
if ~AFM
    spinfile = [spinfile,'_FM'];
end
load([spinfile,'.mat'],'-mat','spinMat');

%% solve EOM
ticMD = tic;
[ T, spinSol ] = IntegrateEOM( spinMat, Jmat, alpha, K, L, M, AFM, Tmax, Tmesh );
[ T, Asites, energy, totalM ] = ...
	AutocorrelationFunction( T, spinSol, spinMat, Jmat, alpha, sigma, K, L, M, AFM );
timeMD = toc(ticMD);
disp(['=======================']);
disp(['Time taken on Molecular Dynamics = ',num2str(timeMD./3600),' hrs']);

%% save full spinSol results (no need to use '-v7.3' as files are NOT large)
spinSolfile = [ID,'/spinSol-a',num2str(alpha),'-b',num2str(beta),...
    '-s',num2str(sigma), ...
    '-',num2str(K),'x',num2str(L),'x',num2str(M),'_basic'];
if ~AFM
    spinSolfile = [spinSolfile,'_FM'];
end
save([spinSolfile,'.mat'],'T','spinSol','-mat');

%% filenames to save data and plots
obsfile = [ID,'/observables-a',num2str(alpha),'-b',num2str(beta),...
    '-s',num2str(sigma), ...
    '-',num2str(K),'x',num2str(L),'x',num2str(M),'_basic'];
if ~AFM
    obsfile = [obsfile,'_FM'];
end

%% check energy drift
energydiff = max(energy) - min(energy);
relErg = energydiff/energy(1);
h = figure;
plot(T,energy);
print(h,'-dpdf',[obsfile,'-energy(',num2str(relErg),').pdf']);
close;

%% check drift in totalM
tm1 = totalM(:,1)/totalM(1,1);
tm2 = totalM(:,2)/totalM(1,2);
tm3 = totalM(:,3)/totalM(1,3);
mdiff1 = max(tm1) - min(tm1);
mdiff2 = max(tm2) - min(tm2);
mdiff3 = max(tm3) - min(tm3);
h = figure;
plot(T,horzcat(tm1,tm2,tm3));
print(h,'-dpdf',[obsfile,'-m(',num2str(max([mdiff1,mdiff2,mdiff3])),').pdf']);
close;

%% 2 kinds of sites
Akagome = horzcat(Asites(:,1),Asites(:,2),Asites(:,3),...
	Asites(:,5),Asites(:,6),Asites(:,7));
Atriangular = horzcat(Asites(:,4),Asites(:,8));
% averaged autocorrelations
Akag = sum(Akagome,2)/6;
Atri = sum(Atriangular,2)/2;
Aavg = sum(Asites,2)/8;

%% save all results (no need to use '-v7.3' as files are NOT large)
save([obsfile,'.mat'],'T','energy','totalM', ...
	'Asites','Akagome','Atriangular','Akag','Atri','Aavg','-mat');

%% make plots
% 8 plots
h = figure;
plot(T,Asites);
set(gca,'FontSize',20);
axis([0,Inf,Amin,1])
legend('1','2','3','4','5','6','7','8',...
	'Location','NorthEastOutside')
print(h,'-dpdf',[obsfile,'-Asites.pdf']);
close;

% 6 plots (kagome sites)
h = figure;
plot(T,Akagome);
set(gca,'FontSize',20);
axis([0,Inf,Amin,1])
legend('1','2','3','5','6','7',...
	'Location','NorthEastOutside')
print(h,'-dpdf',[obsfile,'-Akagome.pdf']);
close;

% 2 plots (triangular sites)
h = figure;
plot(T,Atriangular);
set(gca,'FontSize',20);
axis([0,Inf,Amin,1])
legend('4','8',...
	'Location','NorthEast')
print(h,'-dpdf',[obsfile,'-Atriangular.pdf']);
close;

% plot the various averages
h = figure;
plot(T,horzcat(Akag,Atri,Aavg));
set(gca,'FontSize',20);
axis([0,Inf,Amin,1])
legend('kagome','triangular','average',...
	'Location','NorthEast')
print(h,'-dpdf',[obsfile,'-Akagtriavg.pdf']);
close;

h = figure;
plot(T,horzcat(Akag,Atri));
set(gca,'FontSize',20);
axis([0,Inf,Amin,1])
legend('kagome','triangular',...
	'Location','NorthEast')
print(h,'-dpdf',[obsfile,'-Akagtri.pdf']);
close;

h = figure;
plot(T,Akag);
axis([0,Inf,Amin,1])
set(gca,'FontSize',20);
print(h,'-dpdf',[obsfile,'-Akag.pdf']);
close;

h = figure;
plot(T,Atri);
axis([0,Inf,Amin,1])
set(gca,'FontSize',20);
print(h,'-dpdf',[obsfile,'-Atri.pdf']);
close;

h = figure;
plot(T,Aavg);
axis([0,Inf,Amin,1])
set(gca,'FontSize',20);
print(h,'-dpdf',[obsfile,'-Aavg.pdf']);
close;

