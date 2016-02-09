%% Script loads solutions and does averaging over different spin states
% iterates through subfolders
%% Also produces log(autocorrelation) plots

% cut-off in plots
Amin = 0;

%% automatic initialisation
sigma = 1;
AFM = true;
IDroot = ['a',num2str(alpha),'-b',num2str(beta),'-s',num2str(sigma),...
    '-',num2str(K),'x',num2str(L),'x',num2str(M),'_','*'];

% loop through subfolders
IDs = dir(IDroot);	% match wildcard
count = size(IDs,1);
%count = 0;
AkagTot = [];	% start with empty matrix and build up
AtriTot = [];	% start with empty matrix and build up
for sys = 1:count
%for ID = IDs
	%obsfile = [ID{1},'/observables-a',num2str(alpha),'-b',num2str(beta),...
	obsfile = [IDs(sys).name,'/observables-a',num2str(alpha),'-b',num2str(beta),...
		'-s',num2str(sigma), ...
		'-',num2str(K),'x',num2str(L),'x',num2str(M),'_basic'];
	if ~AFM
		obsfile = [obsfile,'_FM'];
	end
	load([obsfile,'.mat'],'-mat','T','Akag','Atri');
	%count = count + 1;
	AkagTot = horzcat(AkagTot,Akag);
	AtriTot = horzcat(AtriTot,Atri);
end

% average over spin states
Akag = sum(AkagTot,2)/count;
Atri = sum(AtriTot,2)/count;
Aavg = sum(horzcat(3*Akag,1*Atri),2)/4;

% create directory for saving
folder = ['AVG_a',num2str(alpha),'-b',num2str(beta),'-s',num2str(sigma),...
    '-',num2str(K),'x',num2str(L),'x',num2str(M)];
[s,mess,messid] = mkdir(folder);
if s == 1
    % successful in creating directory
    if strcmp(messid,'MATLAB:MKDIR:DirectoryExists')
        % no problem either way
    end
else
    disp(['== ABORT: failed to create directory! ==']);
    return
end
filename = [folder,'/observables-a',num2str(alpha),'-b',num2str(beta),...
	'-s',num2str(sigma), ...
	'-',num2str(K),'x',num2str(L),'x',num2str(M),'_basic_',num2str(count),'states'];
if ~AFM
    filename = [filename,'_FM'];
end
save([filename,'.mat'],'T','Aavg','Akag','Atri','-mat');
Amatrix = horzcat(T,Akag,Atri,Aavg);
save([filename,'.txt'],'Amatrix','-ascii','-double');
%save([filename,'.csv'],'Amatrix','-ascii','-double','-tabs');
dlmwrite([filename,'.csv'],Amatrix,'precision',20);
AmatrixLOG = horzcat(T,log(Akag),log(Atri),log(Aavg));
dlmwrite([filename,'-LOG.csv'],AmatrixLOG,'precision',20);

% plot the various averages
h = figure;
plot(T,horzcat(Akag,Atri,Aavg));
set(gca,'FontSize',20);
axis([0,Inf,Amin,1])
legend('kagome','triangular','average',...
	'Location','NorthEast')
print(h,'-dpdf',[filename,'-Akagtriavg.pdf']);
close;

h = figure;
plot(T,horzcat(Akag,Atri));
set(gca,'FontSize',20);
axis([0,Inf,Amin,1])
legend('kagome','triangular',...
	'Location','NorthEast')
print(h,'-dpdf',[filename,'-Akagtri.pdf']);
close;

h = figure;
plot(T,Akag);
axis([0,Inf,Amin,1])
set(gca,'FontSize',20);
print(h,'-dpdf',[filename,'-Akag.pdf']);
close;

h = figure;
plot(T,Atri);
axis([0,Inf,Amin,1])
set(gca,'FontSize',20);
print(h,'-dpdf',[filename,'-Atri.pdf']);
close;

h = figure;
plot(T,Aavg);
axis([0,Inf,Amin,1])
set(gca,'FontSize',20);
print(h,'-dpdf',[filename,'-Aavg.pdf']);
close;

% plot the various averages (LOG plots)
h = figure;
semilogy(T,Akag);
axis([0,Inf,Amin,1])
set(gca,'FontSize',20);
print(h,'-dpdf',[filename,'-Akag_LOG.pdf']);
close;

h = figure;
semilogy(T,Atri);
axis([0,Inf,Amin,1])
set(gca,'FontSize',20);
print(h,'-dpdf',[filename,'-Atri_LOG.pdf']);
close;

h = figure;
semilogy(T,Aavg);
axis([0,Inf,Amin,1])
set(gca,'FontSize',20);
print(h,'-dpdf',[filename,'-Aavg_LOG.pdf']);
close;

h = figure;
semilogy(T,horzcat(Akag,Atri,Aavg));
set(gca,'FontSize',20);
axis([0,Inf,Amin,1])
legend('kagome','triangular','average',...
	'Location','NorthEast')
print(h,'-dpdf',[filename,'-Akagtriavg_LOG.pdf']);
close;

h = figure;
semilogy(T,horzcat(Akag,Atri));
set(gca,'FontSize',20);
axis([0,Inf,Amin,1])
legend('kagome','triangular',...
	'Location','NorthEast')
print(h,'-dpdf',[filename,'-Akagtri_LOG.pdf']);
close;

