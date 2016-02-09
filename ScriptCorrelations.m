%% Calculate correlations for triangular lattice spins in interplane direction

% automatic initialisation
ID = ['a',num2str(alpha),'-b',num2str(beta),'-s',num2str(sigma),...
    '-',num2str(K),'x',num2str(L),'x',num2str(M),'_',UID];
AFM = true;

folder = [ID,'/correlations2'];
%% "create" directory for saving
[s,mess,messid] = mkdir(folder);
if s == 1
    % successful in creating directory
    if strcmp(messid,'MATLAB:MKDIR:DirectoryExists')
        % directory already exists
    else
        %disp(['== ABORT: directory does not previously exist! ==']);
        %return
    end
else
    disp(['== ABORT: failed to create directory! ==']);
    return
end

%% load final spin configuration
spinfile = [ID,'/spins-a',num2str(alpha),'-b',num2str(beta),...
    '-s',num2str(sigma), ...
    '-',num2str(K),'x',num2str(L),'x',num2str(M)];
if ~AFM
    spinfile = [spinfile,'_FM'];
end
load([spinfile,'.mat'],'-mat','spinMat');

%%% calculate chiralities
%[ bicap, uncap, bsort, usort, b2avg, u2avg ] = Chirality( spinMat, K, L, M, folder );
%[ tsort, t2avg ] = ChiralityTri( spinMat, K, L, M, folder );

%% calculate correlations
CorrelationsAvg( spinMat, K, L, M, 'api-col', folder );
%for code = {'kag-dir-a','kag-dir-b','kag-dir-b-a',...
%	'tri-dir-a','tri-dir-b','tri-dir-b-a','api-col','kag-col'}
%	Correlations( spinMat, K, L, M, code{1}, folder );
%end
%metaCODE = 'kag-dir'; CorrelationsCombine( metaCODE, folder );
%metaCODE = 'tri-dir'; CorrelationsCombine( metaCODE, folder );
