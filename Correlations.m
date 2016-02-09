function [ sites, dotp, avgdotp ] = Correlations( spinMat, K, L, M, CODE, folder )
%Correlations - Calculate spin scalar products and plot
%	sites is the x-axis, dotp is the y-axis
% save in 'folder'

%% create directory for details
folderDetail = [folder,'/detail']; %[CODE,'_detail'];
[s,mess,messid] = mkdir(folderDetail);
if s == 1
    % successful in creating directory
    if strcmp(messid,'MATLAB:MKDIR:DirectoryExists')
        %disp(['== ABORT: directory already exists! ==']);
        %return
    else
        %disp(['== ABORT: directory does not exist! ==']);
        %return
    end
else
    disp(['== ABORT: failed to create directory! ==']);
    return
end

%% system size
N = 8*K*L*M; % 8 sites per unit cell

%% sanity check on size of spins
if size(spinMat,1) ~= N
    disp('== ERROR: Inconsistency in size of spin matrix! ==');
    return
end

%% figure out which correlation to calculate
switch CODE
    case 'kag-dir-a'
        % average over kagome lattice layers - in 'a'-dir
        sites = [1:2*K];
        samples = [1:L*2*M];
        dotp = zeros(length(sites),length(samples));
        % settle sample by sample
        for j = [1:L]
            for k = [1:M]
                tempspins1 = zeros(length(sites),3);
                tempspins2 = zeros(length(sites),3);
                for i = [1:K]
                    tempspins1(2*i-1,:) = spinMat(SiteIndex1D(1,i,j,k,K,L,M),:);
                    tempspins1(2*i,:) = spinMat(SiteIndex1D(2,i,j,k,K,L,M),:);
                    tempspins2(2*i-1,:) = spinMat(SiteIndex1D(7,i,j,k,K,L,M),:);
                    tempspins2(2*i,:) = spinMat(SiteIndex1D(6,i,j,k,K,L,M),:);
                end
                dotp(:,j+L*(2*(k-1))) = dot(tempspins1,repmat(tempspins1(1,:),length(sites),1),2);
                dotp(:,j+L*(2*(k-1)+1)) = dot(tempspins2,repmat(tempspins2(1,:),length(sites),1),2);
            end
        end
    case 'kag-dir-b'
        % average over kagome lattice layers - in 'b'-dir
        sites = [1:2*L];
        samples = [1:K*2*M];
        dotp = zeros(length(sites),length(samples));
        % settle sample by sample
        for i = [1:K]
            for k = [1:M]
                tempspins1 = zeros(length(sites),3);
                tempspins2 = zeros(length(sites),3);
                for j = [1:L]
                    tempspins1(2*j-1,:) = spinMat(SiteIndex1D(1,i,j,k,K,L,M),:);
                    tempspins1(2*j,:) = spinMat(SiteIndex1D(3,i,j,k,K,L,M),:);
                    tempspins2(2*j-1,:) = spinMat(SiteIndex1D(5,i,j,k,K,L,M),:);
                    tempspins2(2*j,:) = spinMat(SiteIndex1D(6,i,j,k,K,L,M),:);
                end
                dotp(:,i+K*(2*(k-1))) = dot(tempspins1,repmat(tempspins1(1,:),length(sites),1),2);
                dotp(:,i+K*(2*(k-1)+1)) = dot(tempspins2,repmat(tempspins2(1,:),length(sites),1),2);
            end
        end
    case 'kag-dir-b-a'
        % average over kagome lattice layers - in 'b-a'-dir
        sites = [1:2*L];
        samples = [1:K*2*M];
        dotp = zeros(length(sites),length(samples));
        % settle sample by sample
        for i = [1:K]
            for k = [1:M]
                tempspins1 = zeros(length(sites),3);
                tempspins2 = zeros(length(sites),3);
                newi = i; % start with this position and gradually shift to go in 'b-a'
                for j = [1:L]
                    tempspins1(2*j-1,:) = spinMat(SiteIndex1D(2,newi,j,k,K,L,M),:);
                    tempspins1(2*j,:) = spinMat(SiteIndex1D(3,newi,j,k,K,L,M),:);
                    tempspins2(2*j-1,:) = spinMat(SiteIndex1D(5,newi,j,k,K,L,M),:);
                    tempspins2(2*j,:) = spinMat(SiteIndex1D(7,newi,j,k,K,L,M),:);
                    if newi == 1
                        newi = K;
                    else
                        newi = newi-1;
                    end
                end
                dotp(:,i+K*(2*(k-1))) = dot(tempspins1,repmat(tempspins1(1,:),length(sites),1),2);
                dotp(:,i+K*(2*(k-1)+1)) = dot(tempspins2,repmat(tempspins2(1,:),length(sites),1),2);
            end
        end
    case 'tri-dir-a'
        % average over triangular lattice layers - in 'a'-dir
        sites = [1:K];
        samples = [1:L*2*M];
        dotp = zeros(length(sites),length(samples));
        % settle sample by sample
        for j = [1:L]
            for k = [1:M]
                tempspins1 = zeros(length(sites),3);
                tempspins2 = zeros(length(sites),3);
                for i = [1:K]
                    tempspins1(i,:) = spinMat(SiteIndex1D(4,i,j,k,K,L,M),:);
                    tempspins2(i,:) = spinMat(SiteIndex1D(8,i,j,k,K,L,M),:);
                end
                dotp(:,j+L*(2*(k-1))) = dot(tempspins1,repmat(tempspins1(1,:),length(sites),1),2);
                dotp(:,j+L*(2*(k-1)+1)) = dot(tempspins2,repmat(tempspins2(1,:),length(sites),1),2);
            end
        end
    case 'tri-dir-b'
        % average over triangular lattice layers - in 'b'-dir
        sites = [1:L];
        samples = [1:K*2*M];
        dotp = zeros(length(sites),length(samples));
        % settle sample by sample
        for i = [1:K]
            for k = [1:M]
                tempspins1 = zeros(length(sites),3);
                tempspins2 = zeros(length(sites),3);
                for j = [1:L]
                    tempspins1(j,:) = spinMat(SiteIndex1D(4,i,j,k,K,L,M),:);
                    tempspins2(j,:) = spinMat(SiteIndex1D(8,i,j,k,K,L,M),:);
                end
                dotp(:,i+K*(2*(k-1))) = dot(tempspins1,repmat(tempspins1(1,:),length(sites),1),2);
                dotp(:,i+K*(2*(k-1)+1)) = dot(tempspins2,repmat(tempspins2(1,:),length(sites),1),2);
            end
        end
    case 'tri-dir-b-a'
        % average over triangular lattice layers - in 'b-a'-dir
        sites = [1:L];
        samples = [1:K*2*M];
        dotp = zeros(length(sites),length(samples));
        % settle sample by sample
        for i = [1:K]
            for k = [1:M]
                tempspins1 = zeros(length(sites),3);
                tempspins2 = zeros(length(sites),3);
                newi = i; % start with this position and gradually shift to go in 'b-a'
                for j = [1:L]
                    tempspins1(j,:) = spinMat(SiteIndex1D(4,newi,j,k,K,L,M),:);
                    tempspins2(j,:) = spinMat(SiteIndex1D(8,newi,j,k,K,L,M),:);
                    if newi == 1
                        newi = K;
                    else
                        newi = newi-1;
                    end
                end
                dotp(:,i+K*(2*(k-1))) = dot(tempspins1,repmat(tempspins1(1,:),length(sites),1),2);
                dotp(:,i+K*(2*(k-1)+1)) = dot(tempspins2,repmat(tempspins2(1,:),length(sites),1),2);
            end
        end
    case 'api-col'
        % average over apical columns
        sites = [1:2*M];
        samples = [1:K*L];
        dotp = zeros(length(sites),length(samples));
        % settle sample by sample
        for i = [1:K]
            for j = [1:L]
                tempspins = zeros(length(sites),3);
                for k = [1:M]
                    tempspins(2*(k-1)+1,:) = spinMat(SiteIndex1D(4,i,j,k,K,L,M),:);
                    tempspins(2*k,:) = spinMat(SiteIndex1D(8,i,j,k,K,L,M),:);
                end
                dotp(:,i+(j-1)*K) = dot(tempspins,repmat(tempspins(1,:),length(sites),1),2);
            end
        end
    case 'kag-col'
        % average over kagome columns (need to skip layers)
        sites = [1:M];
        samples = [1:K*L];
        dotp1 = zeros(length(sites),length(samples));
        dotp2 = zeros(length(sites),length(samples));
        dotp3 = zeros(length(sites),length(samples));
        dotp5 = zeros(length(sites),length(samples));
        dotp6 = zeros(length(sites),length(samples));
        dotp7 = zeros(length(sites),length(samples));
        % settle sample by sample
        for i = [1:K]
            for j = [1:L]
                tempspins1 = zeros(length(sites),3);
                tempspins2 = zeros(length(sites),3);
                tempspins3 = zeros(length(sites),3);
                tempspins5 = zeros(length(sites),3);
                tempspins6 = zeros(length(sites),3);
                tempspins7 = zeros(length(sites),3);
                for k = [1:M]
                    tempspins1(k,:) = spinMat(SiteIndex1D(1,i,j,k,K,L,M),:);
                    tempspins2(k,:) = spinMat(SiteIndex1D(2,i,j,k,K,L,M),:);
                    tempspins3(k,:) = spinMat(SiteIndex1D(3,i,j,k,K,L,M),:);
                    tempspins5(k,:) = spinMat(SiteIndex1D(5,i,j,k,K,L,M),:);
                    tempspins6(k,:) = spinMat(SiteIndex1D(6,i,j,k,K,L,M),:);
                    tempspins7(k,:) = spinMat(SiteIndex1D(7,i,j,k,K,L,M),:);
                end
                dotp1(:,i+(j-1)*K) = dot(tempspins1,repmat(tempspins1(1,:),length(sites),1),2);
                dotp2(:,i+(j-1)*K) = dot(tempspins2,repmat(tempspins2(1,:),length(sites),1),2);
                dotp3(:,i+(j-1)*K) = dot(tempspins3,repmat(tempspins3(1,:),length(sites),1),2);
                dotp5(:,i+(j-1)*K) = dot(tempspins5,repmat(tempspins5(1,:),length(sites),1),2);
                dotp6(:,i+(j-1)*K) = dot(tempspins6,repmat(tempspins6(1,:),length(sites),1),2);
                dotp7(:,i+(j-1)*K) = dot(tempspins7,repmat(tempspins7(1,:),length(sites),1),2);
            end
        end
        dotp = (dotp1 + dotp2 + dotp3 + dotp5 + dotp6 + dotp7)/6;
end

X = 0:1:(length(sites)-1);

%% plots (show average)
% average over samples
avgdotp = sum(dotp,2)/length(samples);
h = figure;
plot(X,avgdotp,'--ko');
if strcmp(CODE,'api-col')
    axis([-Inf,Inf,min(avgdotp),1])
    set(gca,'YTick',linspace(min(avgdotp),1,5))
else
    axis([-Inf,Inf,-1,1])
end
%semilogy(avgdotp);
set(gca,'FontSize',20);
xlabel('i','FontSize',30)
ylabel('\langle S_0 \cdot S_i \rangle','FontSize',30)
grid on
print(h,'-dpng',[folder,'/correlations-avg_',CODE,'.png']);
close

%% plots (show spread)
h = figure;
plot(X,dotp);
if strcmp(CODE,'api-col')
    axis([-Inf,Inf,min(min(dotp)),1])
    set(gca,'YTick',linspace(min(min(dotp)),1,5))
else
    axis([-Inf,Inf,-1,1])
end
%semilogy(dotp);
set(gca,'FontSize',20);
xlabel('i','FontSize',30)
ylabel('S_0 \cdot S_i','FontSize',30)
print(h,'-dpng',[folder,'/correlations-spread_',CODE,'.png']);
close

%% plots (nearest-neighbour distribution)
% define what nearest-neighbour means
% if strcmp(CODE,'tri-dir-b-a') || strcmp(CODE,'kag-dir-b-a')
%     dotpnn = dotp(2,:);
% else
%     dotpnn = horzcat(dotp(2,:),dotp(length(sites),:));
% end
% dotpnn = dotp(2,:);
% dotpnn = dotp(length(sites),:); % this seems problematic... NO! it's correct!
if length(sites) > 1
    if length(sites) == 2
        dotpnn = dotp(2,:);
    else
        dotpnn = horzcat(dotp(2,:),dotp(length(sites),:));
    end
    % only plot if dotpnn defined
    h = figure;
    plot(sort(dotpnn));
    if strcmp(CODE,'api-col')
        axis([-Inf,Inf,min(dotpnn),1])
        set(gca,'YTick',linspace(min(dotpnn),1,5))
    else
        axis([-Inf,Inf,-1,1])
    end
    set(gca,'FontSize',20);
    xlabel('nearest neighbours','FontSize',30)
    ylabel('S_0 \cdot S_1','FontSize',30)
    print(h,'-dpng',[folder,'/correlations-nn_',CODE,'.png']);
    close
else
    dotpnn = 0; % avoid index out of bounds errors
end

% %% plots (show for individual samples)
% for samp = 1:length(samples)
%     h = figure;
%     plot(dotp(:,samp));
%     if strcmp(CODE,'api-col')
%         axis([-Inf,Inf,min(dotp(:,samp)),1])
%         set(gca,'YTick',linspace(min(dotp(:,samp)),1,5))
%     else
%         axis([-Inf,Inf,-1,1])
%     end
%     %semilogy(dotp);
%     set(gca,'FontSize',20);
%     print(h,'-dpng',[folderDetail,'/correlations-_',CODE,'(',num2str(samp),').png']);
%     close
% end

%% save all data
save([folder,'/correlations-_',CODE,'.mat'], ...
    'dotp','avgdotp','dotpnn','sites','K','L','M','-mat');

end
