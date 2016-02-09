function [ Jmat ] = AdjacencyMatrix( alpha, K, L, M, AFM )
%AdjacencyMatrix - Constructs N x N adjacency matrix J_{ij} in SPARSE format.
%       OPPOSITE sign to Walker-Walstedt convention, 
%		as we construct a POSITIVE matrix for AFM 
%		and NEGATIVE matrix for FM.
%	Fixed pathology for K, L = 1 case by using X = X + 1 syntax
%       (just need to change for kagome neighbours, 
%       as unit cell already has two layers of apical spins).

%% system size
N = 8*K*L*M; % 8 sites per unit cell

%% setup matrix - useful to check that it is symmetric!
Jmat = zeros(N,N);
for site = 1:8
    for i = 1:K
        for j = 1:L
            for k = 1:M
                n = SiteIndex1D( site, i, j, k, K, L, M );
                INN = 1;
                JNN = 1;
                KNN = 1;
                switch site
                    case 1
                        % kagome neighbours
                        Jmat(n,n+1) = Jmat(n,n+1) + 1;
                        Jmat(n,n+2) = Jmat(n,n+2) + 1;
                        nni = i - INN;
                        if nni < 1
                            nni = nni + K*INN;
                        end
                        nnj = j - JNN;
                        if nnj < 1
                            nnj = nnj + L*JNN;
                        end
                        Jmat(n,SiteIndex1D( site+1, nni, j, k, K, L, M )) = ...
                            Jmat(n,SiteIndex1D( site+1, nni, j, k, K, L, M )) + 1;
                        Jmat(n,SiteIndex1D( site+2, i, nnj, k, K, L, M )) = ...
                            Jmat(n,SiteIndex1D( site+2, i, nnj, k, K, L, M )) + 1;
                        % apical neighbours
                        nnk = k - KNN;
                        if nnk < 1
                            nnk = nnk + M*KNN;
                        end
                        Jmat(n,n+3) = alpha;
                        Jmat(n,SiteIndex1D( site+7, i, j, nnk, K, L, M )) = alpha;
                    case 2
                        % kagome neighbours
                        Jmat(n,n-1) = Jmat(n,n-1) + 1;
                        Jmat(n,n+1) = Jmat(n,n+1) + 1;
                        nni = i + INN;
                        if nni > K*INN
                            nni = nni - K*INN;
                        end
                        nnj = j - JNN;
                        if nnj < 1
                            nnj = nnj + L*JNN;
                        end
                        Jmat(n,SiteIndex1D( site-1, nni, j, k, K, L, M )) = ...
                            Jmat(n,SiteIndex1D( site-1, nni, j, k, K, L, M )) + 1;
                        Jmat(n,SiteIndex1D( site+1, nni, nnj, k, K, L, M )) = ...
                            Jmat(n,SiteIndex1D( site+1, nni, nnj, k, K, L, M )) + 1;
                        % apical neighbours
                        nnk = k - KNN;
                        if nnk < 1
                            nnk = nnk + M*KNN;
                        end
                        Jmat(n,n+2) = alpha;
                        Jmat(n,SiteIndex1D( site+6, i, j, nnk, K, L, M )) = alpha;
                    case 3
                        % kagome neighbours
                        Jmat(n,n-1) = Jmat(n,n-1) + 1;
                        Jmat(n,n-2) = Jmat(n,n-2) + 1;
                        nni = i - INN;
                        if nni < 1
                            nni = nni + K*INN;
                        end
                        nnj = j + JNN;
                        if nnj > L*JNN
                            nnj = nnj - L*JNN;
                        end
                        Jmat(n,SiteIndex1D( site-1, nni, nnj, k, K, L, M )) = ...
                            Jmat(n,SiteIndex1D( site-1, nni, nnj, k, K, L, M )) + 1;
                        Jmat(n,SiteIndex1D( site-2, i, nnj, k, K, L, M )) = ...
                            Jmat(n,SiteIndex1D( site-2, i, nnj, k, K, L, M )) + 1;
                        % apical neighbours
                        nnk = k - KNN;
                        if nnk < 1
                            nnk = nnk + M*KNN;
                        end
                        Jmat(n,n+1) = alpha;
                        Jmat(n,SiteIndex1D( site+5, i, j, nnk, K, L, M )) = alpha;
                    case 4
                        % neighbours are all J2-type
                        Jmat(n,n-1) = alpha;
                        Jmat(n,n-2) = alpha;
                        Jmat(n,n-3) = alpha;
                        Jmat(n,n+1) = alpha;
                        Jmat(n,n+2) = alpha;
                        Jmat(n,n+3) = alpha;
                    case 5
                        % kagome neighbours
                        Jmat(n,n+1) = Jmat(n,n+1) + 1;
                        Jmat(n,n+2) = Jmat(n,n+2) + 1;
                        nni = i + INN;
                        if nni > K*INN
                            nni = nni - K*INN;
                        end
                        nnj = j - JNN;
                        if nnj < 1
                            nnj = nnj + L*JNN;
                        end
                        Jmat(n,SiteIndex1D( site+1, i, nnj, k, K, L, M )) = ...
                            Jmat(n,SiteIndex1D( site+1, i, nnj, k, K, L, M )) + 1;
                        Jmat(n,SiteIndex1D( site+2, nni, nnj, k, K, L, M )) = ...
                            Jmat(n,SiteIndex1D( site+2, nni, nnj, k, K, L, M )) + 1;
                        % apical neighbours
                        Jmat(n,n-1) = alpha;
                        Jmat(n,n+3) = alpha;
                    case 6
                        % kagome neighbours
                        Jmat(n,n+1) = Jmat(n,n+1) + 1;
                        Jmat(n,n-1) = Jmat(n,n-1) + 1;
                        nni = i + INN;
                        if nni > K*INN
                            nni = nni - K*INN;
                        end
                        nnj = j + JNN;
                        if nnj > L*JNN
                            nnj = nnj - L*JNN;
                        end
                        Jmat(n,SiteIndex1D( site+1, nni, j, k, K, L, M )) = ...
                            Jmat(n,SiteIndex1D( site+1, nni, j, k, K, L, M )) + 1;
                        Jmat(n,SiteIndex1D( site-1, i, nnj, k, K, L, M )) = ...
                            Jmat(n,SiteIndex1D( site-1, i, nnj, k, K, L, M )) + 1;
                        % apical neighbours
                        Jmat(n,n-2) = alpha;
                        Jmat(n,n+2) = alpha;
                    case 7
                        % kagome neighbours
                        Jmat(n,n-1) = Jmat(n,n-1) + 1;
                        Jmat(n,n-2) = Jmat(n,n-2) + 1;
                        nni = i - INN;
                        if nni < 1
                            nni = nni + K*INN;
                        end
                        nnj = j + JNN;
                        if nnj > L*JNN
                            nnj = nnj - L*JNN;
                        end
                        Jmat(n,SiteIndex1D( site-1, nni, j, k, K, L, M )) = ...
                            Jmat(n,SiteIndex1D( site-1, nni, j, k, K, L, M )) + 1;
                        Jmat(n,SiteIndex1D( site-2, nni, nnj, k, K, L, M )) = ...
                            Jmat(n,SiteIndex1D( site-2, nni, nnj, k, K, L, M )) + 1;
                        % apical neighbours
                        Jmat(n,n-3) = alpha;
                        Jmat(n,n+1) = alpha;
                    case 8
                        % neighbours are all J2-type
                        nnk = k + KNN;
                        if nnk > M*KNN
                            nnk = nnk - M*KNN;
                        end
                        Jmat(n,n-1) = alpha;
                        Jmat(n,n-2) = alpha;
                        Jmat(n,n-3) = alpha;
                        Jmat(n,SiteIndex1D( site-5, i, j, nnk, K, L, M )) = alpha;
                        Jmat(n,SiteIndex1D( site-6, i, j, nnk, K, L, M )) = alpha;
                        Jmat(n,SiteIndex1D( site-7, i, j, nnk, K, L, M )) = alpha;
                end
            end
        end
    end
end

%% OPPOSITE to Walker-Walstedt sign convention for matrix
if AFM
    % positive matrix for antiferromagnet
else
    Jmat = -Jmat;
end

%% check for symmetry
if all(all(Jmat == Jmat.'))
    % disp('Adjacency matrix is symmetric.');
else
    disp('== ERROR: Adjacency matrix is NOT symmetric! ==');
end

%% return sparse matrix
Jmat = sparse(Jmat);

end
