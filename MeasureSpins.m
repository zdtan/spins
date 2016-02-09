function [ multipliers, angles ] = MeasureSpins( spinMat, Jmat, alpha, K, L, M, AFM )
%MeasureSpins - Calculate effective (unshattered) Lagrange multipliers and 
%       angular deviation (in degrees) of spin from local field at every site.
%   Expect multipliers to be non-negative and angles to be near zero,
%   for both AFM and FM cases.
%   ** Distinguish analysis for kagome and triangular sites. **
% Technical note: multipliers (unshattered) defined to be the
%   length of the local field "L_i S_i", as the scalar product gives "L_i S_i^2".

%% system size
N = 8*K*L*M; % 8 sites per unit cell

%% sanity check on size of spins
if size(spinMat,1) ~= N
    disp('== ERROR: Inconsistency in size of spin matrix! ==');
    return
end

%% calculate local field
LFmat = -Jmat*spinMat;

%% return values
multipliers = zeros(1,N);
angles = zeros(1,N);
for n = 1:N
    multipliers(n) = norm(LFmat(n,:),2);
    % find angle using "a.b = |a||b|cos(x)"
    angles(n) = real(acosd( ...
        dot(spinMat(n,:),LFmat(n,:),2)/(norm(spinMat(n,:),2).*norm(LFmat(n,:),2)) ...
        )); % take real part because sometimes get small imaginary part...
end

end