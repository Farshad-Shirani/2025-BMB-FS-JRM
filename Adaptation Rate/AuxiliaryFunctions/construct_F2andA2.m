function [F, A] = construct_F2andA2(solution, variableLocation)
%--- This function returnes the value of the ODE vector field at the given solution matrices. The ordering is in x1 direction

global V_s V_u U kappa K Q_opt R

segmentSize = numel(solution(:,1)); % number of (A)ij segments
numNonZero = 3^2 * segmentSize; %number of non-zero elements in A.
odeSystemSize = numel(solution);

N = solution(:,1);
Q = solution(:,2);
V = solution(:,3);

%% Calculating nonlinear terms =============================================
Qm2kVu = Q - 2 * kappa * V_u;

V_squared = V.^2;
Vp2Vu = V + 2 * V_u;
Vp2Vu_squared = Vp2Vu.^2;
VpVp2Vu = V + Vp2Vu;
VpVp2Vu_squared = VpVp2Vu.^2;

Q_opt_squared = Q_opt(:).^2;
Q_squared = Q.^2;
Q_cubed = Q_squared .* Q;
Q_optQ = Q_opt(:) .* Q;
twoQ_optV = 2 * Q_opt(:) .* V;
Q_optQ_squared = Q_opt(:) .* Q_squared;
Q_opt_squaredQ = Q_opt_squared .* Q;
QmQ_opt = Q - Q_opt(:);
VmQ_squared = V - Q_squared;

C = sqrt(2*V_u) * exp( kappa^2 * V_u ) ./ sqrt(VpVp2Vu);
M = exp( -(2 * kappa * V_u)^2 ./ (2 * VpVp2Vu) );
L = ( V .* Qm2kVu + Vp2Vu .* Q ) ./ VpVp2Vu;
S = V .* Vp2Vu ./ VpVp2Vu + L .* (L - 2*Q);
E = ( twoQ_optV + 2 * Q_optQ_squared - Q_opt_squaredQ - 3 * V .* Q - Q_cubed ) / (2*V_s);
Y = U + ( twoQ_optV .* Q - 2 * Q_opt(:) .* Q_cubed - Q_opt_squared .* VmQ_squared - 3 * V_squared + Q_cubed .* Q ) / (2*V_s);

R_K =  R(:) ./ K(:);
RN_K = R_K .* N;
RCN_K = RN_K .* C;
RMN_K = RN_K .* M;
RMCN_K = RCN_K .* M;
RMC_K = R_K .* M .* C;

G = R(:) - RMCN_K - ( QmQ_opt.^2 + V ) / (2 * V_s);
RmG =  R(:) - G;
H = RmG .* Q - RMCN_K .* L + E;
W = RmG .* VmQ_squared - RMCN_K .* S + Y;


%% Calculating F3 ==========================================================
F = zeros(odeSystemSize, 1);
F(variableLocation(1,:)) = G .* N;
F(variableLocation(2,:)) = H;
F(variableLocation(3,:)) = W;

%% Calculating A3 ==========================================================
dC_dV = -C ./ VpVp2Vu;
dM_dV = (2 * kappa * V_u)^2 ./ VpVp2Vu_squared .* M;
dL_dQ = 1;
dL_dV = -(4 * kappa * V_u^2) ./ VpVp2Vu_squared;
LmQ = L - Q;
dS_dQ = 2 * ( LmQ .* dL_dQ - L );
dS_dV = (V_squared + Vp2Vu_squared) ./ VpVp2Vu_squared + 2 * LmQ .* dL_dV;
dE_dQ = ( 4 * Q_optQ - Q_opt_squared - 3 * V - 3 * Q_squared ) / (2 * V_s); 
dE_dV = ( 2 * Q_opt(:) - 3 * Q ) / (2 * V_s); 
dY_dQ = ( twoQ_optV - 6 * Q_optQ_squared + 2 * Q_opt_squaredQ + 4 * Q_cubed ) / (2 * V_s); 
dY_dV = ( 2 * Q_optQ - Q_opt_squared - 6 * V ) / (2 * V_s); 

dG_dN = -RMC_K;
dG_dQ = - QmQ_opt / V_s;
dG_dV = -RCN_K .* dM_dV - RMN_K .* dC_dV - 1 / (2 * V_s);

RN_KdMC_dV = RCN_K .* dM_dV + RMN_K .* dC_dV;

dH_dN = -dG_dN .* Q - RMC_K .* L;
dH_dQ = -dG_dQ .* Q + RmG - RMCN_K .* dL_dQ + dE_dQ;
dH_dV = -dG_dV .* Q - RMCN_K .* dL_dV - L .* RN_KdMC_dV + dE_dV;

dW_dN = -dG_dN .* VmQ_squared - RMC_K .* S;
dW_dQ = -dG_dQ .* VmQ_squared - 2 * Q .* RmG - RMCN_K .* dS_dQ + dY_dQ;
dW_dV = -dG_dV .* VmQ_squared + RmG - RMCN_K .* dS_dV - S .* RN_KdMC_dV + dY_dV;

rows = zeros(numNonZero, 1); % row indices of nonzero elements of A
columns = zeros(numNonZero, 1); % column indicis of nonzero elements of A
elements = zeros(numNonZero, 1); % nonzero elements of A
segmentCounter = 0;

%---constructing non-zero elements at 1st column of each block--------------------------------------
currentColumn = variableLocation(1,:);

segmentBegining = segmentCounter * segmentSize + 1;
segmentEnding = (segmentCounter + 1) * segmentSize;
rows(segmentBegining : segmentEnding) = variableLocation(1,:); % 1st row in each block
columns(segmentBegining : segmentEnding) = currentColumn; 
elements(segmentBegining : segmentEnding) = G + N .* dG_dN;
segmentCounter = segmentCounter + 1;

segmentBegining = segmentCounter * segmentSize + 1;
segmentEnding = (segmentCounter + 1) * segmentSize;
rows(segmentBegining : segmentEnding) = variableLocation(2,:); % 2nd row in each block
columns(segmentBegining : segmentEnding) = currentColumn; 
elements(segmentBegining : segmentEnding) = dH_dN;
segmentCounter = segmentCounter + 1;

segmentBegining = segmentCounter * segmentSize + 1;
segmentEnding = (segmentCounter + 1) * segmentSize;
rows(segmentBegining : segmentEnding) = variableLocation(3,:); % 3rd row in each block
columns(segmentBegining : segmentEnding) = currentColumn; 
elements(segmentBegining : segmentEnding) = dW_dN;
segmentCounter = segmentCounter + 1;

%---constructing non-zero elements at 2nd column of each block--------------------------------------
currentColumn = variableLocation(2,:);

segmentBegining = segmentCounter * segmentSize + 1;
segmentEnding = (segmentCounter + 1) * segmentSize;
rows(segmentBegining : segmentEnding) = variableLocation(1,:); % 1st row in each block
columns(segmentBegining : segmentEnding) = currentColumn; 
elements(segmentBegining : segmentEnding) = N .* dG_dQ;
segmentCounter = segmentCounter + 1;

segmentBegining = segmentCounter * segmentSize + 1;
segmentEnding = (segmentCounter + 1) * segmentSize;
rows(segmentBegining : segmentEnding) = variableLocation(2,:); % 2nd row in each block
columns(segmentBegining : segmentEnding) = currentColumn; 
elements(segmentBegining : segmentEnding) = dH_dQ;
segmentCounter = segmentCounter + 1;

segmentBegining = segmentCounter * segmentSize + 1;
segmentEnding = (segmentCounter + 1) * segmentSize;
rows(segmentBegining : segmentEnding) = variableLocation(3,:); % 3rd row in each block
columns(segmentBegining : segmentEnding) = currentColumn; 
elements(segmentBegining : segmentEnding) = dW_dQ;
segmentCounter = segmentCounter + 1;

%---constructing non-zero elements at 3rd column of each block--------------------------------------
currentColumn = variableLocation(3,:);

segmentBegining = segmentCounter * segmentSize + 1;
segmentEnding = (segmentCounter + 1) * segmentSize;
rows(segmentBegining : segmentEnding) = variableLocation(1,:); % 1st row in each block
columns(segmentBegining : segmentEnding) = currentColumn; 
elements(segmentBegining : segmentEnding) = N .* dG_dV;
segmentCounter = segmentCounter + 1;

segmentBegining = segmentCounter * segmentSize + 1;
segmentEnding = (segmentCounter + 1) * segmentSize;
rows(segmentBegining : segmentEnding) = variableLocation(2,:); % 2nd row in each block
columns(segmentBegining : segmentEnding) = currentColumn; 
elements(segmentBegining : segmentEnding) = dH_dV;
segmentCounter = segmentCounter + 1;

segmentBegining = segmentCounter * segmentSize + 1;
segmentEnding = (segmentCounter + 1) * segmentSize;
rows(segmentBegining : segmentEnding) = variableLocation(3,:); % 3rd row in each block
columns(segmentBegining : segmentEnding) = currentColumn; 
elements(segmentBegining : segmentEnding) = dW_dV;
segmentCounter = segmentCounter + 1;

%---constructing sparse A--------------------------------------------------------------------
A = sparse(rows, columns, elements, odeSystemSize, odeSystemSize);

end
