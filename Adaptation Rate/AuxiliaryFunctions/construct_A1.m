function A = construct_A1(solution, variableLocation, Dx)
%--- This function gives the sparse matrix resulted from discritization in x1 direction & linearization

global D V_u Q_opt A_max dQ_tilde
sigma_squared = D;

segmentSize = numel(solution(:,1)); % number of (A)ij segments
numNonZero = 1 * (9 + 15 + 15) * segmentSize; %number of non-zero elements in A.
odeSystemSize = numel(solution);

%-----------------------------------------------
solution_minusTwo = reflect(solution,[2 0]); % solution(i-2)  reflect also applies the reflecting boundary condition
solution_minusOne = reflect(solution,[1 0]); % solution(i-1)  
solution_plusOne = reflect(solution,[-1 0]); % solution(i+1)  
solution_plusTwo = reflect(solution,[-2 0]); % solution(i+2)  
%----------------------------------------------

dS = ( -solution_plusTwo + 8*solution_plusOne - 8*solution_minusOne + solution_minusTwo ) / (12 * Dx);

%---creating linear arrays
N = solution(:,1);
Q = solution(:,2);
V = solution(:,3);
dN = dS(:,1);
dQ = dS(:,2);
dV = dS(:,3);
dN_N = dN ./ N;
dQ_N = dQ ./ N;
dV_N = dV./ N;

Ad = (A_max ./ V_u) .* dQ_tilde;
ANd = Ad .* N;
AVd = Ad .* V;
AQd = Ad .* (Q - Q_opt);

AVd_N = AVd ./ N;
AddN_N = Ad .* dN_N;
AddQ = Ad .* dQ;
AddV = Ad .* dV;


rows = zeros(numNonZero, 1); % row indices of nonzero elements of A
columns = zeros(numNonZero, 1); % column indicis of nonzero elements of A
elements = zeros(numNonZero, 1); % nonzero elements of A

segmentCounter = 0;
%% diagonal blocks of A (center block of (A)ij)===========================================================
firstColumn = variableLocation(1,:); % indices of the first (left) columns of the center block of (A)ij. 

segmentBegining = segmentCounter * segmentSize + 1;
segmentEnding = (segmentCounter + 1) * segmentSize;
rows(segmentBegining : segmentEnding) = variableLocation(1,:); % 1st row in each block
columns(segmentBegining : segmentEnding) = firstColumn; % 1st column in each block
elements(segmentBegining : segmentEnding) = -30 / (12*Dx^2) * sigma_squared;
segmentCounter = segmentCounter + 1;

%------------------

segmentBegining = segmentCounter * segmentSize + 1;
segmentEnding = (segmentCounter + 1) * segmentSize;
rows(segmentBegining : segmentEnding) = variableLocation(2,:); % 2nd row in each block
columns(segmentBegining : segmentEnding) = firstColumn; % 1st column in each block
elements(segmentBegining : segmentEnding) = -2 * dN_N .* dQ_N * sigma_squared   + dN_N .* AVd_N;
segmentCounter = segmentCounter + 1;

segmentBegining = segmentCounter * segmentSize + 1;
segmentEnding = (segmentCounter + 1) * segmentSize;
rows(segmentBegining : segmentEnding) = variableLocation(2,:); % 2nd row in each block
columns(segmentBegining : segmentEnding) = firstColumn + 1; % 2nd column in each block
elements(segmentBegining : segmentEnding) = -30 / (12*Dx^2) * sigma_squared   - AddQ;
segmentCounter = segmentCounter + 1;

segmentBegining = segmentCounter * segmentSize + 1;
segmentEnding = (segmentCounter + 1) * segmentSize;
rows(segmentBegining : segmentEnding) = variableLocation(2,:); % 2nd row in each block
columns(segmentBegining : segmentEnding) = firstColumn + 2; % 3rd column in each block
elements(segmentBegining : segmentEnding) = - AddN_N;
segmentCounter = segmentCounter + 1;

%-------------------

segmentBegining = segmentCounter * segmentSize + 1;
segmentEnding = (segmentCounter + 1) * segmentSize;
rows(segmentBegining : segmentEnding) = variableLocation(3,:); % 3rd row in each block
columns(segmentBegining : segmentEnding) = firstColumn; % 1st column in each block
elements(segmentBegining : segmentEnding) = -2 * dN_N .* dV_N * sigma_squared;
segmentCounter = segmentCounter + 1;

segmentBegining = segmentCounter * segmentSize + 1;
segmentEnding = (segmentCounter + 1) * segmentSize;
rows(segmentBegining : segmentEnding) = variableLocation(3,:); % 3rd row in each block
columns(segmentBegining : segmentEnding) = firstColumn + 1; % 2nd column in each block
elements(segmentBegining : segmentEnding) = - AddV;
segmentCounter = segmentCounter + 1;

segmentBegining = segmentCounter * segmentSize + 1;
segmentEnding = (segmentCounter + 1) * segmentSize;
rows(segmentBegining : segmentEnding) = variableLocation(3,:); % 3rd row in each block
columns(segmentBegining : segmentEnding) = firstColumn + 2; % 3rd column in each block
elements(segmentBegining : segmentEnding) = -30 / (12*Dx^2) * sigma_squared   - 2 * AddQ;
segmentCounter = segmentCounter + 1;


%% 1st lower diagonal of A (first left block of (A)ij)====================================================
firstColumn = reflect(variableLocation(1,:), 1); % indices of the first (left) columns of the first left block of (A)ij. reflect applies the reflecting boundary condition

segmentBegining = segmentCounter * segmentSize + 1;
segmentEnding = (segmentCounter + 1) * segmentSize;
rows(segmentBegining : segmentEnding) = variableLocation(1,:); % 1st row in each block
columns(segmentBegining : segmentEnding) = firstColumn; % 1st column in each block
elements(segmentBegining : segmentEnding) = 16 / (12*Dx^2) * sigma_squared   + 8 / (12*Dx) * AQd;
segmentCounter = segmentCounter + 1;

segmentBegining = segmentCounter * segmentSize + 1;
segmentEnding = (segmentCounter + 1) * segmentSize;
rows(segmentBegining : segmentEnding) = variableLocation(1,:); % 1st row in each block
columns(segmentBegining : segmentEnding) = firstColumn + 1; % 2nd column in each block
elements(segmentBegining : segmentEnding) = 8 / (12*Dx) * ANd;
segmentCounter = segmentCounter + 1;

%----------------------

segmentBegining = segmentCounter * segmentSize + 1;
segmentEnding = (segmentCounter + 1) * segmentSize;
rows(segmentBegining : segmentEnding) = variableLocation(2,:); % 2nd row in each block
columns(segmentBegining : segmentEnding) = firstColumn; % 1st column in each block
elements(segmentBegining : segmentEnding) = 2 * (-8) / (12*Dx) * dQ_N * sigma_squared   + 8 / (12*Dx) .* AVd_N;
segmentCounter = segmentCounter + 1;

segmentBegining = segmentCounter * segmentSize + 1;
segmentEnding = (segmentCounter + 1) * segmentSize;
rows(segmentBegining : segmentEnding) = variableLocation(2,:); % 2nd row in each block
columns(segmentBegining : segmentEnding) = firstColumn + 1; % 2nd column in each block
elements(segmentBegining : segmentEnding) = ( 16 / (12*Dx^2) + 2 * (-8) / (12*Dx) * dN_N ) * sigma_squared ...
    + 8 / (12*Dx) * AQd - AddQ;
segmentCounter = segmentCounter + 1;

segmentBegining = segmentCounter * segmentSize + 1;
segmentEnding = (segmentCounter + 1) * segmentSize;
rows(segmentBegining : segmentEnding) = variableLocation(2,:); % 2nd row in each block
columns(segmentBegining : segmentEnding) = firstColumn + 2; % 3rd column in each block
elements(segmentBegining : segmentEnding) = 8 / (12*Dx) * Ad - AddN_N;
segmentCounter = segmentCounter + 1;

%------------------------------------

segmentBegining = segmentCounter * segmentSize + 1;
segmentEnding = (segmentCounter + 1) * segmentSize;
rows(segmentBegining : segmentEnding) = variableLocation(3,:); % 3rd row in each block
columns(segmentBegining : segmentEnding) = firstColumn; % 1st column in each block
elements(segmentBegining : segmentEnding) = 2 * (-8) / (12*Dx) * dV_N * sigma_squared;
segmentCounter = segmentCounter + 1;

segmentBegining = segmentCounter * segmentSize + 1;
segmentEnding = (segmentCounter + 1) * segmentSize;
rows(segmentBegining : segmentEnding) = variableLocation(3,:); % 3rd row in each block
columns(segmentBegining : segmentEnding) = firstColumn + 1; % 2nd column in each block
elements(segmentBegining : segmentEnding) = 4 * (-8) / (12*Dx) * dQ * sigma_squared ...
    + 2 * 8 / (12*Dx) * AVd - AddV;
segmentCounter = segmentCounter + 1;

segmentBegining = segmentCounter * segmentSize + 1;
segmentEnding = (segmentCounter + 1) * segmentSize;
rows(segmentBegining : segmentEnding) = variableLocation(3,:); % 3rd row in each block
columns(segmentBegining : segmentEnding) = firstColumn + 2; % 3rd column in each block
elements(segmentBegining : segmentEnding) = ( 16 / (12*Dx^2) + 2 * (-8) / (12*Dx) * dN_N ) * sigma_squared ...
    -2 * AddQ + 8 / (12*Dx) * AQd;
segmentCounter = segmentCounter + 1;


%% 2nd lower diagonal of A (secend left block of (A)ij)======================================================
firstColumn = reflect(variableLocation(1,:), 2); % indices of the first (left) columns of the second left block of (A)ij. 

segmentBegining = segmentCounter * segmentSize + 1;
segmentEnding = (segmentCounter + 1) * segmentSize;
rows(segmentBegining : segmentEnding) = variableLocation(1,:); % 1st row in each block
columns(segmentBegining : segmentEnding) = firstColumn; % 1st column in each block
elements(segmentBegining : segmentEnding) = -1 / (12*Dx^2) * sigma_squared   - 1 / (12*Dx) * AQd;
segmentCounter = segmentCounter + 1;

segmentBegining = segmentCounter * segmentSize + 1;
segmentEnding = (segmentCounter + 1) * segmentSize;
rows(segmentBegining : segmentEnding) = variableLocation(1,:); % 1st row in each block
columns(segmentBegining : segmentEnding) = firstColumn + 1; % 2nd column in each block
elements(segmentBegining : segmentEnding) = -1 / (12*Dx) * ANd;
segmentCounter = segmentCounter + 1;

%---------------------

segmentBegining = segmentCounter * segmentSize + 1;
segmentEnding = (segmentCounter + 1) * segmentSize;
rows(segmentBegining : segmentEnding) = variableLocation(2,:); % 2nd row in each block
columns(segmentBegining : segmentEnding) = firstColumn; % 1st column in each block
elements(segmentBegining : segmentEnding) = 2 / (12*Dx) * dQ_N * sigma_squared   - 1 / (12*Dx) * AVd_N;
segmentCounter = segmentCounter + 1;

segmentBegining = segmentCounter * segmentSize + 1;
segmentEnding = (segmentCounter + 1) * segmentSize;
rows(segmentBegining : segmentEnding) = variableLocation(2,:); % 2nd row in each block
columns(segmentBegining : segmentEnding) = firstColumn + 1; % 2nd column in each block
elements(segmentBegining : segmentEnding) = ( -1 / (12*Dx^2) + 2 / (12*Dx) * dN_N ) * sigma_squared ...
    - 1 / (12*Dx) * AQd - AddQ;
segmentCounter = segmentCounter + 1;

segmentBegining = segmentCounter * segmentSize + 1;
segmentEnding = (segmentCounter + 1) * segmentSize;
rows(segmentBegining : segmentEnding) = variableLocation(2,:); % 2nd row in each block
columns(segmentBegining : segmentEnding) = firstColumn + 2; % 3rd column in each block
elements(segmentBegining : segmentEnding) = - 1 / (12*Dx) * Ad - AddN_N;
segmentCounter = segmentCounter + 1;

%---------------------

segmentBegining = segmentCounter * segmentSize + 1;
segmentEnding = (segmentCounter + 1) * segmentSize;
rows(segmentBegining : segmentEnding) = variableLocation(3,:); % 3rd row in each block
columns(segmentBegining : segmentEnding) = firstColumn; % 1st column in each block
elements(segmentBegining : segmentEnding) = 2 / (12*Dx) * dV_N * sigma_squared;
segmentCounter = segmentCounter + 1;

segmentBegining = segmentCounter * segmentSize + 1;
segmentEnding = (segmentCounter + 1) * segmentSize;
rows(segmentBegining : segmentEnding) = variableLocation(3,:); % 3rd row in each block
columns(segmentBegining : segmentEnding) = firstColumn + 1; % 2nd column in each block
elements(segmentBegining : segmentEnding) = 4 / (12*Dx) * dQ * sigma_squared   - 2 / (12*Dx)* AVd - AddV;
segmentCounter = segmentCounter + 1;

segmentBegining = segmentCounter * segmentSize + 1;
segmentEnding = (segmentCounter + 1) * segmentSize;
rows(segmentBegining : segmentEnding) = variableLocation(3,:); % 3rd row in each block
columns(segmentBegining : segmentEnding) = firstColumn + 2; % 3rd column in each block
elements(segmentBegining : segmentEnding) = ( -1 / (12*Dx^2) + 2 / (12*Dx) * dN_N ) * sigma_squared ...
    - 2 * AddQ - 1 / (12*Dx) * AQd;
segmentCounter = segmentCounter + 1;


%% 1st upper diagonal of A (first right block of (A)ij)======================================================
firstColumn = reflect(variableLocation(1,:), -1); % indices of the first (left) columns of the first right block of (A)ij. 

segmentBegining = segmentCounter * segmentSize + 1;
segmentEnding = (segmentCounter + 1) * segmentSize;
rows(segmentBegining : segmentEnding) = variableLocation(1,:); % 1st row in each block
columns(segmentBegining : segmentEnding) = firstColumn; % 1st column in each block
elements(segmentBegining : segmentEnding) = 16 / (12*Dx^2) * sigma_squared    - 8 / (12*Dx) * AQd;
segmentCounter = segmentCounter + 1;

segmentBegining = segmentCounter * segmentSize + 1;
segmentEnding = (segmentCounter + 1) * segmentSize;
rows(segmentBegining : segmentEnding) = variableLocation(1,:); % 1st row in each block
columns(segmentBegining : segmentEnding) = firstColumn + 1; % 2nd column in each block
elements(segmentBegining : segmentEnding) = -8 / (12*Dx) * ANd;
segmentCounter = segmentCounter + 1;

%---------------------------

segmentBegining = segmentCounter * segmentSize + 1;
segmentEnding = (segmentCounter + 1) * segmentSize;
rows(segmentBegining : segmentEnding) = variableLocation(2,:); % 2nd row in each block
columns(segmentBegining : segmentEnding) = firstColumn; % 1st column in each block
elements(segmentBegining : segmentEnding) = 2 * 8 / (12*Dx) * dQ_N * sigma_squared   - 8 / (12*Dx) .* AVd_N;
segmentCounter = segmentCounter + 1;

segmentBegining = segmentCounter * segmentSize + 1;
segmentEnding = (segmentCounter + 1) * segmentSize;
rows(segmentBegining : segmentEnding) = variableLocation(2,:); % 2nd row in each block
columns(segmentBegining : segmentEnding) = firstColumn + 1; % 2nd column in each block
elements(segmentBegining : segmentEnding) = ( 16 / (12*Dx^2) + 2 * 8 / (12*Dx) * dN_N ) * sigma_squared ...
    - 8 / (12*Dx) * AQd - AddQ;
segmentCounter = segmentCounter + 1;

segmentBegining = segmentCounter * segmentSize + 1;
segmentEnding = (segmentCounter + 1) * segmentSize;
rows(segmentBegining : segmentEnding) = variableLocation(2,:); % 2nd row in each block
columns(segmentBegining : segmentEnding) = firstColumn + 2; % 3rd column in each block
elements(segmentBegining : segmentEnding) = -8 / (12*Dx) * Ad - AddN_N;
segmentCounter = segmentCounter + 1;

%----------------------------------

segmentBegining = segmentCounter * segmentSize + 1;
segmentEnding = (segmentCounter + 1) * segmentSize;
rows(segmentBegining : segmentEnding) = variableLocation(3,:); % 3rd row in each block
columns(segmentBegining : segmentEnding) = firstColumn; % 1st column in each block
elements(segmentBegining : segmentEnding) = 2 * 8 / (12*Dx) * dV_N * sigma_squared;
segmentCounter = segmentCounter + 1;

segmentBegining = segmentCounter * segmentSize + 1;
segmentEnding = (segmentCounter + 1) * segmentSize;
rows(segmentBegining : segmentEnding) = variableLocation(3,:); % 3rd row in each block
columns(segmentBegining : segmentEnding) = firstColumn + 1; % 2nd column in each block
elements(segmentBegining : segmentEnding) = 4 * 8 / (12*Dx) * dQ * sigma_squared ...
    - 2 * 8 / (12*Dx) * AVd - AddV;
segmentCounter = segmentCounter + 1;

segmentBegining = segmentCounter * segmentSize + 1;
segmentEnding = (segmentCounter + 1) * segmentSize;
rows(segmentBegining : segmentEnding) = variableLocation(3,:); % 3rd row in each block
columns(segmentBegining : segmentEnding) = firstColumn + 2; % 3rd column in each block
elements(segmentBegining : segmentEnding) = ( 16 / (12*Dx^2) + 2 * 8 / (12*Dx) * dN_N ) * sigma_squared ...
    -2 * AddQ - 8 / (12*Dx) * AQd;
segmentCounter = segmentCounter + 1;


%% 2nd upper diagonal of A (secend right block of (A)ij)=====================================================
firstColumn = reflect(variableLocation(1,:), -2); % indices of the first (left) columns of the second right block of (A)ij. 

segmentBegining = segmentCounter * segmentSize + 1;
segmentEnding = (segmentCounter + 1) * segmentSize;
rows(segmentBegining : segmentEnding) = variableLocation(1,:); % 1st row in each block
columns(segmentBegining : segmentEnding) = firstColumn; % 1st column in each block
elements(segmentBegining : segmentEnding) = -1 / (12*Dx^2) * sigma_squared   + 1 / (12*Dx) * AQd;
segmentCounter = segmentCounter + 1;

segmentBegining = segmentCounter * segmentSize + 1;
segmentEnding = (segmentCounter + 1) * segmentSize;
rows(segmentBegining : segmentEnding) = variableLocation(1,:); % 1st row in each block
columns(segmentBegining : segmentEnding) = firstColumn + 1; % 2nd column in each block
elements(segmentBegining : segmentEnding) = 1 / (12*Dx) * ANd;
segmentCounter = segmentCounter + 1;

%---------------------

segmentBegining = segmentCounter * segmentSize + 1;
segmentEnding = (segmentCounter + 1) * segmentSize;
rows(segmentBegining : segmentEnding) = variableLocation(2,:); % 2nd row in each block
columns(segmentBegining : segmentEnding) = firstColumn; % 1st column in each block
elements(segmentBegining : segmentEnding) = 2 * (-1) / (12*Dx) * dQ_N * sigma_squared   + 1 / (12*Dx) * AVd_N;
segmentCounter = segmentCounter + 1;

segmentBegining = segmentCounter * segmentSize + 1;
segmentEnding = (segmentCounter + 1) * segmentSize;
rows(segmentBegining : segmentEnding) = variableLocation(2,:); % 2nd row in each block
columns(segmentBegining : segmentEnding) = firstColumn + 1; % 2nd column in each block
elements(segmentBegining : segmentEnding) = ( -1 / (12*Dx^2) + 2 * (-1) / (12*Dx) * dN_N ) * sigma_squared ...
    + 1 / (12*Dx) * AQd - AddQ;
segmentCounter = segmentCounter + 1;

segmentBegining = segmentCounter * segmentSize + 1;
segmentEnding = (segmentCounter + 1) * segmentSize;
rows(segmentBegining : segmentEnding) = variableLocation(2,:); % 2nd row in each block
columns(segmentBegining : segmentEnding) = firstColumn + 2; % 3rd column in each block
elements(segmentBegining : segmentEnding) = 1 / (12*Dx) * Ad - AddN_N;
segmentCounter = segmentCounter + 1;

%-----------------------

segmentBegining = segmentCounter * segmentSize + 1;
segmentEnding = (segmentCounter + 1) * segmentSize;
rows(segmentBegining : segmentEnding) = variableLocation(3,:); % 3rd row in each block
columns(segmentBegining : segmentEnding) = firstColumn; % 1st column in each block
elements(segmentBegining : segmentEnding) = 2 * (-1) / (12*Dx) * dV_N * sigma_squared;
segmentCounter = segmentCounter + 1;

segmentBegining = segmentCounter * segmentSize + 1;
segmentEnding = (segmentCounter + 1) * segmentSize;
rows(segmentBegining : segmentEnding) = variableLocation(3,:); % 3rd row in each block
columns(segmentBegining : segmentEnding) = firstColumn + 1; % 2nd column in each block
elements(segmentBegining : segmentEnding) = 4 * (-1) / (12*Dx) * dQ * sigma_squared   + 2 / (12*Dx)* AVd - AddV;
segmentCounter = segmentCounter + 1;

segmentBegining = segmentCounter * segmentSize + 1;
segmentEnding = (segmentCounter + 1) * segmentSize;
rows(segmentBegining : segmentEnding) = variableLocation(3,:); % 3rd row in each block
columns(segmentBegining : segmentEnding) = firstColumn + 2; % 3rd column in each block
elements(segmentBegining : segmentEnding) = ( -1 / (12*Dx^2) + 2 * (-1) / (12*Dx) * dN_N ) * sigma_squared ...
    - 2 * AddQ + 1 / (12*Dx) * AQd;
segmentCounter = segmentCounter + 1;


%% Constructing sparse A ==================================================================================

A = sparse(rows, columns, elements, odeSystemSize, odeSystemSize);

end

