function F = construct_F2(solution, variableLocation, Dx)
%--- This function returns the value of the F2 component of ODE vector field. The ordering is in x1 direction.

global D V_u Q_opt A_max dQ_tilde boundaryCondition_x2
sigma_squared = D(2,2);
dQ_tilde_here = dQ_tilde(:,:,2); %2nd component of dQ_tilde, since we are calulating F2

Ad = (A_max ./ V_u) .* dQ_tilde_here;
AQd = Ad .* ( solution(:,:,2) - Q_opt );
AVd = Ad .* solution(:,:,3);

if strcmp( boundaryCondition_x2, 'periodic' )
    solution_minusTwo = circshift( solution, [0 2 0] ); % solution(i,j-2)  circshift also applies the periodic boundary condition
    solution_minusOne = circshift( solution, [0 1 0] ); % solution(i,j-1)  
    solution_plusOne = circshift( solution, [0 -1 0] ); % solution(i,j+1)  
    solution_plusTwo = circshift( solution, [0 -2 0] ); % solution(i,j+2)

    Ad_minusTwo = circshift( Ad, [0 2] ); % Ad(i,j-2)  circshift also applies the periodic boundary condition
    Ad_minusOne = circshift( Ad, [0 1] ); % Ad(i,j-1)  
    Ad_plusOne = circshift( Ad, [0 -1] ); % Ad(i,j+1)  
    Ad_plusTwo = circshift( Ad, [0 -2] ); % Ad(i,j+2)

    AQd_minusTwo = circshift( AQd, [0 2] ); % AQd(i,j-2)  circshift also applies the periodic boundary condition
    AQd_minusOne = circshift( AQd, [0 1] ); % AQd(i,j-1)  
    AQd_plusOne = circshift( AQd, [0 -1] ); % AQd(i,j+1)  
    AQd_plusTwo = circshift( AQd, [0 -2] ); % AQd(i,j+2)
else
    solution_minusTwo = reflect( solution, [0 2 0] ); % solution(i,j-2)  reflect also applies the reflecting boundary condition
    solution_minusOne = reflect( solution, [0 1 0] ); % solution(i,j-1)  
    solution_plusOne = reflect( solution, [0 -1 0] ); % solution(i,j+1)  
    solution_plusTwo = reflect( solution, [0 -2 0] ); % solution(i,j+2)
 
    Ad_minusTwo = reflect( Ad, [0 2] ); % Ad(i,j-2)  reflect also applies the periodic boundary condition
    Ad_minusOne = reflect( Ad, [0 1] ); % Ad(i,j-1)  
    Ad_plusOne = reflect( Ad, [0 -1] ); % Ad(i,j+1)  
    Ad_plusTwo = reflect( Ad, [0 -2] ); % Ad(i,j+2)
 
    AQd_minusTwo = reflect( AQd, [0 2] ); % AQd(i,j-2)  refelect also applies the periodic boundary condition
    AQd_minusOne = reflect( AQd, [0 1] ); % AQd(i,j-1)  
    AQd_plusOne = reflect( AQd, [0 -1] ); % AQd(i,j+1)  
    AQd_plusTwo = reflect( AQd, [0 -2] ); % AQd(i,j+2)
end
%------------------------------------

dS = ( -solution_plusTwo + 8*solution_plusOne - 8*solution_minusOne + solution_minusTwo ) / (12 * Dx);
dN_N = dS(:,:,1) ./ solution(:,:,1);
dQ = dS(:,:,2);
dV = dS(:,:,3);

d2S = ( -solution_plusTwo + 16*solution_plusOne - 30*solution + 16*solution_minusOne - solution_minusTwo ) / (12 * Dx^2);
d2S = permute(d2S,[3,1,2]); % reorders 3-dim array d2S so that linear indexing d2S(:) corresponds to current ordeing in U

AVd_minusTwo = Ad_minusTwo .* solution_minusTwo(:,:,3);
AVd_minusOne = Ad_minusOne .* solution_minusOne(:,:,3);
AVd_plusOne = Ad_plusOne .* solution_plusOne(:,:,3);
AVd_plusTwo = Ad_plusTwo .* solution_plusTwo(:,:,3);
dAVd = ( -AVd_plusTwo + 8 * AVd_plusOne - 8 * AVd_minusOne + AVd_minusTwo ) / (12 * Dx);

ANQd_minusTwo = AQd_minusTwo .* solution_minusTwo(:,:,1);
ANQd_minusOne = AQd_minusOne .* solution_minusOne(:,:,1);
ANQd_plusOne = AQd_plusOne .* solution_plusOne(:,:,1);
ANQd_plusTwo = AQd_plusTwo .* solution_plusTwo(:,:,1);
dANQd = ( -ANQd_plusTwo + 8 * ANQd_plusOne - 8 * ANQd_minusOne + ANQd_minusTwo ) / (12 * Dx);

dN_N = dN_N(:);
dQ = dQ(:);
dV = dV(:);
AQd = AQd(:);
AVd = AVd(:);
dAVd = dAVd(:);
dANQd = dANQd(:);

F = d2S(:);
F(variableLocation(1,:)) = F(variableLocation(1,:)) * sigma_squared   - dANQd;
F(variableLocation(2,:)) = ( F(variableLocation(2,:)) + 2 * dN_N .* dQ ) * sigma_squared ...
    - dAVd - dN_N .* AVd - dQ .* AQd;
F(variableLocation(3,:)) = ( F(variableLocation(3,:)) + 2 * dN_N .* dV + 2 * dQ.^2 ) * sigma_squared ...
    - 2 * dQ .* AVd - dV .* AQd;

end