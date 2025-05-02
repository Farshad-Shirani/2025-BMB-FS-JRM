function F = construct_F1(solution, variableLocation, Dx)
%--- This function returnes the value of the F1 component of ODE vector field. The ordering is in x1 direction.

global D V_u Q_opt A_max dQ_tilde
sigma_squared = D;

Ad = (A_max ./ V_u) .* dQ_tilde;
AQd = Ad .* ( solution(:,2) - Q_opt );
AVd = Ad .* solution(:,3);
%--------------------------------------
solution_minusTwo = reflect(solution,[2 0]); % solution(i-2)  reflect also applies the reflecting boundary condition
solution_minusOne = reflect(solution,[1 0]); % solution(i-1)  
solution_plusOne = reflect(solution,[-1 0]); % solution(i+1)  
solution_plusTwo = reflect(solution,[-2 0]); % solution(i+2)

Ad_minusTwo = reflect(Ad,[2 0]); % Ad(i-2)  reflect also applies the reflecting boundary condition
Ad_minusOne = reflect(Ad,[1 0]); % Ad(i-1)  
Ad_plusOne = reflect(Ad,[-1 0]); % Ad(i+1)  
Ad_plusTwo = reflect(Ad,[-2 0]); % Ad(i+2)

AQd_minusTwo = reflect(AQd,[2 0]); % AQd(i-2)  reflect also applies the reflecting boundary condition
AQd_minusOne = reflect(AQd,[1 0]); % AQd(i-1)  
AQd_plusOne = reflect(AQd,[-1 0]); % AQd(i+1)  
AQd_plusTwo = reflect(AQd,[-2 0]); % AQd(i+2)
%----------------------------------------

dS = ( -solution_plusTwo + 8 * solution_plusOne - 8 * solution_minusOne + solution_minusTwo ) / (12 * Dx);
dN_N = dS(:,1) ./ solution(:,1);
dQ = dS(:,2);
dV = dS(:,3);

d2S = ( -solution_plusTwo + 16 * solution_plusOne - 30 * solution + 16 * solution_minusOne - solution_minusTwo ) / (12 * Dx^2);
d2S = d2S'; % reorders 2-dim array d2S so that linear indexing d2S(:) corresponds to current ordeing in U

AVd_minusTwo = Ad_minusTwo .* solution_minusTwo(:,3);
AVd_minusOne = Ad_minusOne .* solution_minusOne(:,3);
AVd_plusOne = Ad_plusOne .* solution_plusOne(:,3);
AVd_plusTwo = Ad_plusTwo .* solution_plusTwo(:,3);
dAVd = ( -AVd_plusTwo + 8 * AVd_plusOne - 8 * AVd_minusOne + AVd_minusTwo ) / (12 * Dx);

ANQd_minusTwo = AQd_minusTwo .* solution_minusTwo(:,1);
ANQd_minusOne = AQd_minusOne .* solution_minusOne(:,1);
ANQd_plusOne = AQd_plusOne .* solution_plusOne(:,1);
ANQd_plusTwo = AQd_plusTwo .* solution_plusTwo(:,1);
dANQd = ( -ANQd_plusTwo + 8 * ANQd_plusOne - 8 * ANQd_minusOne + ANQd_minusTwo ) / (12 * Dx);

F = d2S(:);
F(variableLocation(1,:)) = F(variableLocation(1,:)) * sigma_squared   - dANQd;
F(variableLocation(2,:)) = ( F(variableLocation(2,:)) + 2 * dN_N .* dQ ) * sigma_squared ...
    - dAVd - dN_N .* AVd - dQ .* AQd;
F(variableLocation(3,:)) = ( F(variableLocation(3,:)) + 2 * dN_N .* dV + 2 * dQ.^2 ) * sigma_squared ...
    - 2 * dQ .* AVd - dV .* AQd;

end