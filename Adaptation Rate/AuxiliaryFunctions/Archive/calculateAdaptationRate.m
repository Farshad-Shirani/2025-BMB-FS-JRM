function [random,  nonRandom] = calculateAdaptationRate(solution, variableLocation, Dx)
%--- This function is the same as construct_F1. It just returns the two components of F1 associated
%--- with random and non-random gene flow separaterly. F1 is the some of the two components. 


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

random = d2S(:);
random(variableLocation(1,:)) = random(variableLocation(1,:)) * sigma_squared;
random(variableLocation(2,:)) = ( random(variableLocation(2,:)) + 2 * dN_N .* dQ ) * sigma_squared;
random(variableLocation(3,:)) = ( random(variableLocation(3,:)) + 2 * dN_N .* dV + 2 * dQ.^2 ) * sigma_squared;

nonRandom = zeros(size(random));
nonRandom(variableLocation(1,:)) = - dANQd;
nonRandom(variableLocation(2,:)) = - dAVd - dN_N .* AVd - dQ .* AQd;
nonRandom(variableLocation(3,:)) = - 2 * dQ .* AVd - dV .* AQd;

end