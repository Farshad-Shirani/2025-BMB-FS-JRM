function habitat = createPatchyHabitat(simulationParameters, discretizationParamaters)

habitat_min = 0.05; % minimum value of the habitat carrying capacity. maximum is 1.
L_c = 10; % Center patch length 
W_c = 10; % Center patch wodth 
fixedSize = 2; % fixed value used for width of patches on the right, length of patches on the left
variableSize = [10, 15]; % range of values for the variable length/width of patches

x1_0 = simulationParameters.x1_0;
x1_I = simulationParameters.x1_I;
x2_0 = simulationParameters.x2_0;
x2_J = simulationParameters.x2_J;

Dx1 = discretizationParamaters.Dx1;
Dx2 = discretizationParamaters.Dx2;

x1 = x1_0 : Dx1 : x1_I;
x2 = x2_0 : Dx2 : x2_J;
I = length(x1);
J = length(x2);

if strcmp( simulationParameters.boundaryCondition_x2, 'periodic' )
    J_hat = J-1; % if periodic boundary condition is used in x2-direction
else
    J_hat = J; % if reflecting boundary condition is used in x2-direction
end

%% Creating the Habitat ============================================================================
[X1, X2] = meshgrid(x1(1:I), x2(1:J_hat));
habitatLength = x1_I - x1_0;
habitatWidth = x2_J - x2_0;
habitat = zeros(size(X1));

%---Center Patch-------------------------------------------
habitat = habitat + makePatch( X1, X2, L_c, W_c);

%---Rigth segment------------------------------------------
W = fixedSize;
for i = 1 : ceil((habitatWidth/2) / W)
    j = 0; % to generate new seeds for random generator
    
    %---top segment---
    c2 = i * W - W/2; % x_2 coordinate of the center of the patch
    if abs(c2) < W_c / 2
        totalLength = L_c / 2; % not to override the center patch
    else
        totalLength = 0;
    end
    while 1
        rng(i + j,'twister');
        L = randi(variableSize); % random patch lengths
        c1 = L/2 + totalLength; % x_1 coordinate of the center of the patch

        habitat = habitat + makePatch( X1 - c1, X2 - c2,  L,  W );
        
        totalLength = totalLength + L;
        j = j+1;
        
        if totalLength > habitatLength / 2
            break;
        end
    end  
    
    %---bottom segment---
    c2 = W/2 - i * W; % x_2 coordinate of the center of the patch
    if abs(c2) < W_c / 2
        totalLength = L_c / 2; % not to override the center patch
    else
        totalLength = 0;
    end
    while 1
        rng(i + j,'twister');
        L = randi(variableSize); % random patch lengths
        c1 = L/2 + totalLength; % x_1 coordinate of the center of the patch

        habitat = habitat + makePatch( X1 - c1, X2 - c2,  L,  W );
        
        totalLength = totalLength + L;
        j = j+1;
        
        if totalLength > habitatLength / 2
            break;
        end
    end       

end


%---Left segment------------------------------------------
L = fixedSize;
for i = 1 : ceil((habitatLength/2) / L)
    j = 0; % to generate new seeds for random generator
    c1 = L/2 - i * L; % x_1 coordinate of the center of the patch
    
    if abs(c1) < L_c / 2
        totalWidth = W_c / 2; % not to override the center patch
        while 1 % constructing patches above the center
            rng(i + j,'twister');
            W = randi(variableSize); % random patch widths
            c2 = W/2 + totalWidth; % x_2 coordinate of the center of the patch 
            habitat = habitat + makePatch( X1 - c1, X2 - c2,  L,  W );
            totalWidth = totalWidth + W;
            j = j+1;
            if totalWidth > habitatWidth / 2
                break;
            end
        end
        totalWidth = W_c / 2; % not to override the center patch
        while 1 % constructing patches below the center
            rng(i + j,'twister');
            W = randi(variableSize); % random patch widths
            c2 = -W/2 - totalWidth; % x_2 coordinate of the center of the patch 
            habitat = habitat + makePatch( X1 - c1, X2 - c2,  L,  W );
            totalWidth = totalWidth + W;
            j = j+1;
            if totalWidth > habitatWidth / 2
                break;
            end
        end
    else
        totalWidth = 0; 
        while 1 % constructing patches right to the center
            rng(i + j,'twister');
            W = randi(variableSize); % random patch widths
            c2 = x2_0 + W/2 + totalWidth; % x_2 coordinate of the center of the patch 
            habitat = habitat + makePatch( X1 - c1, X2 - c2,  L,  W );
            totalWidth = totalWidth + W;
            j = j+1;
            if totalWidth > habitatWidth 
                break;
            end
        end 
    end    
end

%---Making bounadry values compatible with periodic boundary condition in x2-direction---------
if strcmp( simulationParameters.boundaryCondition_x2, 'periodic' )
    delta = 2; %width of the smooth cut-off function used at vicinity of the boundary
    chi = @(x, a) ( exp( (x-a).^2 ./ ((x-a).^2 - delta^2 - eps) ) .* (1 - heaviside(abs(x-a) - delta) ) ); % smooth cut-off function
    habitat = habitat - chi(X2, x2_0) .* habitat  - chi(X2, x2_J) .* habitat; % cutting the values of habitat at the vicinity of the x2-boundaries
end

%---Adjusting the minimum value-----------------------
habitat = habitat_min + (1 - habitat_min) *  habitat;

end




function patch = makePatch(X1, X2, L1, L2)
%---makes a patch, centered at the origin, that is tricubic in each direction (maximum is normalized to 1)
%---'patch' is a matrix giving values of the function (usually carrying capasity) 
%---L1 = patch length in x1-direction
%---L2 = patch length in x2-direction
%---minimum value of patch is 0 and maximum is 1

R1 = L1/2;
R2 = L2/2;
patch = (1 - abs(X1./R1).^3).^3 .* heaviside(1-abs(X1./R1)) ...
    .* (1 - abs(X2./R2).^3).^3 .* heaviside(1-abs(X2./R2)); % 0 <= patch <=1
end

