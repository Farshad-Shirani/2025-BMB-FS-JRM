function Y = reflect(X, varargin)
%   My REFLECT reflects positions of elements on the boundary
%   Y = REFLECT(X,K) where K is an integer scalar reflects 
%   the elements in the array X by K positions. If X is a vector and K is
%   positive, then the values of X are shifted from the beginning 
%   to the end and positions released at the begining are filled by the mirror of the entries at the beginig.
%   If K is negative, they are shifted from the end to the 
%   beginning and positions released at the end are filled by the mirror of the entries at the end. 
%   If X is a matrix, REFLECT reflects along columns. If X is an
%   N-D array, REFLECT refelcts along the first nonsingleton dimension.
%   
%   Y = REFLECT(X,K,DIM) reflects along the dimension DIM.
%
%   Y = REFLECT(X,V) REFLECTS the values in the array X
%   by V elements. V is a vector of integers where the N-th element 
%   specifies the shift amount along the N-th dimension of
%   array X.

reflections = zeros(1, ndims(X));
if nargin == 2
    if length(varargin{1}) == 1
        ind = find(size(X) > 1); % Y = REFLECT(X,K)
        reflections(ind(1)) = varargin{1}; 
    else
        if length(varargin{1}) ~= ndims(X)
            disp('Error in REFLECT: Size of the vector V must match the dimension of the array X!')
            Y = NaN;    
            return
        end            
        reflections = varargin{1}; % Y = REFLECT(X,V)
    end       
elseif nargin == 3
    reflections(varargin{2}) = varargin{1}; %Y = REFLECT(X,K,DIM)
else
    disp('Error in REFLECT: Invalid number of input arguments!')
    Y = NaN;
    return
end


for i = 1 : length(reflections)
    if abs(reflections(i)) >= size(X,i)
        disp(['Error in REFLECT: Reflection size is larger than array dimension ', num2str(i)])
        Y = NaN;
        return
    end
    colons_begining = repmat({':'}, 1, i-1);
    colons_end = repmat({':'}, 1, ndims(X)-i);
    last = size(X,i);
    if reflections(i) > 0
        k = reflections(i);
        mirror = X(colons_begining{:}, k+1:-1:2, colons_end{:});
        X(colons_begining{:}, 1+k : last, colons_end{:}) = X(colons_begining{:}, 1 : last-k, colons_end{:}); %shift
        X(colons_begining{:}, 1 : k, colons_end{:}) = mirror; %fill by mirror       
    elseif reflections(i) < 0
        k = -reflections(i);
        mirror = X(colons_begining{:}, last-1:-1:last-k, colons_end{:});
        X(colons_begining{:}, 1 : last-k, colons_end{:}) = X(colons_begining{:}, 1+k : last, colons_end{:}); %shift
        X(colons_begining{:}, last-k+1:last, colons_end{:}) = mirror; %fill by mirror 
    end  
end 

Y = X;
