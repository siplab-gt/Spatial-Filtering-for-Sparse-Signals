function s_a = sangle(X, varargin)

% Calculate spectral angle
%
% s_a(x,y) = <x,y>/(||x||_2 * ||y||_2)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Parse inputs

if nargin > 1
    compare_mode = 'cross';
    Y = varargin{1};
    
    if size(X, 2) ~= size(Y,2)
        error('Number of vectors must be the same!')
    end
    
    if size(X, 1) ~= size(Y,1)
        error('Size of vectors must be the same!')
    end
    
else
    compare_mode = 'self';
    Y = [];
end


if nargin > 2
    comp_type = varargin{2};
else
    comp_type = 'vec';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Compute spectral angles

if strcmp(compare_mode, 'self')
%     num_vec = size(X, 2);
    norm_vals = sqrt(sum(X.^2, 1));
    inner_prod = (X')*X;
    s_a = acosd(inner_prod./((norm_vals')*norm_vals));
elseif strcmp(compare_mode, 'cross')
    norm_X = sqrt(sum(X.^2, 1));
    norm_Y = sqrt(sum(Y.^2, 1));
    if strcmp(comp_type, 'vec')
        inner_prod = sum((Y).*X, 1);
        s_a = acosd(inner_prod./((norm_X).*norm_Y));
    else
        inner_prod = (Y')*X;
        s_a = acosd(inner_prod./((norm_X')*norm_Y));
    end
else
    error('Bad comparison mode!')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
