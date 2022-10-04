function H_temp = HS_image(index_struct, data, HS_size)

% H_temp = function HS_image(index_struct, data)
% 
% Generates a Hyperspectral Image true to the scene.
% index_struct.x - x indecies in order of where they go
% index_struct.y - y indecies in order of where they go
% data - vector of same length as index.x and index.y
% 
% 
% Last Modified 3/16/2010 - Adam Charles

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Error Checking

if (size(data, 1) ~= 1) && (size(data, 2) ~= 1)
    error('data must be a vector!')
end

if size(data, 1) == length(index_struct.x)
    data = data.';
end

if (length(index_struct.x) ~= size(data, 2)) || (length(index_struct.y) ~= size(data, 2))
    error('Index numbers don''t match!')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot

if nargin < 3
    H_temp = zeros(max(index_struct.x), max(index_struct.y));
    H_temp = reshape(H_temp, 1, []);
    H_temp(index_struct.x + max(index_struct.x)*(index_struct.y-1)) = data;
    H_temp = reshape(H_temp, max(index_struct.x), max(index_struct.y));
else
    H_temp = zeros(HS_size);
    H_temp = reshape(H_temp, 1, []);
    H_temp(index_struct.x + HS_size(1)*(index_struct.y-1)) = data;
    H_temp = reshape(H_temp, HS_size(1), HS_size(1));
end
% Actually Plot?
if nargout == 0
    imagesc(H_temp); colormap('gray')
    H_temp = 1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%