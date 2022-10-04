function coef_cube = rwl1sf_infer(data_cube, spec_dict, im_mask, opts)

% Decompose a Data cube using the reweighted-l1 spatial filtering scheme
% 
% data_mat   = cube of data indexed by (x, y, l) - first two dimensions are
%              spacial, and the 3rd dimension is spectral
% spec_dict  = spectral dictionary - MxN matrix where M is the number of
%              spectral bands in the HSI image and N is the number of
%              dictionary elements
% im_mask    = (x,y) indexed mask representing correlations between pixels
% opts       = options for optimization program
% 
% 1/26/2012 - Adam Charles


%% Some initial calculations
%
[X, Y, M] = size(data_cube);           % Get image cube sizes
N = size(spec_dict, 2);                % Get number of dictionary elements
num_pix = numel(data_cube(:, :, 1));   % Total number of pixels

% Initialize Matrices
coef_mat = zeros(N, num_pix);          % Matrix of sparse coefficients
weight_mat = ones(size(coef_mat));     % Matrix of coefficent variances (weights)
weight_cube = ones(X, Y, N);           % Cube of sparsity weights

% Set algorithm to use l1ls_nneg for sparse inference
infer_handle = @l1ls_nneg_wrapper;


%% Perform the Decomposition

data_mat = reshape(permute(data_cube, [3, 1, 2]), M, []);
if (numel(im_mask) == 1)&&(im_mask(1) == 0)
    % Perform the sparse decomposition at every pixel
    coef_mat = gen_multi_infer(spec_dict, data_mat, infer_handle, opts);
    % reshape the coefficient matrix to match the HSI data cube
    coef_cube = permute(reshape(coef_mat, [N, X, Y]), [2, 3, 1]);
    fprintf('Finished decomposition.\n')
else
    % Perform initial decomposition 
    coef_mat = gen_rw_multi_infer(spec_dict, data_mat, weight_mat, ...
        infer_handle, opts);
    % reshape the coefficient matrix to match the HSI data cube
    coef_cube = permute(reshape(coef_mat, [N, X, Y]), [2, 3, 1]);
    fprintf('Finished iteration 1.\n')
    for kk = 2:4
        % Update weights. First get the convolution of the coefficient spread with the 
        % image mask kernel
        parfor ll = 1:N
            weight_cube(:, :, ll) = conv2(squeeze(coef_cube(:, :, ll)), im_mask, 'same');
        end
        % Reshape the weight cube to match the HSI size
        weight_mat = reshape(permute(weight_cube, [3, 1, 2]), N, []);
        weight_mat = 0.75./(abs(weight_mat) + 0.2);
        % Decompose every spectrum in the data cube
        coef_mat = gen_rw_multi_infer(spec_dict, data_mat, weight_mat,...
            infer_handle, opts);
        coef_cube = permute(reshape(coef_mat, [N, X, Y]), [2, 3, 1]);
        fprintf('Finished iteration %d.\n', kk)
    end
end

end

