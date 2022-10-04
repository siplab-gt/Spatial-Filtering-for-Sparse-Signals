function C_out = match_vecs(C_in, C_ref, varargin)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Parse Inputs

if nargin < 2
    error('Need at least two inputs!')
end

if nargin > 2
    greed_type = varargin{1};
else
    greed_type = 'greedy';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Run comparison

C_tmp = C_in;
C_out = zeros(size(C_in));
C_tmp_norms = sqrt(sum(C_tmp.^2, 1));
C_ref_norms = sqrt(sum(C_ref.^2, 1));

if strcmp(greed_type, 'greedy') == 1
    for kk = 1:size(C_ref, 2)
        comp_vals = abs((C_tmp')*C_ref(:, kk))./(sum(C_ref_norms(kk).^2)*C_tmp_norms');
        ix = find(comp_vals == max(comp_vals), 1);
        C_out(:, kk) = C_tmp(:, ix);
        C_tmp = [C_tmp(:, 1:ix-1), C_tmp(:, ix+1:end)];
        C_tmp_norms = [C_tmp_norms(1:ix-1), C_tmp_norms(ix+1:end)];
    end
elseif strcmp(greed_type, 'global') == 1
    norms_mat = (C_tmp_norms')*C_ref_norms;
    C_ref_tmp = C_ref;
    ind_tmp = 1:size(C_ref, 2);
    ind_ref_tmp = 1:size(C_ref, 2);
    
    for kk = 1:size(C_ref, 2)
        comp_vals = abs((C_tmp')*C_ref_tmp./norms_mat);
        [ix, jx] = find(comp_vals == max(min(comp_vals)), 1);
        C_out(:, ind_ref_tmp(jx)) = C_in(:, ind_tmp(ix));
        C_tmp = [C_tmp(:, 1:ix-1), C_tmp(:, ix+1:end)];
        ind_tmp = [ind_tmp(:, 1:ix-1), ind_tmp(:, ix+1:end)];
        C_ref_tmp = [C_ref_tmp(:, 1:jx-1), C_ref_tmp(:, jx+1:end)];
        ind_ref_tmp = [ind_ref_tmp(:, 1:jx-1), ind_ref_tmp(:, jx+1:end)];
        norms_mat = norms_mat([1:ix-1, ix+1:end], [1:jx-1, jx+1:end]);
    end
    
else
    error('Unknown search type!')
end



end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%