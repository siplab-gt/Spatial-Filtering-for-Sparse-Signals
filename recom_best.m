function [c_out, p_outs] = recom_best(c_in)

% [c_out, p_outs] = recom_best(c_in)
% 
% recom_best searches for the permutation of rows and columns of the matrix
% c_in that yields the maximum value for 
%
%            sum_i ( | X_{i,i} | )
%
% The structure p_outs contains the the value for the optimal sum, as well 
% as the value for the missed prediction
%
%            sum_{i!=j} ( | X_{i,j} | )
%
% Written by Adam Charles
% 2/1/2013

%% Initialization

size_ix = size(c_in, 1);

if size_ix ~= size(c_in, 2)
    error('Need to work on a square matrix!')
end

ix_perms = perms(1:size_ix);

c_out = c_in;
score_out = sum(abs(diag(c_in)));

%% Search for max

for kk = 1:size(ix_perms, 1)
    for ll = 1:size(ix_perms, 1)
        c_tmp = c_in(ix_perms(kk, :), ix_perms(ll, :));
        score_tmp = sum(abs(diag(c_tmp)));
        if  score_tmp > score_out
            c_out = c_tmp;
            score_out = score_tmp;
        end
    end
end

p_outs.p_same = score_out;
p_outs.p_diff = sum(sum(c_out)) - score_out;

end
