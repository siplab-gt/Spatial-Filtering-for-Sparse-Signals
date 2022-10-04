%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% HSI_RWinfer_test
% 
% 01/27/2012 - Adam Charles
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load Data
load('Results/2010_09_20_HSIdata_20010822southend')
load('fixed_known_data.mat')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load Dictionaries
clear basis_cell

comp_bases_path = '';
comp_bases_names = {'2010_07_13_22hyper_01.mat', ...
    '2010_07_17_22hyper_01.mat', '2010_07_14_44hyper_01.mat',...
    '2010_07_19_44hyper_01b.mat'};
c_opts.bnames = {'N = 22 Learned 1', 'N = 22 Learned 2', ...
    'N = 44 Learned 1', 'N = 44 Learned 1'};

for index = 1:length(comp_bases_names)
    load([comp_bases_path, comp_bases_names{index}]);
    basis_cell{index} = basic_cell.dictionary;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Partition a cube to work on

x_choices = 518:585;
y_choices = 271:558;

HS1 = HS_image(hyper_index, data_obj(1, :));
figure
Hpart1 = HS1(x_choices, y_choices);
imagesc(Hpart1);
colormap gray
axis off
axis image

Cube_Seg = zeros(size(Hpart1, 1), size(Hpart1, 2), size(data_obj, 1));
Cube_Seg(:, :, 1) = Hpart1;
parfor kk = 2:size(data_obj, 1)
    temp = HS_image(hyper_index, data_obj(kk, :));
    Cube_Seg(:, :, kk) = temp(x_choices, y_choices);
end

clear HS1 temp

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Compress Cube (To MSI data)

% Make MSI testing data
[MSI_testing_data4, Xform_mat4, MSI_wmat4] = HSI2MSI(testing_data, ...
    wvl_X(wvl_X(:, end) == 1, 2), 4);
Cube_SegM_fs = zeros(size(Hpart1, 1), size(Hpart1, 2), size(MSI_testing_data4, 1));

[MSI_testing_data5, Xform_mat5, MSI_wmat5] = HSI2MSI(testing_data, ...
    wvl_X(wvl_X(:, end) == 1, 2), 5);
Cube_SegM_ms = zeros(size(Hpart1, 1), size(Hpart1, 2), size(MSI_testing_data5, 1));

for kk = 1:size(Hpart1, 1)
    parfor ll = 1:size(Hpart1, 2)
        Cube_SegM_fs(kk, ll, :) = Xform_mat4*squeeze(Cube_Seg(kk, ll, :));
        Cube_SegM_ms(kk, ll, :) = Xform_mat5*squeeze(Cube_Seg(kk, ll, :));
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Recover Cube Without Re-Weighting

opts.lambda = 0.001;
opts.tol = 1e-3;
im_mask = 0;
coef_cube_ms = rwl1sf_infer(Cube_SegM_ms, Xform_mat5*basis_cell{4}, im_mask, opts);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Recover Cube With Reweighting

opts.lambda = 0.001;
im_mask = 1;
coef_cube_ms1 = rwl1sf_infer(Cube_SegM_ms, Xform_mat5*basis_cell{4}, im_mask, opts);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Recover Cube With Correlated Reweighting

opts.lambda = 0.001;
% No mask
IMW = 5; 
IMV = linspace(-2, 2, IMW);
IMvar = 0.8;
im_mask = exp((-(IMV.'*ones(1, IMW)).^2 - (ones(IMW, 1)*IMV).^2)/IMvar);
coef_cube_ms2 = rwl1sf_infer(Cube_SegM_ms, Xform_mat5*basis_cell{4}, im_mask, opts);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Translate results

dictionary4 = basis_cell{4};

clear HRec_ms4 HRec_rw_ms4 Cube_Rec_ms4
for kk = 1:size(Hpart1, 1)
    parfor ll = 1:size(Hpart1, 2)
        HRec_ms4(kk, ll, :) = dictionary4*squeeze(coef_cube_ms(kk, ll, :));
        HRec_rw_ms4(kk, ll, :) = dictionary4*squeeze(coef_cube_ms1(kk, ll, :));
        Cube_Rec_ms4(kk, ll, :) = dictionary4*squeeze(coef_cube_ms2(kk, ll, :));
    end
end

err_mat_ms = sum((HRec_ms4 - Cube_Seg).^2, 3)./sum(Cube_Seg.^2, 3);
err_mat_ms1 = sum((HRec_rw_ms4 - Cube_Seg).^2, 3)./sum(Cube_Seg.^2, 3);
err_mat_ms2 = sum((Cube_Rec_ms4 - Cube_Seg).^2, 3)./sum(Cube_Seg.^2, 3);

%% 
err_mat_ms(err_mat_ms == Inf) = 0;
err_mat_ms1(err_mat_ms1 == Inf) = 0;
err_mat_ms2(err_mat_ms2 == Inf) = 0;
mean_err0 = mean(mean(err_mat_ms));
mean_err1 = mean(mean(err_mat_ms1));
mean_err2 = mean(mean(err_mat_ms2));
median_err0 = median(reshape(err_mat_ms, [], 1));
median_err1 = median(reshape(err_mat_ms1, [], 1));
median_err2 = median(reshape(err_mat_ms2, [], 1));

c_range = [min([min(min(err_mat_ms)), min(min(err_mat_ms1)), ...
    min(min(err_mat_ms2))]), max([max(max(err_mat_ms)),max(max(err_mat_ms1)), ...
    max(max(err_mat_ms2))])];

figure(1)
subplot(3, 1, 1), imagesc(err_mat_ms, c_range)
axis image
axis off
title('BPDN Recovery', 'FontSize', 22)
set(gca, 'FontSize', 18)
h1 = colorbar;
set(h1, 'FontSize', 18)
subplot(3, 1, 2), imagesc(err_mat_ms1, c_range)
% colormap gray
axis image
axis off
title('RWL1 Recovery', 'FontSize', 22)
set(gca, 'FontSize', 18)
h2 = colorbar;
set(h2, 'FontSize', 18)
subplot(3, 1, 3), imagesc(err_mat_ms2, c_range)
axis image
axis off
title('RWL1-SF Recovery', 'FontSize', 22)
set(gca, 'FontSize', 18)
h3 = colorbar;
set(h3, 'FontSize', 18)

figure(2);

j_val = [22, 59]; %58, 59, 30, 27
k_val = [60, 78]; %281, 225, 156, 277

for kk = 1:2
    clear tmp1 tmp2 tmp3 tmp4
    tmp1 = zeros(size(wvl_X, 1), 1);
    tmp1(wvl_X(:, end) == 1) = squeeze(Cube_Seg(j_val(kk), k_val(kk), :));
    tmp2 = zeros(size(wvl_X, 1), 1);
    tmp2(wvl_X(:, end) == 1) = squeeze(HRec_ms4(j_val(kk), k_val(kk), :));
    tmp3 = zeros(size(wvl_X, 1), 1);
    tmp3(wvl_X(:, end) == 1) = squeeze(HRec_rw_ms4(j_val(kk), k_val(kk), :));
    tmp4 = zeros(size(wvl_X, 1), 1);
    tmp4(wvl_X(:, end) == 1) = squeeze(Cube_Rec_ms4(j_val(kk), k_val(kk), :));
    
    max_yval = max([max(tmp1), max(tmp2), max(tmp3), max(tmp4)]);
    subplot(2, 1, kk), plot(wvl_X(:, 2), tmp1(:), '-b', 'LineWidth', 3)
    subplot(2, 1, kk), hold on,
    plot(wvl_X(:, 2), tmp2(:), ':k', 'LineWidth', 3)
    plot(wvl_X(:, 2), tmp3(:), '-.r', 'LineWidth', 3)
    plot(wvl_X(:, 2), tmp4(:), '--c', 'LineWidth', 3)
    if kk == 2
        legend('Actual Spectrum', 'BPDN', 'RWL1', 'RWL1-SF')
    end
    xlabel('Wavelength ({\mu}m)', 'FontSize', 25)
    ylabel('Reflectance', 'FontSize', 25)
    set(gca, 'FontSize', 20, 'Xlim', [min(wvl_X(:, 2)), max(wvl_X(:, 2))], ...
        'Ylim', [0, 1.02*max_yval])
    hold off
end

%%

wave_range = wvl_X(:, 2);

figure, hold on;
plot(wave_range(Xform_mat5(1, :)==1), ones(sum(Xform_mat5(1, :)==1)), 'r', 'Linewidth', 35)
color_cell = {'y', 'm', 'b', 'c', 'k', 'g', 'b'};
for kk = 2:8
    temp = wave_range(Xform_mat5(kk-1, :)==1);
    w_start = temp(end);
    temp = wave_range(Xform_mat5(kk, :)==1);
    w_end = temp(end);
    plot(wave_range((wave_range >= w_start)&(wave_range<=w_end)), kk*ones(sum(Xform_mat5(kk, :)==1)+1), ...
        color_cell{kk-1}, 'Linewidth', 35)
end
temp = wave_range(Xform_mat5(8, :)==1);
w_end = temp(end);
plot(wave_range(61:68), 10*ones(numel(61:68)), 'k', 'Linewidth', 35)
plot(wave_range(91:95), 10*ones(numel(91:95)), 'k', 'Linewidth', 35)

plot(wave_range(68:91), 9*ones(numel(68:91)), 'b', 'Linewidth', 35)
plot(wave_range(95:end), 9*ones(numel(wave_range(95:end))), 'b', 'Linewidth', 35)
plot(wave_range(59:61), 9*ones(numel(59:61)), 'b', 'Linewidth', 35)
set(gca, 'Ylim', [0,11], 'Xlim', [wave_range(1), wave_range(end)], 'YTick', 1:10, ...
    'YTickLabel', {'Band 1', 'Band 2', 'Band 3', 'Band 4', 'Band 5', 'Band 6', ...
    'Band 7', 'Band 8', 'HSI Only', 'Water Ab.'}, 'FontSize', 20)

%%

c_range = [min([min(min(err_mat_ms)), min(min(err_mat_ms1)), ...
    min(min(err_mat_ms2))]), max([max(max(err_mat_ms)),max(max(err_mat_ms1)), ...
    max(max(err_mat_ms2))])];

figure(2)
subplot(2, 2, 1), imagesc(sum((Cube_Seg).^2, 3))
axis image
axis off
title('Histogram - Direct CS', 'FontSize', 22)
set(gca, 'FontSize', 18)
colorbar
subplot(2, 2, 2), imagesc(err_mat_ms, c_range)
axis image
axis off
title('Histogram - Direct CS', 'FontSize', 22)
set(gca, 'FontSize', 18)
colorbar
subplot(2, 2, 3), imagesc(err_mat_ms1, c_range)
axis image
axis off
title('Histogram - RW-CS', 'FontSize', 22)
set(gca, 'FontSize', 18)
colorbar
subplot(2, 2, 4), imagesc(err_mat_ms2, c_range)
axis image
axis off
title('Histogram - RW-Corr CS', 'FontSize', 22)
set(gca, 'FontSize', 18)
colorbar

%% 
mean_err0 = mean(mean(err_mat_ms));
mean_err1 = mean(mean(err_mat_ms1));
mean_err2 = mean(mean(err_mat_ms2));
median_err0 = median(reshape(err_mat_ms, [], 1));
median_err1 = median(reshape(err_mat_ms1, [], 1));
median_err2 = median(reshape(err_mat_ms2, [], 1));

max_err = max([max(reshape(err_mat_ms, [], 1)), max(reshape(err_mat_ms1, [], 1)), max(reshape(err_mat_ms2, [], 1))]);

c_range = [min([min(min(err_mat_ms)), min(min(err_mat_ms1)), ...
    min(min(err_mat_ms2))]), max([max(max(err_mat_ms)),max(max(err_mat_ms1)), ...
    max(max(err_mat_ms2))])];



figure(1)
subplot(2, 1, 1), imagesc(err_mat_ms, c_range)
axis image
axis off
title('BPDN Recovery Errors', 'FontSize', 22)
set(gca, 'FontSize', 18)
h1 = colorbar;
set(get(h1, 'ylabel'), 'String', 'rMSE', 'FontSize', 18)
set(h1, 'FontSize', 18)


RGB_image = zeros(size(Cube_Seg, 1), size(Cube_Seg, 2), 3);
RGB_image(:, :, 1) = sum(Cube_Seg(:, :, 12:18), 3);
RGB_image(:, :, 2) = sum(Cube_Seg(:, :, 6:11), 3);
RGB_image(:, :, 3) = sum(Cube_Seg(:, :, 1:5), 3);
RGB_image = RGB_image/max(max(max(RGB_image)));

subplot(2, 1, 2), imagesc(RGB_image)
axis image
axis off
title('RGB Image', 'FontSize', 22)
set(gca, 'FontSize', 18)

[~, x_bins] = hist(reshape(err_mat_ms, [], 1), 100);


CDF_vals = [cumsum(hist(reshape(err_mat_ms, [], 1), x_bins)).', ...
    cumsum(hist(reshape(err_mat_ms1, [], 1), x_bins)).', ...
    cumsum(hist(reshape(err_mat_ms2, [], 1), x_bins)).']/numel(err_mat_ms);

figure, hold on; % subplot(2, 2, [2,4]),
plot([0, x_bins], [0; CDF_vals(:, 1)], ':k', 'LineWidth', 5)
plot([0, x_bins], [0; CDF_vals(:, 2)], '-.r', 'LineWidth', 5)
plot([0, x_bins], [0; CDF_vals(:, 3)], '--c', 'LineWidth', 5)
% title('CDF of rMSE', 'FontSize', 22)
legend('BPDN', 'RWL1', 'RWL1-SF')
xlabel('rMSE', 'FontSize', 20)
axis square
set(gca, 'FontSize', 18, 'XLim', [0, max_err], 'YLim', [0, 1.03])

%%

x_range_plot4 = (22:32) + 6;
y_range_plot4 = (59:69) + 9;

point_vec_p4 = -2:0.1:2;
Xplot4 = ones(numel(point_vec_p4), 1)*point_vec_p4;
Yplot4 = point_vec_p4'*ones(1, numel(point_vec_p4));
Gfunc =  0.1*exp(-(Yplot4.^2 + Xplot4.^2)/2);

G_func_superX = -5:0.1:5;
G_func_superY = -5:0.1:5;
G_func_super = zeros(numel(G_func_superX), numel(G_func_superY));
G_func_super(15+(1:numel(point_vec_p4)), 55+(1:numel(point_vec_p4))) = Gfunc;

el = 22;
az = -42;

figure(1)
subplot(1, 2, 1), surf(-5:5, -5:5, 0*ones(size(err_mat_ms(x_range_plot4 , y_range_plot4))), ...
    err_mat_ms(x_range_plot4 , y_range_plot4),'EdgeColor','black'); 
axis off
colormap gray
axis square
set(gca, 'XLim', [-5,5], 'YLim', [-5,5], 'ZLim', 0.25*[-1,1])
view(az,el)


subplot(1, 2, 2), hold on;
surf(G_func_superX, G_func_superY, G_func_super, ...
    max(max(err_mat_ms(x_range_plot4 , y_range_plot4)))*G_func_super/max(max(G_func_super)),...
    'EdgeColor', 'none', 'FaceColor', 0.7*[1, 1, 1]); 
surf(G_func_superX(1:10:end), G_func_superY, G_func_super(:, 1:10:end), ...
    max(max(err_mat_ms(x_range_plot4 , y_range_plot4)))*G_func_super(:, 1:10:end)/max(max(G_func_super)),...
    'EdgeColor', 'black', 'FaceColor', 'none', 'MeshStyle', 'column', ...
    'LineWidth', 2); 
surf(G_func_superX, G_func_superY(1:10:end), G_func_super(1:10:end, :), ...
    max(max(err_mat_ms(x_range_plot4 , y_range_plot4)))*G_func_super(1:10:end, :)/max(max(G_func_super)),...
    'EdgeColor', 'black', 'FaceColor',  'none', 'MeshStyle', 'row', ...
    'LineWidth', 2);
axis off
colormap gray
axis square
set(gca, 'XLim', [-5,5], 'YLim', [-5,5], 'ZLim', 0.25*[-1,1])
view(az,el)
hold off;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
