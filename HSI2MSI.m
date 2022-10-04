function [MSI_data, Xform_mat, MSI_wmat] = HSI2MSI(varargin)

% [MSI_data, Xform_mat, MSI_wmat] = HSI2MSI(varargin)
% 
% Function to calculate equivalent Multispectral Image (MSI) data from
% Hyperspectral Image (HSI) data. The transform may either be input, or the
% function will default to binning based on the WorldViewII Sensor.
% 
% Inputs: 
%       HSI_data    - Hyperspectral data t transform (NxK matrix)
%       w_vec       - wavelenghs for HSI data bins (Nx1 vector)
%       band_bins   - bin banding scheme (Mx2 matrix)
% 
% Outputs: 
%       MSI_data    - Multispectral data (MxK matrix)
%       Xform_mat   - transform for MSI to MSI (MxN matrix)
%       MSI_wmat    - wavelength centers for MSI data (Mx1 vector)
% 
% Uses:
%   [MSI_data, Xform_mat, MSI_wmat] = HSI2MSI(HSI_data, w_vec) - turn
%       HSI_data into MSI data based on the WorldViewII binning and the
%       association to w_vec, the wavelength vector corresponding to the
%       HSI data.
%   [MSI_data, Xform_mat, MSI_wmat] = HSI2MSI(HSI_data, w_vec, band_bins) -
%       turn HSI_data into MSI data based on the binning given by band_bins
%       and its association to w_vec, the wavelength vector corresponding 
%       to the HSI data.
% 
% Last Modified 07/14/2010 - Adam Charles

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Input Checks

if nargin == 3
    HSI_data =  varargin{1};
    w_vec = varargin{2};
    band_bins = varargin{3};
    if numel(band_bins) == 1
        switch band_bins
            case 1
                band_bins = [0.4, 0.45;
                             0.45, 0.51; 
                             0.51, 0.58; 
                             0.585, 0.625;
                             0.63, 0.69; 
                             0.705, 0.745; 
                             0.77, 0.895; 
                             0.860, 1.04
                             1.04, 1.4;
                             1.46, 1.618;
                             1.618, 1.77
                             1.99, 2.24;
                             2.24, 2.49];
            case 2
                band_bins = [0.4, 0.45;
                             0.45, 0.51; 
                             0.51, 0.58; 
                             0.585, 0.625;
                             0.63, 0.69; 
                             0.705, 0.745; 
                             0.77, 0.895; 
                             0.860, 1.04
                             1.04, 1.4;
                             1.46, 1.77
                             1.99, 2.49];
            case 3
                band_bins = [0.435, 0.5472;
                             0.5472, 0.6587; 
                             0.6587, 0.7702; 
                             0.7702, 0.8817;
                             0.8817, 0.9932; 
                             0.9932, 1.1047; 
                             1.1047, 1.2162; 
                             1.2162, 1.3278;
                             1.46, 1.618;
                             1.618, 1.77
                             1.99, 2.24;
                             2.24, 2.49];
            case 4
                band_bins = [0.435, 0.6587; 
                             0.6587, 0.8817;
                             0.8817, 1.1047; 
                             1.1047, 1.3278;
                             1.46, 1.618;
                             1.618, 1.77
                             1.99, 2.24;
                             2.24, 2.49];
            case 5
                band_bins = [0.435, 0.5472;
                             0.5472, 0.6587; 
                             0.6587, 0.7702; 
                             0.7702, 0.8817;
                             0.8817, 0.9932; 
                             0.9932, 1.1047; 
                             1.1047, 1.2162; 
                             1.2162, 1.3278];
        end
    else
        
    end
elseif nargin == 2
    HSI_data =  varargin{1};
    w_vec = varargin{2};
    band_bins = [0.4, 0.45; 0.45, 0.51; 0.51, 0.58; 0.585, 0.625;...
        0.63, 0.69; 0.705, 0.745; 0.77, 0.895; 0.860, 1.04];
else
    error('Bad number of inputs!')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Calculate MSI Data and params

Xform_mat = zeros(size(band_bins, 1), size(HSI_data, 1));
for jj = 1:size(band_bins, 1)
    Xform_mat(jj, :) = (w_vec > band_bins(jj, 1))&(w_vec <= band_bins(jj, 2));
end
MSI_data = Xform_mat*HSI_data;
MSI_wmat = mean(band_bins, 2);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
