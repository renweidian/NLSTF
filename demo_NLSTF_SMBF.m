
% NLSTF and NLSTF_SMBF for Hyperspectral image and multispectral image fusion, Version 2.0
% Copyright(c) 2019 Renwei Dian
% All Rights Reserved.
%
% ----------------------------------------------------------------------
% Permission to use, copy, or modify this software and its documentation
% for educational and research purposes only and without fee is here
% granted, provided that this copyright notice and the original authors'
% names appear on all copies and supporting documentation. This program
% shall not be used, rewritten, or adapted as the basis of a commercial
% software or hardware product without first obtaining permission of the
% authors. The authors make no representations about the suitability of
% this software for any purpose. It is provided "as is" without express
% or implied warranty.
%----------------------------------------------------------------------
%
% This is an implementation of the algorithm for Hyperspectral image super-
% resolution from a pair of low-resolution hyperspectral image and a high-
% resolution multispectral image.
% 
% if you use this code, Please cite the following paper:
%
%  R. Dian, L. Fang,and S. Li,  Hyperspectral Image Super-Resolution via Non-local Sparse Tensor Factorization(NLSTF), CVPR, 2017
%  R. Dian, S. Li, L. Fang. T. Lu, and J. Bioucas-Dias,  Non-local Sparse Tensor Factorization for Semi-Blind Hyperspectral
%  and Multispectral Image Fusion (NLSTF_SMBF), IEEE TCYB, 2020

clear
clc
addpath(genpath('tensor_toolbox_2.6'))
addpath(genpath('NLSTF'))
pathstr = fileparts('.\data');
dirname = fullfile(pathstr, 'data','*.mat');
imglist = dir(dirname);
% F = [0.005 0.007 0.012 0.015 0.023 0.025 0.030 0.026 0.024 0.019 0.010 0.004 0     0      0    0     0     0     0     0     0     0     0     0     0     0     0     0    0     0       0  
%     0.000 0.000 0.000 0.000 0.000 0.001 0.002 0.003 0.005 0.007 0.012 0.013 0.015 0.016 0.017 0.02 0.013 0.011 0.009 0.005  0.001  0.001  0.001 0.001 0.001 0.001 0.001 0.001 0.002 0.002 0.003
%     0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.001 0.003 0.010 0.012  0.013  0.022  0.020 0.020 0.018 0.017 0.016 0.016 0.014 0.014 0.013];
F=create_F();
downsampling_scale = 32;
% psf        =    ones(downsampling_scale)/(downsampling_scale^2);

M=512;
N=512;
% par.fft_B      =    psf2otf(psf,[M N]);
% par.fft_BT     =    conj(par.fft_B);
s0=downsampling_scale/2;
% par.H          =    @(z)H_z(z, par.fft_B, downsampling_scale, [M N],s0 );
% par.HT         =    @(y)HT_y(y, par.fft_BT, downsampling_scale,  [M N],s0);




%% NLSTF
for yy=1:length(imglist)
    im_structure =load(fullfile(pathstr, 'data', imglist(yy).name));
    S = im_structure.b;           
    [M,N,L] = size(S);
    S_bar = hyperConvert2D(S);
    par             =    Parameters_setting( downsampling_scale, 'Gaussian_blur', [M N],s0 );
    Y_h_bar            =    par.H(S_bar);
    [Y_h] =hyperConvert3D(Y_h_bar,M/downsampling_scale, N/downsampling_scale);
    Y = hyperConvert3D((F*S_bar), M, N);
    para.W=10; para.H=10;  para.S=14; para.lambda=1e-6; para.K=151; para.lambda1=1e-5; para.lambda2=1e-5;para.lambda3=1e-5;
    t0=clock;
    Z = NLSTF_SMBF(Y_h,Y,F,downsampling_scale,para);
    t5(yy)=etime(clock,t0)
    [psnr5(yy),rmse5(yy), ergas5(yy), sam5(yy), uiqi5(yy),ssim5(yy),DD5(yy),CC5(yy)] = quality_assessment(double(im2uint8(S)), double(im2uint8(Z)), 0, 1.0/downsampling_scale);
end


