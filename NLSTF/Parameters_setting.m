function   par   =  Parameters_setting( sf, kernel_type, sz,s0 )
par.h          =    sz(1);
par.w          =    sz(2);
if strcmp(kernel_type, 'Uniform_blur')
    psf        =    ones(sf)/(sf^2);
elseif strcmp(kernel_type, 'Gaussian_blur')
    psf        =    fspecial('gaussian',5,2);
end
  par.fft_B      =    psf2otf(psf,sz);
 



% middlel = round((sz(1)+1)/2);
% middlec = round((sz(2)+1)/2);
% B = zeros(sz(1),sz(2));
% % Starck-Murtagh filter
% B(middlel-2:middlel+2, middlec-2:middlec+2) = [1 4 6 4 1; 4 16 24 16 4; 6 24 36 24 6; 4 16 24 16 4; 1 4 6 4 1];
% % Circularly center B
% B = ifftshift(B);
% % Normalize
% B = B/sum(sum(B));
% % Fourier transform of the filters
% FB = fft2(B);
%    par.fft_B      =    FB;

par.fft_BT     =    conj(par.fft_B);
par.H          =    @(z)H_z(z, par.fft_B, sf, sz,s0 );
par.HT         =    @(y)HT_y(y, par.fft_BT, sf, sz,s0);


  
    