% test if spectral density of alpha=0 2D RPW is flat inside the h-ball.
% barnett 8/30/17

clear; verb = 0;
M = 4000; ppw = 10; alpha = 0;
[u x] = rpw2dnufft(M, ppw, alpha);

if verb, figure; imagesc(x,x,u); caxis([-2 2]); axis equal; title('u'); end
  
% check spectral density:
Fu = fftshift(fft2(u)); clear u;
m = max(abs(Fu(:)));
xi = ppw/M*[-M/2:M/2-1];   % corresp wavenumber grid xi. should have unit ball:
ii = find(abs(xi)<=1.2);  % indices to include, slightly larger than unit ball

figure; imagesc(xi(ii),xi(ii),abs(Fu(ii,ii))); caxis([0 m]); axis equal;
title('F');

% dump spectral dens into radial bins:
redges = ppw/M*(0:M/2); redges = redges(redges<1.5); % keep central only
nb = numel(redges)-1;  % # bins
[xxi yyi] = meshgrid(xi,xi);
rxi = sqrt(xxi.^2+yyi.^2);
bdens = zeros(1,nb);
FFu = abs(Fu).^2;
for i=1:nb, bdens(i) = mean(FFu(rxi>redges(i) & rxi<redges(i+1))); end % few s
rcens = (redges(1:end-1)+redges(2:end))/2;
figure; plot(rcens, bdens, '+-'); xlabel('k radius'); ylabel('mean spec dens');

% (alpha=1 should give uniform-density shell.. to do)
