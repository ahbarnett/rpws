% test if spectral density of alpha=0 3D RPW is flat inside the h-ball.
% barnett 8/21/14

clear; verb = 0;
M = 300; ppw = 16; alpha = 0;
[u x] = rpw3dnufft(M, ppw, alpha);

if verb % show u:
  figure; for z=1:M, imagesc(x,x,squeeze(u(:,:,z))); caxis([-2 2]); axis equal; title(sprintf('z=%.3g',x(z))); drawnow; end
end
  
% check spectral density:
Fu = fftshift(fftn(u)); clear u;                 % takes 1 sec (4 cores)
m = max(abs(Fu(:)));
xi = ppw/M*[-M/2:M/2-1];   % corresp wavenumber grid xi. should have unit ball:
ii = find(abs(xi)<=1.1);  % indices to include, slightly larger than unit ball

if verb % animate a slice:
  figure; for z=ii, imagesc(xi(ii),xi(ii),squeeze(abs(Fu(ii,ii,z)))); caxis([0 m]); axis equal; title(sprintf('z=%.3g',x(z))); drawnow; pause(0.1); end
end
  
% dump spectral dens into radial bins:
redges = ppw/M*(0:M/2); redges = redges(redges<1.5); % keep central only
nb = numel(redges)-1;  % # bins
[xxi yyi zzi] = ndgrid(xi,xi,xi);
rxi = sqrt(xxi.^2+yyi.^2+zzi.^2);
bdens = zeros(1,nb);
FFu = abs(Fu).^2;
for i=1:nb, bdens(i) = mean(FFu(rxi>redges(i) & rxi<redges(i+1))); end % few s
rcens = (redges(1:end-1)+redges(2:end))/2;
figure; plot(rcens, bdens, '+-'); xlabel('k radius'); ylabel('mean spec dens');

% (alpha=1 should give uniform-density shell.. to do)
