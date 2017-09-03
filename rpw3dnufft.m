function [u x] = rpw3dnufft(M, ppw, alpha, eps, task, nopts)
% RPW3DNUFFT - evaluate real 3D random plane wave for alpha=1 or 0, via NUFFT
%
% [u x] = rpw3dnufft(M, ppw, alpha)
% where alpha=1 (monochromatic) or alpha=0 (Fubini-Study ensemble),
% returns a 3d M-by-M-by-M grid of real values u, sampled on the grid x on
% each of the three axes, from the unit-wavenumber random plane wave ensemble.
% ppw controls the number of points per wavelength; the cube is thus
% M/ppw wavelengths on a side. The output field u has unit variance.
%
% The FINUFFT library is used, which is multi-core, calls FFTW, and is usually
% several times faster than the single-core Greengard-Lee CMCL library.  As
% with any routine based on the FFT, if possible the box size M should be
% chosen to have small prime factors (eg a power of 2) for efficiency.
%
% [u x] = rpw3dnufft(M, ppw, alpha, eps) also controls overall accuracy
% (default eps = 1e-6, also used if eps = [])
%
% [u x] = rpw3dnufft(M, ppw, alpha, eps, task) allows testing:
%  task = 'b' generates error plot for spherical-harmonic (local expansion)
%             m=n=0 basis function (alpha=1 case), or j_1(r)/r (alpha=0).
%             Note: slow since has to evaluate Bessel; keep M smaller here.
%  task = '1' generates one z-traveling plane wave
%  task = 'p' (default) generates RPW as above.
%
% [u x] = rpw3dnufft(M, ppw, alpha, eps, task, nopts) passes nopts to nufft
%  library. Eg: nopts.fftw = 0 (ESTIMATE), 1 (MEASURE) controls FFTW's plan
%  (which is saved within a MATLAB session).
%  nopts.debug=1 shows FINUFFT internal info.
%  nopts.nthreads controls number of threads used by FINUFFT.
%  nopts.cmcl = 0 (use FINUFFT, default), 1 (use CMCL, assuming in path).
%
% Examples: (with timings on 2017 i7-7700HQ laptop)
%
% 20 wavelengths cube @ 16 ppw:
% M = 320; [u x] = rpw3dnufft(M, 16, 1);       % takes 6 sec, uses around 5 GB
%
% animate slices through it:
% figure; for z=1:M, imagesc(x,x,squeeze(u(:,:,z))); caxis([-2 2]); axis equal; title(sprintf('z=%.3g',x(z))); drawnow; end
%
% 50 wavelengths cube @ 10 ppw, showing FINUFFT info:
% M=500; o.debug=1; [u x] = rpw3dnufft(M, 10, 1,[],[],o);   % 11 s, around 20 GB
% M=500; o.debug=1; [u x] = rpw3dnufft(M, 10, 0,[],[],o);   % 14 s, around 20 GB
%
% Note that here alpha=0, despite needing 1e7 quadr pts, takes only 3 s more.
%
% Strangely, using FFTW_ESTIMATE (as default above) is *slower* for powers of
% two, against the usual wisdom. Thus, unless you use FFTW_MEASURE, avoid
% powers of two. See also examples in rpw2dnufft.m
%
% Dependencies:
%
% https://github.com/ahbarnett/finufft must be installed
% and built, and its matlab directory must be in the matlab path, eg via:
% appdath('/path/to/finufft/matlab');
%
% If nopts.cmcl used, http://www.cims.nyu.edu/cmcl/nufft/nufft.html must be
% installed and built, and its top directory in the matlab path.
%
% Notes:
% 8/30/17, added alpha=1 'b' case.  (c) Alex Barnett 2006-2017.
% Supercedes rpw3dnufft.m of 8/16/14 which used CMCL NUFFT.
% For work w/ collaborators Sarnak, Canzani, Jin, Konrad, Decourcy-Ireland.
% Includes routines by Andras Pataki, Greg von Winckel.

if nargin<4 | isempty(eps), eps=1e-6; end  % nufft tolerance; gives 8-9 digits
if nargin<5 | isempty(task), task='p'; end
if nargin<6, nopts=[]; end
if alpha~=1 & alpha~=0, error('alpha must be 0 or 1'); end % monochromatic

fprintf('box is %.3g wavelengths on a side\n',M/ppw)
fprintf('predicted RAM usage is %.2g GB\n',1.6e-7*M^3)
h = 2*pi/ppw; % coord grid spacing = radius in k-space = Kyle's alpha

if task=='1'  % one horiz-traveling plane wave
  N = 1;
  xj = 0; yj = 0; zj = h; % wavevector (0,0,h)
  fj = 1;
  
elseif task=='b'
  if alpha==1      % spherical j_0(r) bessel func is FT of const func on S^2
    Np = ceil(1.0*h*M);   % # uniform sample pts in 0<=phi<pi
    pj = pi*(1:Np)/Np;    % half circle
    Nz = ceil(1.0*h*M);   % # sample pts in -1<z<1
    [z w] = lgwt(Nz,-1,1);  % vertical quadr scheme
    N = Np*Nz; xj = zeros(1,N); yj = xj; zj = xj; fj = xj; % alloc arrays
    for j=1:Nz
      jj = (j-1)*Np + (1:Np); % indices
      rho = sqrt(1-z(j)^2);   % circle radius in xy plane
      xj(jj) = h*rho*cos(pj); yj(jj) = h*rho*sin(pj); zj(jj) = h*z(j);
      fj(jj) = (0.5/Np)*w(j);   % norm factor 0.5 found empirically
    end
  else      % spherical j_1(r)/r bessel func is FT of char func of unit ball
    Np = ceil(1.0*h*M);   % # uniform sample pts in 0<=phi<pi
    pj = pi*(1:Np)/Np;    % half circle
    Nz = ceil(1.0*h*M);   % # sample pts in -1<z<1
    [z w] = lgwt(Nz,-1,1);  % vertical quadr scheme
    hxi = 2*pi/(2.0*M*h);   % grid spacing in \xi if from 0 to 1 in radius
    [r wr] = QuadNodesInterval(0,1,0,hxi,1,1,16); % smooth Alpert-corr rule
    Nr = numel(r); N = Np*Nz*Nr;   % # radial pts, total # pts
    xj = zeros(1,N); yj = xj; zj = xj; fj = xj; % alloc arrays
    for i=1:Nr  % loop over radial quadr pts
      for j=1:Nz  % loop over vertical quadr pts
        jj = (i-1)*(Nz*Np) + (j-1)*Np + (1:Np); % indices
        rho = sqrt(1-z(j)^2);   % circle radius in xy plane
        xj(jj) = h*r(i)*rho*cos(pj); yj(jj) = h*r(i)*rho*sin(pj);
        zj(jj) = h*r(i)*z(j);
        v = r(i)^2*wr(i)*w(j);   % but why is the 3/2 needed below.. ?
        fj(jj) = (0.5/Np) * v;    % prefactor guessed
      end
    end
  end
    
elseif task=='p' % ...... random plane waves .....................
  if alpha==1            % random data on S^2
    Np = ceil(1.0*h*M);   % # uniform sample pts in 0<=phi<pi
    pj = pi*(1:Np)/Np;    % half circle
    Nz = ceil(1.0*h*M);   % # sample pts in -1<z<1
    [z w] = lgwt(Nz,-1,1);  % vertical quadr scheme
    N = Np*Nz; xj = zeros(1,N); yj = xj; zj = xj; fj = xj; % alloc arrays
    for j=1:Nz
      jj = (j-1)*Np + (1:Np); % indices
      rho = sqrt(1-z(j)^2);   % circle radius in xy plane
      xj(jj) = h*rho*cos(pj); yj(jj) = h*rho*sin(pj); zj(jj) = h*z(j);
      v = 0.5*w(j);               % quadr wei applied to variance
      fj(jj) = sqrt(v/Np)*(randn(Np,1)+1i*randn(Np,1)); % unit variance
    end
  else           % alpha=0: integrate over radii, ie solid sphere radius h
    Np = ceil(1.0*h*M);   % # uniform sample pts in 0<=phi<pi
    pj = pi*(1:Np)/Np;    % half circle
    Nz = ceil(1.0*h*M);   % # sample pts in -1<z<1
    [z w] = lgwt(Nz,-1,1);  % vertical quadr scheme
    hxi = 2*pi/(2.0*M*h);   % grid spacing in \xi if from 0 to 1 in radius
    [r wr] = QuadNodesInterval(0,1,0,hxi,1,1,16); % smooth Alpert-corr rule
    Nr = numel(r); N = Np*Nz*Nr;   % # radial pts, total # pts
    xj = zeros(1,N); yj = xj; zj = xj; fj = xj; % alloc arrays
    for i=1:Nr  % loop over radial quadr pts
      for j=1:Nz  % loop over vertical quadr pts
        jj = (i-1)*(Nz*Np) + (j-1)*Np + (1:Np); % indices
        rho = sqrt(1-z(j)^2);   % circle radius in xy plane
        xj(jj) = h*r(i)*rho*cos(pj); yj(jj) = h*r(i)*rho*sin(pj);
        zj(jj) = h*r(i)*z(j);
        v = 3/2*r(i)^2*wr(i)*w(j);  % NB r^2; quadr wei applied to variance
        fj(jj) = sqrt(v/Np)*(randn(Np,1)+1i*randn(Np,1)); % unit variance
      end
    end
  end
end
iflag = 1; % choose +i in exp sum
t=tic;
if ~isfield(nopts,'cmcl')
  %addpath('/home/alex/numerics/finufft/matlab');
  u = finufft3d1(xj,yj,zj,fj,iflag,eps,M,M,M,nopts);
else
  %addpath('/home/alex/numerics/nufftall-1.33/');
  u = nufft3d1(N,xj,yj,zj,fj,iflag,eps,M,M,M); u = u*N; % old interface
end
fprintf('nufft time = %g s\n',toc(t))
u = real(u);  % makes real and unit-variance
if task=='p'
  fprintf('mean sq u = %.6g (should be close to 1)\n', mean(u(:).^2));
end

x = h*(-M/2:M/2-1); % coord space grid (must match NUFFT's definition of k grid)
u =reshape(u, [M M M]); % output index ordering is (x,y,z)

if task=='b' % make error plot
  [xx yy zz] = meshgrid(x); rr = sqrt(xx.^2+yy.^2+zz.^2);
  if alpha==1
    uex = sin(rr)./rr;  % spherical Bessel j_0(r), roundoff around r->0
  else
    uex = (sin(rr)./rr - cos(rr))./(rr.^2);  % j_1(r)/r, bad roundoff r->0
  end
  err = u-uex; fprintf('task=b, j_0 test: max err = %.3g\n',max(err(:)))
  figure; for z=1:M, imagesc(x,x,squeeze(log10(abs(err(:,:,z)))));
    caxis([-16 0]); axis equal; title(sprintf('z=%.3g',x(z)));
    drawnow; end
%  figure; imagesc(x,x,log10(abs(squeeze(err(1,:,:))))); title('slice error vs j_0(r)'); caxis([-16 0]); colorbar; axis xy equal tight;
%  figure; imagesc(x,x,squeeze(u(1,:,:)));
end


%================== helper functions: quadrature ===========================

function [x w] = lgwt(N,a,b)
% [x w] = lgwt(N,a,b)
%
% This script is for computing definite integrals using Legendre-Gauss 
% Quadrature. Computes the N Legendre-Gauss nodes x and weights w on an
% interval [a,b]
%
% Suppose you have a continuous function f(x) which is defined on [a,b]
% which you can evaluate at any x in [a,b]. Simply evaluate it at all of
% the values contained in the x vector to obtain a vector f. Then compute
% the definite integral using sum(f.*w);
%
% Written by Greg von Winckel - 02/25/2004
N=N-1;
N1=N+1; N2=N+2;

xu=linspace(-1,1,N1)';

% Initial guess
y=cos((2*(0:N)'+1)*pi/(2*N+2))+(0.27/N1)*sin(pi*xu*N/N2);

% Legendre-Gauss Vandermonde Matrix
L=zeros(N1,N2);

% Derivative of LGVM
Lp=zeros(N1,N2);

% Compute the zeros of the N+1 Legendre Polynomial
% using the recursion relation and the Newton-Raphson method

y0=2;

% Iterate until new points are uniformly within epsilon of old points
while max(abs(y-y0))>eps
    
    L(:,1)=1;
    Lp(:,1)=0;
    
    L(:,2)=y;
    Lp(:,2)=1;
    
    for k=2:N1
        L(:,k+1)=( (2*k-1)*y.*L(:,k)-(k-1)*L(:,k-1) )/k;
    end
 
    Lp=(N2)*( L(:,N1)-y.*L(:,N2) )./(1-y.^2);   
    
    y0=y;
    y=y0-L(:,N2)./Lp;
    
end

% Linear map from[-1,1] to [a,b]
x=(a*(1-y)+b*(1+y))/2;      

% Compute the weights
w=(b-a)./((1-y.^2).*Lp.^2)*(N2/N1)^2;

x = x(end:(-1):1);


%----------------------------------------------------------------------
function [Ax, Aw] = QuadNodesInterval(a, b, N, h, corra, corrb, order)
% Return the quadrature nodes and weights on an interval [a,b]
% with N points or step size approximately h (if N is set to 0)
% corra, corrb determine the endpoint corrections used at the two endpoints:
%        0: none
%        1: smooth function
%        2: square root singularity
%        3: log singularity
% order: is the order of the endpoint corrections used
% The nodes (Ax) and the weights (Aw) are returned
%
% Andras Pataki 2008
    if (corra == 0)
	NodesToSkipL = 0;
    elseif (corra == 1)
	[ExtraNodesL, ExtraWeightsL, NodesToSkipL] = QuadSmoothExtraPtNodes(order);
    else
	error(sprintf('Unknown endpoint corretion method corra=%d', corra));
    end

    if (corrb == 0)
	NodesToSkipR = 0;
    elseif (corrb == 1)
	[ExtraNodesR, ExtraWeightsR, NodesToSkipR] = QuadSmoothExtraPtNodes(order);
    else
	error(sprintf('Unknown endpoint corretion method corrb=%d', corrb));
    end

    if (N < 2)
	N = max(2, ceil(abs(b-a)/h));
    end
    
    % We have to have enough nodes to skip
    N = max(NodesToSkipL + NodesToSkipR, N);

    % Generate the regular nodes
    N1 = N-1; 
    h = (b-a)/N1;
    Ax = a + (0:N1)' * h;
    Aw = ones(size(Ax)) * h;
    Aw(1) = 0.5*h;
    Aw(length(Ax)) = 0.5*h;
    
    % Add the left endpoint corrections
    if (NodesToSkipL > 0)
	Ax = [a+ExtraNodesL*h; Ax(NodesToSkipL+1:length(Aw)) ];
	Aw = [ExtraWeightsL*h; Aw(NodesToSkipL+1:length(Aw)) ];
    end

    % Add the right endpoint corrections
    if (NodesToSkipR > 0)
	Ax = [ Ax(1:length(Ax)-NodesToSkipR); flipud(b-ExtraNodesR*h) ];
	Aw = [ Aw(1:length(Aw)-NodesToSkipR); flipud(ExtraWeightsR*h) ];
    end
    
%----------------------------------------------------------------------
function [ExtraNodes, ExtraWeights, NodesToSkip] = QuadSmoothExtraPtNodes(order)
% Andras Pataki 2008

    if (order <= 3)
	M = [
	    1.666666666666667E-01 5.0000000000 00000E-01
	];
	NodesToSkip = 1;

    elseif (order <= 4)
	M = [
	    2.000000000000000E-01 5.208333333333333E-01
	    1.000000000000000E+00 9.791666666666667E-01
	];
	NodesToSkip = 2;

    elseif (order <= 5)
	M = [
	    2.245784979812614E-01 5.540781643606372E-01
	    1.013719374359164E+00 9.459218356393628E-01
	];
	NodesToSkip = 2;

    elseif (order <= 6)
	M = [
	    2.250991042610971E-01 5.549724327164180E-01
	    1.014269060987992E+00 9.451317411845473E-01
	    2.000000000000000E+00 9.998958260990347E-01
	];
	NodesToSkip = 3;

    elseif (order <= 7)
	M = [
	    2.180540672543505E-01 5.408088967208193E-01
	    1.001181873031216E+00 9.516615045823566E-01
	    1.997580526418033E+00 1.007529598696824E+00
	];
	NodesToSkip = 3;

    elseif (order <= 8)
	M = [
	    2.087647422032129E-01 5.207988277246498E-01
	    9.786087373714483E-01 9.535038018555888E-01
	    1.989541386579751E+00 1.024871626402471E+00
	    3.000000000000000E+00 1.000825744017291E+00
	];
	NodesToSkip = 4;

    elseif (order <= 12)
	M = [
	    7.023955461621939E-02 1.922315977843698E-01
	    4.312297857227970E-01 5.348399530514687E-01
	    1.117752734518115E+00 8.170209442488760E-01
	    2.017343724572518E+00 9.592111521445966E-01
	    3.000837842847590E+00 9.967143408044999E-01
	    4.000000000000000E+00 9.999820119661890E-01
	];
	NodesToSkip = 5;

    elseif (order <= 16)
	M = [
	    9.919337841451028E-02 2.528198928766921E-01
	    5.076592669645529E-01 5.550158230159486E-01
	    1.184972925827278E+00 7.852321453615224E-01
	    2.047493467134072E+00 9.245915673876714E-01
	    3.007168911869310E+00 9.839350200445296E-01
	    4.000474996776184E+00 9.984463448413151E-01
	    5.000007879022339E+00 9.999592378464547E-01
	    6.000000000000000E+00 9.999999686258662E-01
	];
	NodesToSkip = 7;

    elseif (order <= 20)
	M = [
	    9.209200446233291E-02 2.351836144643984E-01
	    4.752021947758861E-01 5.248820509085946E-01
	    1.124687945844539E+00 7.634026409869887E-01
	    1.977387385642367E+00 9.284711336658351E-01
	    2.953848957822108E+00 1.010969886587741E+00
	    3.976136786048776E+00 1.024959725311073E+00
	    4.994354281979877E+00 1.010517534639652E+00
	    5.999469539335291E+00 1.001551595797932E+00
	    6.999986704874333E+00 1.000061681794188E+00
	    8.000000000000000E+00 1.000000135843597E+00
	];
	NodesToSkip = 9;

    elseif (order <= 24)
	M = [
	    6.001064731474805E-02 1.538932104518340E-01
	    3.149685016229433E-01 3.551058128559424E-01
	    7.664508240518316E-01 5.449200036280007E-01
	    1.396685781342510E+00 7.104078497715549E-01
	    2.175195903206602E+00 8.398780940253654E-01
	    3.062320575880355E+00 9.272767950890611E-01
	    4.016440988792476E+00 9.750605697371132E-01
	    5.002872064275734E+00 9.942629650823470E-01
	    6.000285453310164E+00 9.992421778421898E-01
	    7.000012964962529E+00 9.999534370786161E-01
	    8.000000175554469E+00 9.999990854912925E-01
	    9.000000000000000E+00 9.999999989466828E-01
	];
	NodesToSkip = 10;

    elseif (order <= 28)
	M = [
	    6.234360533194102E-02 1.595975279734157E-01
	    3.250286721702614E-01 3.637046028193864E-01
	    7.837350794282182E-01 5.498753177297441E-01
	    1.415673112616924E+00 7.087986792086956E-01
	    2.189894250061313E+00 8.335172275501195E-01
	    3.070053877483040E+00 9.204446510608518E-01
	    4.018613756218047E+00 9.710881776552090E-01
	    5.002705902035397E+00 9.933296578555239E-01
	    5.999929741810400E+00 9.994759087910050E-01
	    6.999904720846024E+00 1.000133030254421E+00
	    7.999986894843540E+00 1.000032915011460E+00
	    8.999999373380393E+00 1.000002261653775E+00
	    9.999999992002911E+00 1.000000042393520E+00
	    1.100000000000000E+01 1.000000000042872E+00
	];
	NodesToSkip = 12;

    elseif (order <= 32)
	M = [
	    5.899550614325259E-02 1.511076023874179E-01
	    3.082757062227814E-01 3.459395921169090E-01
	    7.463707253079130E-01 5.273502805146873E-01
	    1.355993726494664E+00 6.878444094543021E-01
	    2.112943217346336E+00 8.210319140034114E-01
	    2.987241496545946E+00 9.218382875515803E-01
	    3.944798920961176E+00 9.873027487553060E-01
	    4.950269202842798E+00 1.018251913441155E+00
	    5.972123043117706E+00 1.021933430349293E+00
	    6.989783558137742E+00 1.012567983413513E+00
	    7.997673019512965E+00 1.004052289554521E+00
	    8.999694932747039E+00 1.000713413344501E+00
	    9.999979225211805E+00 1.000063618302950E+00
	    1.099999938266130E+01 1.000002486385216E+00
	    1.199999999462073E+01 1.000000030404477E+00
	    1.300000000000000E+01 1.000000000020760E+00
	];
	NodesToSkip = 14;

    else
	error(sprintf('QuadSmoothExtraPtNodes: order %d is too high', order));
    end

    ExtraNodes = M(:,1);
    ExtraWeights = M(:,2);
