function [u x] = rpw2dnufft(M, ppw, alpha, eps, task, nopts)
% RPW2DNUFFT - evaluate real 2D random plane wave for alpha=1 or 0, via NUFFT
%
% [u x] = rpw2dnufft(M, ppw, alpha)
% where alpha=1 (monochromatic) or alpha=0 (Fubini-Study ensemble),
% returns a 2d M-by-M grid of real values u, sampled on the grid x in
% each coordinate, from the unit-wavenumber random plane wave ensemble.
% ppw controls the number of points per wavelength; the box size is thus
% M/ppw wavelengths on a side. The output field u has unit variance.
%
% The FINUFFT library is used, which is multi-core, calls FFTW, and is usually
% several times faster than the single-core Greengard-Lee CMCL library.  As
% with any routine based on the FFT, if possible the box size M should be
% chosen to have small prime factors (eg a power of 2) for efficiency.
%
% [u x] = rpw2dnufft(M, ppw, alpha, eps) also controls overall accuracy
% (default eps = 1e-6, also used if eps = [])
%
% [u x] = rpw2dnufft(M, ppw, alpha, eps, task) allows testing:
%  task = 'b' generates error plot for Fourier-Bessel (local expansion)
%             with high order (alpha=1 case), or 0-order aperture func (alpha=1)
%             Note: slow since has to evaluate Bessel; keep M smaller here.
%  task = '1' generates one horiz-traveling plane wave
%  task = 'p' (default) generates RPW as above.
%
% [u x] = rpw2dnufft(M, ppw, alpha, eps, task, nopts) passes nopts to nufft
%  library. Eg: nopts.fftw = 0 (ESTIMATE), 1 (MEASURE) controls FFTW's plan
%  (which is saved within a MATLAB session).
%  nopts.debug=1 shows FINUFFT internal info.
%  nopts.nthreads controls number of threads used by FINUFFT.
%  nopts.cmcl = 0 (use FINUFFT, default), 1 (use CMCL, assuming in path).
%
% Examples: (with timings on 2017 i7-7700HQ laptop)
%
% [u x] = rpw2dnufft(4096, 16, 1);       % uses 1.5 GB
% figure; imagesc(x,x,u); axis xy equal tight;
%
% FFTW timing is strange: this takes 2.5 s, using default FFTW ESTIMATE,
% which is no faster than the single-core CMCL NUFFT library. However,
%
% [u x] = rpw2dnufft(4100, 16, 1);       % takes only 1.5 s.
%
% Better is to use MEASURE to choose the fastest FFTW plan:
% nopts.fftw=1; [u x] = rpw2dnufft(4096, 16, 1, [], [], nopts);
% takes 35 sec when first run, due to FFTW_MEASURE, then repeats are only 0.5 s,
% around 4x faster than CMCL, even with 1 thread. Thus, if you intend to
% generate many samples, first use (the expensive) FFTW_MEASURE.
%
% [u x] = rpw2dnufft(4096, 16, 0);       % takes 0.8 s (after FFTW_MEASURE)
%
% 1000-wavelengths size box @ 16 points per wavelength: (uses ~17 GB)
% [u x] = rpw2dnufft(16348, 16, 1);
% uses FFTW_ESTIMATE and takes 85 sec, but
% nopts.fftw=1; [u x] = rpw2dnufft(16348, 16, 1,[],[],nopts);
% takes 220 s for FFTW_MEASURE, then repeat calls are only 19 s.
%
% At least on my machine, *avoiding* powers of 2 is better (weirdly) if you
% use the default FFTW_ESTIMATE.
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
% 8/30/17, added alpha=1 'b' case.    (c) Alex Barnett 2006-2017.
% Supercedes rpw2dsample.m from 2009, using much less RAM, and including
% case alpha=0. Supercedes rpw2dnufft.m of 8/16/14 which used CMCL NUFFT.
% For work w/ collaborators Sarnak, Canzani, Jin, Konrad, Decourcy-Ireland.
% Includes routines by Andras Pataki.

if nargin<4 | isempty(eps), eps=1e-6; end  % nufft tolerance; gives 8-9 digits
if nargin<5 | isempty(task), task='p'; end
if nargin<6, nopts=[]; end
if alpha~=1 & alpha~=0, error('alpha must be 0 or 1'); end % monochromatic

fprintf('box is %.3g wavelengths on a side\n',M/ppw)
fprintf('predicted RAM usage is %.3g GB\n',9e-8*M^2)
h = 2*pi/ppw;    % coord grid spacing = radius in k-space = Kyle's alpha

if task=='1'  % one horiz-traveling plane wave
  N = 1;
  xj = h; yj = 0; % wavevector (h,0)
  fj = 1;
elseif task=='b'    % high-order J_n bessel test (uex slow for large n)
  if alpha==1
    n = ceil(0.9*h*M/sqrt(2));   % order n to test
    fprintf('testing bessel order n=%d\n',n);
    % (need to check out to p = h*M/sqrt(2), eg 580 for ppw=16, M=2048)
    N = ceil(2.0*h*M);  % for high-n bessel, needs to go up to around predicted
    tj = 2*pi*(1:N)/N;  % must be full circle
    xj = h*cos(tj); yj = h*sin(tj);
    fj = (1/N)*(-1i)^n*exp(1i*n*tj);   % normalized
  else                   % only case n=0 for now
    Nt = ceil(2.0*h*M);  % N per ring (full circle)
    hxi = 2*pi/(2.0*M*h);   % grid spacing in \xi from 0 to 1 in radius
    [r w] = QuadNodesInterval(0,1,0,hxi,1,1,16); % smooth alpert-corr rule
    Nr = numel(r); N = Nt*Nr;   % # radial pts, total # pts
    tj = 2*pi*(1:Nt)/Nt;          % equispaced angles on full circle
    xj = zeros(1,N); yj = xj; fj = xj;
    for i=1:Nr   % loop radially over quadr nodes, with weight 2r
      v = 2*r(i)*w(i);   % quadr wei is applied to variance
      jj = (i-1)*Nt + (1:Nt); % indices
      xj(jj) = h*r(i)*cos(tj); yj(jj) = h*r(i)*sin(tj);
      fj(jj) = 0.5*v/Nt;      % prefactor by trial & error
    end

  end
elseif task=='p'  % Default task: random plane wave sample
  if alpha==1
    N = ceil(1.0*h*M);
    tj = pi*(1:N)/N;   % half circle (ok for Re RPW ensemble)
    xj = h*cos(tj); yj = h*sin(tj);
    fj = (1/sqrt(N))*(randn(N,1)+1i*randn(N,1));  % each plane wave variance = 2
  else
    Nt = ceil(1.0*h*M);  % N per ring (half circle)
    hxi = 2*pi/(2.0*M*h);   % grid spacing in \xi from 0 to 1 in radius
    [r w] = QuadNodesInterval(0,1,0,hxi,1,1,16); % smooth alpert-corr rule
    Nr = numel(r); N = Nt*Nr;   % # radial pts, total # pts
    tj = pi*(1:Nt)/Nt;          % equispaced angles on half-circle
    xj = zeros(1,N); yj = xj; fj = xj;
    for i=1:Nr   % loop radially over quadr nodes, with weight 2r
      v = 2*r(i)*w(i);   % quadr wei is applied to variance
      jj = (i-1)*Nt + (1:Nt); % indices
      xj(jj) = h*r(i)*cos(tj); yj(jj) = h*r(i)*sin(tj);
      fj(jj) = sqrt(v/Nt)*(randn(Nt,1)+1i*randn(Nt,1)); % for unit variance
    end
  end
end
iflag = 1; % choose +i in exp sum
tic(t);
if ~isfield(nopts,'cmcl')
  %addpath('/home/alex/numerics/finufft/matlab');
  u = finufft2d1(xj,yj,fj,iflag,eps,M,M,nopts);
else
  %addpath('/home/alex/numerics/nufftall-1.33/');
  u = nufft2d1(N,xj,yj,fj,iflag,eps,M,M); u = u*N; % old interface
end
fprintf('nufft time = %g s\n',toc(t))
u = real(u)'; % makes real and unit-variance
if task=='p'
  fprintf('mean sq u = %.6g (should be close to 1)\n', mean(u(:).^2));
end

x = h*(-M/2:M/2-1); % coord space grid (must match NUFFT's definition of k grid)

if task=='b' % make error plot
  if M>4096, warning('task=b may be slow; consider shrinking M'); end
  [xx yy] = meshgrid(x); rr = sqrt(xx.^2+yy.^2);
  if alpha==1
    uex = besselj(n,rr).*cos(n*atan2(yy,xx)); % eval Re F-B func
  else
    uex = besselj(1,rr)./rr;                 % rad-symm aperture func
  end
  err = u-uex; fprintf('task=b, J_n test: max err = %.3g\n',max(err(:)))
  figure; imagesc(x,x,log10(abs(err)));
  title('log_{10} error vs J_n(r) e^{in\theta}');
  caxis([-16 0]); colorbar; axis xy equal tight;
end

%================== helper functions: quadrature ===========================

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
% Andras Pataki 2008, edited slightly by Barnett.
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
