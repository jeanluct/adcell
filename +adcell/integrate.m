function varargout = integrate(intdecomp,theta0,tmax,L)
%INTEGRATE   Integrate advection-diffusion equation.
%
%   See also ADCELL.ADFFT, ADCELL.DECOMP.

%
% Copyright (c) 2016-2020 Jean-Luc Thiffeault <jeanluc@mailaps.org>
%
% See the file LICENSE for copying permission.
%

% Record movie frames?
if nargout > 1, domovie = true; else, domovie = false; end

if nargin < 4 || isempty(L), L = 2*pi; end
if nargin < 3 || isempty(tmax), tmax = 20; end    % Integration time

dt = intdecomp.dt;
Lo = intdecomp.L;
Up = intdecomp.U;
P = intdecomp.P;
Q = intdecomp.Q;
R = intdecomp.R;

N = sqrt(size(Lo,1)+1);
if N ~= floor(N)
  error('Something''s wrong... matrices should have size N^2-1 x N^2-1.')
end

thetak = pk(fft2(unpk(theta0)));
thetak0 = thetak(1);
thetak = thetak(2:end);

% Spatial domain grid (for computing covariance).
x = L*(0:N-1)/N; y = x';
[xx,yy] = meshgrid(x,y);

close all

set(0,'Units','pixels')
scnsize = get(0,'ScreenSize');
fig1 = figure(1);
set(fig1,'Position',[20 50 .48*scnsize(3) .4*scnsize(3)])
% Add some points to get a smoother plot.
Nplot = N+100;

% The time integration: iterate, solving the appropriate linear
% system for the chosen implicit scheme.
Niter = ceil(tmax/dt);
covar = zeros(Niter,2,2);
thetavar = zeros(Niter,1);
for it = 1:Niter
  switch lower(intdecomp.method)
   case 'imidpoint'
    thetak = Q*(Up\(Lo\(P*(R\(intdecomp.C*thetak)))));
   case 'ieuler'
    thetak = Q*(Up\(Lo\(P*(R\thetak))));
   otherwise
    error('Unknown integration scheme.')
  end

  % Plot the concentration.
  theta = pk(ifft2(unpk([thetak0;thetak]),'symmetric'));
  thetaplot = refine2(unpk(theta),Nplot);
  xplot = L*(0:Nplot-1)/Nplot; yplot = xplot';
  figure(1)
  set(pcolor(xplot,yplot,thetaplot),'EdgeColor','none')
  axis square, axis xy
  xlabel('x'), ylabel('y')
  title(sprintf('t = %g',it*dt))
  drawnow
  if domovie, Mov(it) = getframe; end

  % Compute covariance.
  int = sum(sum(unpk(theta)));
  meanx = sum(sum(xx.*unpk(theta))) / int;
  meany = sum(sum(yy.*unpk(theta))) / int;
  covar(it,1,1) = sum(sum(xx.^2.*unpk(theta)))/int - meanx^2;
  covar(it,1,2) = sum(sum(xx.*yy.*unpk(theta)))/int - meanx*meany;
  covar(it,2,1) = covar(it,1,2);
  covar(it,2,2) = sum(sum(yy.^2.*unpk(theta)))/int - meany^2;

  % Compute scalar variance.
  thetavar(it) = var(theta);
end

if nargout > 0, varargout{1} = covar; end
if nargout > 1, varargout{2} = Mov; end
if nargout > 2, varargout{3} = thetavar; end
