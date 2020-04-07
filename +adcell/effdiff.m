function [varargout] = effdiff(psi,Diff,L)
%EFFDIFF   Effective diffusivity for a cellular flow.
%   DEFF = EFFDIFF(PSI,D) returns the effective diffusion tensor DEFF for
%   the advection-diffusion equation with streamfunction PSI and molecular
%   diffusivity D (see <strong>adcell.adfft</strong>).
%
%   See also ADCELL.VEL, ADCELL.ADFFT.

if nargin < 3 || isempty(L)
  L = 2*pi;
end

if nargin < 1 || isempty(psi)
  % Sample streamfunction psi.
  N = 31;
  % Spatial domain grid.
  x = 2*pi*(0:N-1)/N; y = x';
  [xx,yy] = meshgrid(x,y);
  % Streamfunction, with cell-averaged energy 1.
  psi = sqrt(2) * sin(xx).*sin(yy);
end

if nargin < 2 || isempty(Diff)
  Diff = .1;
end

% FFT of advection-diffusion operator and velocity components.
[Ak,uxk,uyk] = adcell.adfft(psi,Diff,L);

% Solve the cell problem for chi.
chixk = Ak \ uxk;
chiyk = Ak \ uyk;

chix = ifft2(unpk([0 ; chixk]),'symmetric');
chiy = ifft2(unpk([0 ; chiyk]),'symmetric');

if ~nargout
  figure(1)
  set(pcolor(x,y,chix),'EdgeColor','none')
  axis square, axis xy
  xlabel('x'), ylabel('y')
  title('chi_x')

  figure(2)
  set(pcolor(x,y,chiy),'EdgeColor','none')
  axis square, axis xy
  xlabel('x'), ylabel('y')
  title('chi_y')
end

intfg = @(f,g) real(f'*g / N^4);  % inner product

% Integrate to find the enhanced diffusivity -<u chi>.
Dxx = -intfg(uxk,chixk);
Dxy = -intfg(uxk,chiyk);
Dyx = -intfg(uyk,chixk);
Dyy = -intfg(uyk,chiyk);

% Effective diffusion tensor.
Deff = diag([Diff Diff]) + [Dxx Dxy; Dyx Dyy];

% Enhancement over molecular value.
%fprintf('%g\n',U^2*L^2/16/pi^2 / Diff)

if nargout > 0
  varargout{1} = Deff;

  if nargout > 1
    varargout{2} = chix;
    if nargout > 2
      varargout{3} = chiy;
    end
  end
end
