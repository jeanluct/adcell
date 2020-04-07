function [Ak,uxk,uyk] = adfft(psi,Diff,L)
%ADFFT   Fourier transform of the advection-diffusion operator.
%   AK = ADFFT(PSI) returns the discrete Fourier transform AK of the 2D
%   advection operator A = -U d/dx - V d/dy with velocity streamfunction
%   PSI.  The streamfunction is defined on a spatially-periodic grid of size
%   N x N.  The components U and V are given by
%
%     U = d(PSI)/dy,  V = -d(PSI)/dx.
%
%   The FFT of A is returned as a complex matrix AK of size (N^2-1) x
%   (N^2-1), with the constant (k=0) mode dropped.
%
%   AK = ADFFT(PSI,D) returns the discrete Fourier transform AK of the 2D
%   advection-diffusion operator A = (-U d/dx - V d/dy + D laplacian) with
%   diffusivity D.
%
%   AK = ADFFT(PSI,D,L) also specifies the domain size L (default 2*PI).
%
%   [AK,UK,VK] = ADFFT(...) also returns the discrete Fourier transform
%   UK and VK of the velocity components U and V, with the constant (k=0)
%   mode dropped.  UK and VK are column-vectors of size N^2-1.
%
%   See also ADCELL.VEL, ADCELL.DECOMP, ADCELL.INTEGRATE.

%
% Copyright (c) 2016-2020 Jean-Luc Thiffeault <jeanluc@mailaps.org>
%
% See the file LICENSE for copying permission.
%

if nargin < 1 || isempty(psi)
  % Sample streamfunction psi.
  % Number of gridpoints and number of full cells.
  N = 21; ks = 2;
  % Spatial domain grid.
  x = 2*pi*ks*(0:N-1)/N; y = x';
  [xx,yy] = meshgrid(x,y);
  % Streamfunction, with cell-averaged energy 1.
  psi = sqrt(2) * sin(xx).*sin(yy);
end

if nargin < 2 || isempty(Diff), Diff = 0; end
if nargin < 3 || isempty(L), L = 2*pi; end

N = size(psi,1);

if N ~= size(psi,2)
  % TODO: allow nonsquare.  Need to change unpk.
  error('Streamfunction must be a square matrix.')
end

% Fourier domain grid.
k1 = 2*pi/L;
% If N is even, there is one more negative mode than positive, in
% addition to the zero mode. floor takes care of this.
kmin = floor(-(N-1)/2); kmax = floor((N-1)/2); ik = [0 1:kmax kmin:-1];
k = k1*ik;
[kx,ky] = meshgrid(k,k');

% Velocity components, scaled to domain size L = 2*pi/k1.
[ux,uy] = adcell.vel(psi);
ux = k1*ux; uy = k1*uy;

% This is the function that does all the work: return the FFT of
% the u.grad opetator.
Ak = fft2udotgrad(ux,uy,L);

% Verify non-diffusive part is antihermitian.
aherr = full(max(max(abs(Ak+Ak'))));
if aherr > 1e-10
  warning('Antihermitian error = %g',aherr)
end

if Diff ~= 0
  % Diffusive part of operator.
  Ak = Ak + spconvert([(1:N^2)' (1:N^2)' -Diff*pk(kx.^2+ky.^2)]);
end

% Drop the first row and column (constant mode).
Ak = sparsify(Ak(2:end,2:end));

if nargout > 1
  uxk = pk(fft2(ux)); uxk(1) = [];
  if nargout > 2
    uyk = pk(fft2(uy)); uyk(1) = [];
  end
end
