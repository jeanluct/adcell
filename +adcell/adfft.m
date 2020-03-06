function [Ak,uxk,uyk] = adfft(psi,Diff,L)
%ADFFT   Fourier transform of the advection-diffusion operator.
%   AK = ADFFT(PSI) returns the discrete Fourier transform AK of the 2D
%   advection operator A = -U d/dx - V d/dy with velocity streamfunction
%   PSI.  The streamfunction is defined on a spatial, periodic grid of size
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
%   See also ADCELL.DECOMP, ADCELL.INTEGRATE.

if nargin < 1 || isempty(psi)
  % Sample streamfunction psi.
  % Number of gridpoints and number of full cells.
  N = 21; ks = 2;
  % Spatial domain grid.
  x = 2*pi*ks*(0:N-1)/N; y = x';
  [xx,yy] = meshgrid(x,y);
  psi = sqrt(2) * sin(xx).*sin(yy);
end

if nargin < 2 || isempty(Diff), Diff = 0; end
if nargin < 3 || isempty(L), L = 2*pi; end

N = size(psi,1);

if N ~= size(psi,2)
  % TODO: allow nonsquare.
  error('Streamfunction must be a square matrix.')
end

% Fourier domain grid.
k1 = 2*pi/L;
% If N is even, there is one more negative mode than positive, in
% addition to the zero mode. floor takes care of this.
kmin = floor(-(N-1)/2); kmax = floor((N-1)/2); ik = [0 1:kmax kmin:-1];
k = k1*ik;
[kx,ky] = meshgrid(k,k');

% Fourier differentiation matrices.
[~,D] = fourdif(N,1); D = k1*D;
ux = D*psi;  % ux =  dpsi/dy (act on columns)
uy = psi*D;  % uy = -dpsi/dx (act on rows)

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
