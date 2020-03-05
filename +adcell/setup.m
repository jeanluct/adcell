function [Ak,ux,uy] = setup(psi,Diff,L)

if nargin < 3, L = 2*pi; end

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
if full(max(max(abs(Ak+Ak')))) > 1e-10
  error('Somehow the matrix is not antihermitian...')
end

% Diffusive part of operator.
Ak = Ak + spconvert([(1:N^2)' (1:N^2)' -Diff*pk(kx.^2+ky.^2)]);

% Drop the first row and column (constant mode).
Ak = sparsify(Ak(2:end,2:end));
