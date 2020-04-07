function [ux,uy] = vel(psi)
%VEL   Velocity components associated with a stremfunction.
%   [U,V] = VEL(PSI) returns the velocity components U and V corresponding
%   to the streamfunction PSI, defined on a spatially-periodic grid of size
%   N x N.  The components U and V are given by
%
%     U = d(PSI)/dy,  V = -d(PSI)/dx.
%
%   The physical domain is assumed to be (2*pi) x (2*pi).
%
%   See also ADCELL.ADFFT, ADCELL.DECOMP, ADCELL.INTEGRATE.

%
% Copyright (c) 2016-2020 Jean-Luc Thiffeault <jeanluc@mailaps.org>
%
% See the file LICENSE for copying permission.
%

N = size(psi,1);

if N ~= size(psi,2)
  % TODO: allow nonsquare.  Need to change unpk.
  error('Streamfunction must be a square matrix.')
end

% Fourier differentiation matrices.
[~,D] = fourdif(N,1);
ux = D*psi;  % ux =  dpsi/dy (act on columns)
uy = psi*D;  % uy = -dpsi/dx (act on rows)
