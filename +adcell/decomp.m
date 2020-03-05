function intdecomp = decomp(Ak,dt,intmethod)
%DECOMP   Matrix decomposition of integrator for advection-diffusion equation.
%
%   See also ADCELL.ADFFT, ADCELL.INTEGRATE.

if nargin < 3 || isempty(intmethod)
  % Default integration method.
  intmethod = 'imidpoint';
end

switch lower(intmethod)
 case 'imidpoint'
  % Use implicit midpoint method.
  if nargin < 2 || isempty(dt), dt = .1; end  % default timestep
  B = speye(size(Ak)) - .5*dt*Ak;
  C = speye(size(Ak)) + .5*dt*Ak;
  intdecomp.C = C;
 case 'ieuler'
  % Use implicit Euler method.
  if nargin < 2 || isempty(dt), dt = .025; end  % default timestep
  B = speye(size(Ak)) - dt*Ak;
 otherwise
  error('Unknown integration scheme.')
end

intdecomp.method = intmethod;
intdecomp.dt = dt;

% LU-decomposition of Fourier-space advection-diffusion operator.
[L,U,P,Q,R] = lu(B);  % L and U are the two largest memory hogs.

% Check accuracy of decomposition.
decomperr = full(max(max(abs(B - R*(P'*(L*(U*(Q'))))))));
if decomperr > 1e-10
  warning('Matrix decomposition error = %g',decomperr)
end

intdecomp.L = L;
intdecomp.U = U;
intdecomp.P = P;
intdecomp.Q = Q;
intdecomp.R = R;
