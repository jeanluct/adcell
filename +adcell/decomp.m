function intdecomp = decomp(Ak,dt,intmethod)

if nargin < 3, intmethod = 'imidpoint'; end  % Integration method

switch lower(intmethod)
 case 'imidpoint'
  % Use implicit midpoint method.
  if nargin < 2, dt = .1; end
  B = speye(size(Ak)) - .5*dt*Ak;
  C = speye(size(Ak)) + .5*dt*Ak;
  intdecomp.C = C;
 case 'ieuler'
  % Use implicit Euler method.
  if nargin < 2, dt = .025; end
  B = speye(size(Ak)) - dt*Ak;
 otherwise
  error('Unknown integration scheme.')
end

intdecomp.method = intmethod;
intdecomp.dt = dt;

% LU-decomposition of Fourier-space advection-diffusion operator.
[L,U,P,Q,R] = lu(B);  % L and U are the two largest memory hogs.

% Check accuracy of decomposition.
if full(max(max(abs(B - R*(P'*(L*(U*(Q')))))))) > 1e-10
  error('The LU decomposition seems dubious...')
end

intdecomp.L = L;
intdecomp.U = U;
intdecomp.P = P;
intdecomp.Q = Q;
intdecomp.R = R;
