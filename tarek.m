function tarek

% Parameters.
N = 501; l = 2*pi; ks = 1; L = ks*l; Diff = .01; dt = .05; tmax = 10;

x = L*(0:N-1)/N; y = x'; [xx,yy] = meshgrid(x,y);

% Streamfunction
U = 1; psi = (sqrt(2)*U*l/2/pi) * sin(2*pi*xx/l).*sin(2*pi*yy/l);

Ak = adcell_setup(psi,Diff,L); % fill advection-diffusion sparse matrix
lu = adcell_decomp(Ak,dt);     % LU-decomposition of integrator

% Initial condition.
l0 = .1*l/2; icx = L/2; icy = icx;
theta = 1/(2*pi*l0^2)*pk(exp(-((xx - icx).^2 + (yy - icy).^2)/(2*l0^2)));

% Integrate the advection-diffusion equation up to tmax.
[covar,~,thetavar] = adcell_integrate(lu,theta,tmax,L);

figure(2)
t = lu.dt*(1:size(covar,1));
semilogy(t,thetavar)
