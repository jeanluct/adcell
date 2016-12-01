addpath extern

% Number of gridpoints
N = 201;
% Size of one cell
l = 2*pi;
% Physical size of domain
ks = 10; L = ks*l;
% Diffusivity
Diff = 1;
% Integration time
tmax = 20;

% Spatial domain grid.
x = L*(0:N-1)/N; y = x';
[xx,yy] = meshgrid(x,y);

% Streamfunction
U = 1;
psi = (sqrt(2)*U*l/2/pi) * sin(2*pi*xx/l).*sin(2*pi*yy/l);

fprintf('Peclet number = %g\n',U*l/Diff)

Ak = adcell_setup(psi,Diff,L); % fill advection-diffusion sparse matrix
lu = adcell_decomp(Ak);        % LU-decomposition of integrator

% Initial condition.
% Gaussian centred on the middle cell.  It is wide enough to almost
% fill the cell, but not so wide that it spills over the next cell.
l0 = .12*l/2; icx = L/2; icy = icx;
theta = 1/(2*pi*l0^2)*pk(exp(-((xx - icx).^2 + (yy - icy).^2)/(2*l0^2)));

% Integrate the advection-diffusion equation up to tmax.
covar = adcell_integrate(lu,theta,tmax,L);

% Plot covariance.
t = lu.dt*(1:size(covar,1));
figure(2)
plot(t,covar(:,1,1),'.','MarkerSize',10), hold on
hold on
Deff = Diff + U^2*l^2/16/pi^2/Diff;
fprintf('D=%g  Deff=%g\n',Diff,Deff)
plot(t,2*(t-t(end))*Diff + covar(end,1,1),'g--','LineWidth',2)
plot(t,2*(t-t(end))*Deff + covar(end,1,1),'r','LineWidth',2)
hold off
xlabel('$t$','Interpreter','LaTeX','FontSize',22)
ylabel('$\langle x^2 \rangle$','Interpreter','LaTeX','FontSize',22)
set(gca,'FontSize',18,'FontName','Times')
print -dpdf cellular_D1_var.pdf

figure(1)
%print -dpdf -zbuffer -painters cellular_D1.pdf
print -dpng cellular_D1.png
