function adcell_run(runname)

if nargin < 1, runname = '1'; end

addpath extern

switch runname
 % classic cellular, large diff
 case {'1','run1','low-Pe'}
  N = 201; l = 2*pi; ks = 10; L = ks*l; Diff = 1; tmax = 20;
  x = L*(0:N-1)/N; y = x'; [xx,yy] = meshgrid(x,y);

  % Streamfunction
  U = 1; psi = (sqrt(2)*U*l/2/pi) * sin(2*pi*xx/l).*sin(2*pi*yy/l);

  Ak = adcell.adfft(psi,Diff,L); % fill advection-diffusion sparse matrix
  lu = adcell.decomp(Ak);        % LU-decomposition of integrator

  % Initial condition.
  l0 = .12*l/2; icx = L/2; icy = icx;
  theta = 1/(2*pi*l0^2)*pk(exp(-((xx - icx).^2 + (yy - icy).^2)/(2*l0^2)));

  % Integrate the advection-diffusion equation up to tmax.
  covar = adcell.integrate(lu,theta,tmax,L);

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
 
 % classic cellular, moderate diff
 case {'2'}
  N = 301; l = 2*pi; ks = 8; L = ks*l; Diff = .1; tmax = 40;
  x = L*(0:N-1)/N; y = x'; [xx,yy] = meshgrid(x,y);

  % Streamfunction
  U = 1; psi = (sqrt(2)*U*l/2/pi) * sin(2*pi*xx/l).*sin(2*pi*yy/l);

  Ak = adcell.adfft(psi,Diff,L); % fill advection-diffusion sparse matrix
  lu = adcell.decomp(Ak);        % LU-decomposition of integrator

  % Initial condition.
  l0 = .12*l/2; icx = L/2; icy = icx;
  theta = 1/(2*pi*l0^2)*pk(exp(-((xx - icx).^2 + (yy - icy).^2)/(2*l0^2)));
  adcell.integrate(lu,theta,tmax,L)

 % classic cellular, small diff
 case {'3'}
  N = 501; l = 2*pi; ks = 6; L = ks*l; Diff = .01; tmax = 40;
  x = L*(0:N-1)/N; y = x'; [xx,yy] = meshgrid(x,y);

  % Streamfunction
  U = 1; psi = (sqrt(2)*U*l/2/pi) * sin(2*pi*xx/l).*sin(2*pi*yy/l);

  Ak = adcell.adfft(psi,Diff,L); % fill advection-diffusion sparse matrix
  lu = adcell.decomp(Ak);        % LU-decomposition of integrator

  % Initial condition.
  l0 = .12*l/2; icx = L/2; icy = icx;
  theta = 1/(2*pi*l0^2)*pk(exp(-((xx - icx).^2 + (yy - icy).^2)/(2*l0^2)));

  % Integrate the advection-diffusion equation up to tmax.
  adcell.integrate(lu,theta,tmax,L)

 % classic cellular, very small diff
 case {'4'}
  N = 1201; l = 2*pi; ks = 6; L = ks*l; Diff = .001; tmax = 40;
  x = L*(0:N-1)/N; y = x'; [xx,yy] = meshgrid(x,y);

  % Streamfunction
  U = 1; psi = (sqrt(2)*U*l/2/pi) * sin(2*pi*xx/l).*sin(2*pi*yy/l);

  Ak = adcell.adfft(psi,Diff,L); % fill advection-diffusion sparse matrix
  lu = adcell.decomp(Ak);        % LU-decomposition of integrator

  % Initial condition.
  l0 = .12*l/2; icx = L/2; icy = icx;
  theta = 1/(2*pi*l0^2)*pk(exp(-((xx - icx).^2 + (yy - icy).^2)/(2*l0^2)));
  adcell.integrate(lu,theta,tmax,L)

 % truncated ABC, closed streamlines
 % See Crisanti et al. (1990).
 case {'5'}
  N = 401; l = 2*pi; ks = 6; L = ks*l; Diff = .01; tmax = 40;
  x = L*(0:N-1)/N; y = x'; [xx,yy] = meshgrid(x,y);

  % Streamfunction
  A = -1; B = 1;
  psi = (l/2/pi) * (B*sin(2*pi*yy/l) + A*cos(2*pi*xx/l));

  Ak = adcell.adfft(psi,Diff,L); % fill advection-diffusion sparse matrix
  lu = adcell.decomp(Ak);        % LU-decomposition of integrator

  % Initial condition.
  l0 = .12*l/2; icx = L/2; icy = icx;
  theta = 1/(2*pi*l0^2)*pk(exp(-((xx - icx).^2 + (yy - icy).^2)/(2*l0^2)));
  adcell.integrate(lu,theta,tmax,L)

 % truncated ABC, open streamlines
 % See Crisanti et al. (1990).
 case {'6'}
  N = 401; l = 2*pi; ks = 6; L = ks*l; Diff = .01; tmax = 40;
  x = L*(0:N-1)/N; y = x'; [xx,yy] = meshgrid(x,y);

  % Streamfunction
  A = -1.3; B = 1;
  psi = (l/2/pi) * (B*sin(2*pi*yy/l) + A*cos(2*pi*xx/l));

  Ak = adcell.adfft(psi,Diff,L); % fill advection-diffusion sparse matrix
  lu = adcell.decomp(Ak);        % LU-decomposition of integrator

  % Initial condition.
  l0 = .12*l/2; icx = L/2; icy = icx;
  theta = 1/(2*pi*l0^2)*pk(exp(-((xx - icx).^2 + (yy - icy).^2)/(2*l0^2)));
  adcell.integrate(lu,theta,tmax,L)
end
