function Deff = effdiff

addpath extern

% Parameters.
N = 201; L = 2*pi; Diff = .001;

x = L*(0:N-1)/N; y = x'; [xx,yy] = meshgrid(x,y);

% Streamfunction, with total energy U^2.
U = 1; psi = (sqrt(2)*U*L/2/pi) * sin(2*pi*xx/L).*sin(2*pi*yy/L);

% FFT of advection-diffusion operator and velocity components.
[Ak,uxk,uyk] = adcell.adfft(psi,Diff,L);

chixk = Ak \ uxk;
chiyk = Ak \ uyk;

chix = ifft2(unpk([0 ; chixk]),'symmetric');
chiy = ifft2(unpk([0 ; chiyk]),'symmetric');

figure(1)
set(pcolor(x,y,chix),'EdgeColor','none')
axis square, axis xy
xlabel('x'), ylabel('y')

figure(2)
set(pcolor(x,y,chiy),'EdgeColor','none')
axis square, axis xy
xlabel('x'), ylabel('y')

intfg = @(f,g) real(f'*g / N^4);

intfg(uxk,uxk) + intfg(uyk,uyk)

Dxx = -intfg(uxk,chixk);
Dxy = -intfg(uxk,chiyk);
Dyx = -intfg(uyk,chixk);
Dyy = -intfg(uyk,chiyk);

% Effective diffusion tensor.
Deff = [Dxx Dxy; Dyx Dyy]

U^2*L^2/16/pi^2 / Diff
