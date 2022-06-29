% Code to generate Figure 9 of Moulton et al. (2023).
%
% The net horizontal displacement of particles in linear internal waves is
% calculated using the equation in Franks et al. (2020).
%
% The net horizontal displacement of particles in solitary internal waves
% is calculated by advecting particles in wave flow fields as described in
% Garwood et al. (2020) and coded in the function IWtracks.m included in
% the folder Functions
%
% Dependencies:
%   - IWtracks.m
%   - ddz.m (required by IWtracks.m from Smyth et al. 2011)
%
% References:
%
% Franks et al. 2020. Stokes drift of plankton in linear internal waves: 
%   Cross-shore transport of neutrally buoyant and depth-keeping organisms. 
%   Limnology and Oceanography
%
% Garwood et al. 2020. Life in Internal Waves. Oceanography.
%
% Moulton et al. 2023. Exchange of plankton, pollutants, and particles 
%   across the nearshore region. Annual Review of Marine Science.
%
% Smyth, W.D., J.N. Moum and J.D. Nash. 2011. Narrowband, high-frequency 
% 	oscillations at the equator. Part II: Properties of shear instabilities.
%

%% Addpath to functions
addpath('../Functions/')

%% Set color for third wave
w3col   = [0 0.6 0.6];

%% Wave parameters - 3 different amplitudes
z       = (-1:0.01:0)';
H       = 1;
amp1    = 0.05;
amp2    = 0.1;
amp3    = 0.15;
phi     = -sin(pi*z);
T       = 20*60;
cp      = 0.1;
z_org   = -1:0.01:0;
deltat  = 1;
nsteps  = 7000;  % 5000 steps was not enough
wshape  = 'sech2';

% Linear stratification
rho_bot = 1025;  % bottom density
rho_top = 1024;
mu      = (rho_bot/rho_top - 1)/(H);

rho0    = rho_bot * (1 - mu * z);

omega   = 2 * pi / T;
k       = omega/cp;

%% Linear wave - Stokes drift passive organisms (Franks et al. 2020)
uSt.p1   = (amp1^2 * pi^2 * omega) ./ (2 * k * H^2) .* cos(2 * pi * z / H);
uSt.p2   = (amp2^2 * pi^2 * omega) ./ (2 * k * H^2) .* cos(2 * pi * z / H);
uSt.p3   = (amp3^2 * pi^2 * omega) ./ (2 * k * H^2) .* cos(2 * pi * z / H);

%% Sech2 profile - total transport
[dk1b, p1b, wave1b] = IWtracks(z, H, phi, rho0, amp1, T, cp, z_org, ...
    deltat, nsteps, wshape);
[dk2b, p2b, wave2b] = IWtracks(z, H, phi, rho0, amp2, T, cp, z_org, ...
    deltat, nsteps, wshape);
[dk3b, p3b, wave3b] = IWtracks(z, H, phi, rho0, amp3, T, cp, z_org, ...
    deltat, nsteps, wshape);

%% Plot net horizontal displacement of particles
fig2 = figure(2);
clf
set(gcf, 'Paperunits', 'inches' )
set(gcf, 'PaperSize', [7.5 3.4])
set(gcf,'PaperPosition',[0 0 7.5 3.4])
set(gcf, 'Units', 'inches')

s1 = subplot(141);
% Integrating Eulerian velocities in linear waves leads to 0 at all depths
plot(NaN, 'k', 'linewidth', 1.5), hold on
plot(NaN, 'm', 'linewidth', 1.5), hold on
plot(NaN, 'color', w3col, 'linewidth', 1.5)

    ylim([-1 0])
    xlim([-0.225 0.225])
    grid on
    
    v0 = xline(0);
    set(v0, 'color', w3col, 'linestyle', '-', 'linewidth', 1.5);
    
    set(gca, 'xtick', [-0.2:0.1:0.2])
    
    xlabel('\Deltax_p/\lambda')
    ylabel('z/h')
    title('Sum of velocities at a point')
    
    % Legend
    str_l1 = sprintf('A/h = 0.05\t\t Fr = 0.16');
    str_l2 = sprintf('A/h = 0.10\t\t Fr = 0.32');
    str_l3 = sprintf('A/h = 0.15\t\t Fr = 0.47');

    lg1 = legend(str_l1, str_l2, str_l3, 'box', 'off', ...
        'location', 'northoutside');
    lg1.ItemTokenSize = [15, 20];
    
s2 = subplot(142);
% Plot the Franks et al. (2020) solution
plot(uSt.p1 * T / (T * cp), z_org, 'k', 'linewidth', 1.5), hold on
plot(uSt.p2 * T / (T * cp), z_org, 'm', 'linewidth', 1.5), hold on
plot(uSt.p3 * T / (T * cp), z_org, 'color', w3col, 'linewidth', 1.5)

    xlim([-0.225 0.225])
    grid on
    
    set(gca, 'linewidth', 1, 'fontsize', 9, 'xtick', [-0.2:0.1:0.2])
    
    xlabel('\Deltax_p/\lambda')
    title('Organisms in wave')
        
s3 = subplot(143);
% Plot integrated Eulerian velocities in solitary wave 
dta = mean(diff(wave1b.xgrid))./cp;

plot(sum(wave1b.u .* dta, 2) ./ wave1b.lambda, z_org, 'k', 'linewidth', 1.5), hold on
plot(sum(wave2b.u .* dta, 2) ./ wave2b.lambda, z_org, 'm', 'linewidth', 1.5)
plot(sum(wave3b.u .* dta, 2) ./ wave3b.lambda, z_org, ...
    'color', w3col, 'linewidth', 1.5)

    xlim([-1 3.5])
    grid on
          
    xlabel('\Deltax_p/\lambda')
    ylabel('z/h')
    title('Sum of velocities at a point')
    
    % Legend
    str_l1 = sprintf('A/h = 0.05\t\t Fr = 0.31');
    str_l2 = sprintf('A/h = 0.10\t\t Fr = 0.63');
    str_l3 = sprintf('A/h = 0.15\t\t Fr = 0.94');
    
    lg2 = legend(str_l1, str_l2, str_l3, 'box', 'off', ...
        'location', 'NorthOutside');
    lg2.ItemTokenSize = [15, 20];
      
s4 = subplot(144);
% Plot advected particles in solitary wave as in Garwood et al. (2020)
plot(p1b.deltax ./ wave1b.lambda, z_org, 'k', 'linewidth', 1.5), hold on
plot(p2b.deltax ./ wave2b.lambda, z_org, 'm', 'linewidth', 1.5)
plot(p3b.deltax ./ wave3b.lambda, z_org, ...
    'color', w3col, 'linewidth', 1.5)

    xlim([-1 3.5])

    grid on
     
    xlabel('\Deltax_p/\lambda')
    title('Organisms in wave')
    
% Axes look
set(s1, 'linewidth', 1, 'fontsize', 9)
set(s2, 'linewidth', 1, 'fontsize', 9, 'yticklabels', [])
set(s3, 'linewidth', 1, 'fontsize', 9)
set(s4, 'linewidth', 1, 'fontsize', 9, 'yticklabels', [])

% Layout
posW = 1.45;
posH = 2;

hsp = 0.2;

posx1 = 0.5;
posx2 = posx1 + posW + hsp;
posx3 = posx2 + posW + 3*hsp;
posx4 = posx3 + posW + hsp;

posy = 0.45;

set([s1 s2 s3 s4 lg1 lg2], 'units', 'inches');
set(s1, 'position', [posx1 posy posW posH]);
set(s2, 'position', [posx2 posy posW posH]);
set(s3, 'position', [posx3 posy posW posH]);
set(s4, 'position', [posx4 posy posW posH]);  

    posL1       = get(lg1, 'position');
    posL1(1)    = posx1 + 0.5*(posW + hsp); 
    set(lg1, 'position', posL1, 'box', 'on', 'linewidth', 1);
    
    posL2       = get(lg2, 'position');
    posL2(1)    = posx3 + 0.5*(posW + hsp); 
    set(lg2, 'position', posL2, 'box', 'on', 'linewidth', 1);
    
% Panel IDs
pID_x   = 0.08;
pID_y   = posH*0.92;

text(s1, pID_x, pID_y, 'a', 'units', 'inches', ...
    'fontsize', 12, 'fontname', 'arial', 'fontweight', 'bold')
text(s2, pID_x, pID_y, 'b', 'units', 'inches', ...
    'fontsize', 12, 'fontname', 'arial', 'fontweight', 'bold')
text(s3, pID_x, pID_y, 'c', 'units', 'inches', ...
    'fontsize', 12, 'fontname', 'arial', 'fontweight', 'bold')
text(s4, pID_x, pID_y, 'd', 'units', 'inches', ...
    'fontsize', 12, 'fontname', 'arial', 'fontweight', 'bold')

%% Print
% exportgraphics('-f2','Fig9_IW.pdf','ContentType','vector','Resolution',1200)
% saveas(fig2, 'Fig9_IW.svg')
% saveas(fig2, 'Fig9_IW.pdf')
