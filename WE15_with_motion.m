% This code accompanies the manuscript "How sea ice motion influences sea
% ice extent" by Mason, Wagner, and Eisenman
% submitted for peer-review to Geophys. Rev. Lett. on Aug 25 2020
%
% The model consists of the previously published model in 
% Wagner, Till JW, and Ian Eisenman. "How climate model complexity 
% influences sea ice stability." Journal of Climate 28.10 (2015): 3998-4014.
%
% PLUS
%
% an account of sea ice motion. 
%
% Till Wagner, Aug 2020, wagnert@uncw.edu
%
function [Tfin, Efin] = WE15_with_motion
%%Model parameters (WE15, Table 1 and Section 2d) -------------------------
D  = 0.6;     %diffusivity for heat transport (W m^-2 K^-1)
A  = 193;     %OLR when T = T_m (W m^-2)
B  = 2.1;     %OLR temperature dependence (W m^-2 K^-1)
cw = 9.8;     %ocean mixed layer heat capacity (W yr m^-2 K^-1)
S0 = 420;     %insolation at equator  (W m^-2)
S1 = 338;     %insolation seasonal dependence (W m^-2)
S2 = 240;     %insolation spatial dependence (W m^-2)
a0 = 0.7;     %ice-free co-albedo at equator
a2 = 0.1;     %ice=free co-albedo spatial dependence
ai = 0.4;     %co-albedo where there is sea ice
Fb = 4;       %heat flux from ocean below (W m^-2)
k  = 2;       %sea ice thermal conductivity (W m^2 K^-1)
Lf = 9.5;     %sea ice latent heat of fusion (W yr m^-3)
% Tm = 0;     %melting temp., not included in the eqns below
F  = 0;       %radiative forcing (W m^-2)
cg = 0.01*cw; %ghost layer heat capacity(W yr m^-2 K^-1)
tau = 1e-5;   %ghost layer coupling timescale (yr)
%%The default run in WE15, Fig 2 uses the time-stepping parameters: -------
% n=400; % # of evenly spaced latitudinal gridboxes (equator to pole)
% nt=1e3; % # of timesteps per year (approx lower limit of stability)
% dur=200; % # of years for the whole run
%%For a quicker computation, use the parameters: --------------------------
n  = 100;
nt = 1200;
dur= 30;
dt = 1/nt;
%%Spatial Grid ------------------------------------------------------------
dx = 1/n;     %grid box width
x = (dx/2:dx:1-dx/2)';  %native grid
%%Diffusion Operator (WE15, Appendix A) -----------------------------------
xb = (dx:dx:1.0-dx)';
lambda=D/dx^2*(1-xb.^2); L1=[0; -lambda]; L2=[-lambda; 0]; L3=-L1-L2;
diffop = - diag(L3) - diag(L2(1:n-1),1) - diag(L1(2:n),-1);
%%Definitions for implicit scheme on Tg
cg_tau = cg/tau;
dt_tau = dt/tau;
dc = dt_tau*cg_tau;
kappa = (1+dt_tau)*eye(n)-dt*diffop/cg;
%%Seasonal forcing (WE15 eq.3)
ty = dt/2:dt:1-dt/2;
S=repmat(S0-S2*x.^2,[1,nt])-repmat(S1*cos(2*pi*ty),[n,1]).*repmat(x,[1,nt]);
%%Further definitions
M = B+cg_tau;
aw= a0-a2*x.^2;   %open water albedo
kLf = k*Lf;
% Motion definitions:
speed = 0.20; %ice drift speed in m/s
v_direction = -1;  %1 poleward, -1 equatorward
h_cap = 4;
% Motion cap
E_cap = -1*h_cap*Lf;%*ones(n);
sqx = sqrt((1-x.^2));
%%Initial conditions ------------------------------------------------------
T = 7.5+20*(1-2*x.^2);
Tg = T; E = cw*T;
%%Set up output arrays, saving 100 timesteps/year
vec100 = 1:nt/100:nt*dur;
E100 = zeros(n,dur*100); T100 = E100;
p = 0; m = -1;
save_period = round(nt/100);
%%Integration (see WE15_NumericIntegration.pdf)----------------------------
% Loop over Years ---------------------------------------------------------
month = 9;
for years = 1:dur
    % Loop within One Year-------------------------------------------------
    for i = 1:nt
        %store 100 timesteps per year
        m = m+1;
        if mod(m,save_period)==0
            p = p+1;
            E100(:,p) = E;
            T100(:,p) = T;
        end
        % forcing
        alpha = aw.*(E>0) + ai*(E<0);    % WE15, eq.4
        C =alpha.*S(:,i)+cg_tau*Tg-A+F;
        % surface temperature
        T0 =  C./(M-kLf./E);                 %WE15, eq.A3
        T = E/cw.*(E>=0)+T0.*(E<0).*(T0<0);  %WE15, eq.9
        %%%%----------------------------------------------------------------------
        %setting up monthly velocities
        if i/100 <= month && i/100 > month-1
            v = speed;
        else
            v = 0;
        end
        CFL = v*dt/dx;
        % add Motion
        Etemp = E; 
        Etemp(E>0)=0;
        
        Ebool = Etemp;
        Ebool(Etemp>E_cap) = 1;
        Ebool(Etemp<=E_cap) = 0;
        E_adv2 = Ebool.*Etemp.*v.*v_direction.*dt.*x./sqx;
        Etemp = Etemp.*circshift(Ebool,-1*v_direction);
        E_adv0 = CFL.*Etemp.*sqx;
        E_adv1 = CFL.*circshift(Etemp, v_direction).*sqx;

        if v_direction == 1 %toward pole
            E_adv0(end) = 0; % in the last cell i don't lose anything
            E_adv1(1) = E_adv0(1); % at equator i get what i lose (through symmetry)
        else %equatorward
            E_adv1(end) = 0;
        end
        E_adv = E_adv2 + E_adv1 - E_adv0;
        
        %%%%%----------------------------------------------------------------
        % Forward Euler on E
        E = E+dt*(C-M*T+Fb)+E_adv;%+Eadv2;                 %WE15, eq.A2
        % Implicit Euler on Tg
        Tg = (kappa-diag(dc./(M-kLf./E).*(T0<0).*(E<0)))\(Tg + (dt_tau*(E/cw.*(E>=0)+(ai*S(:,i)-A+F)./(M-kLf./E).*(T0<0).*(E<0))));        %WE15, eq.A1 
    end
    yrs = sprintf('year %d complete',years); disp(yrs)
end
% -------------------------------------------------------------------------
%%output only converged, final year
tfin = linspace(0,1,100);
Efin = E100(:,end-99:end);
Tfin = T100(:,end-99:end);
%
% -------------------------------------------------------------------------
%
%WE15, Figure 2: Default Steady State Climatology -------------------------

% -------------------------------------------------------------------------

set(0,'defaulttextinterpreter','latex')

winter=26;    %time of coldest <T>

summer=76;    %time of warmest <T>

%%compute seasonal ice edge

xi = zeros(1,100);

for j = 1:length(tfin)
    
    Ej = Efin(:,j);
    
    if isempty(find(Ej<0,1))==0
        
        xi(j) = x(find(Ej<0,1,'first'));
        
    else
        
        xi(j) = max(x);
        
    end
    
end

fig = figure(6); clf

% set(fig, 'Position', [2561 361 1920 984]);

%%plot the enthalpy (Fig 2a)

subplot(1,4,1)

clevs = [-40:20:0 50:50:300];

contourf(tfin,x,Efin,clevs)

%%plot ice edge on E

hold on

plot(tfin,xi,'k')

colorbar

%%alternatively, use

% cbarf(Efin,clevs);  %makes nice discrete colorbar

%%(http://www.mathworks.com/matlabcentral/fileexchange/14290-cbarf)

xlabel('$t$ (final year)')

ylabel('$x$')

title('$E$ (J/m$^2$)')

% plot the temperature (Fig 2b)

clevs = -30:5:30;

subplot(1,4,2)

contourf(tfin,x,Tfin,clevs)

% plot ice edge on T

hold on

plot(tfin,xi,'k')

%%plot T=0 contour (the region between ice edge and T=0 contour is the

%%region of summer ice surface melt)

contour(tfin,x,Tfin,[0,0],'r')

colorbar

% cbarf(Tfin,clevs);

xlabel('$t$ (final year)')

ylabel('$x$')

title('$T$ ($^\circ$C)')

%%plot the ice thickness (Fig 2c)

hfin = -(Efin/Lf.*(Efin<0));

subplot(1,4,3)

clevs = 0.0001:.2:4;

contourf(tfin,x,hfin,clevs)

caxis([0 3.5])

%%plot ice edge on h

hold on

contour(tfin,x,hfin,[0,0])

plot([tfin(winter) tfin(winter)], [0 max(x)],'k')

plot([tfin(summer) tfin(summer)], [0 max(x)],'k--')

xlabel('$t$ (final year)')

ylabel('$x$')

title('$h$ (m)')

colorbar

% cbarf(hfin,round(clevs,1));

%%plot temperature profiles (Fig 2d)

subplot(4,4,4)

plot(x,Tfin(:,summer),'k--')

hold on

plot(x,Tfin(:,winter),'k')

plot([0 1], [0 0],'k')

xlabel('$x$')

ylabel('$T$ ($^\circ$C)')

legend('summer','winter','location','southwest')

%%plot ice thickness profiles (Fig 2e)

subplot(4,4,8)

plot(x,hfin(:,summer),'k--')

hold on

plot(x,hfin(:,winter),'k')

plot([0 1], [0 0],'k')

xlim([0.7 1])

xlabel('$x$')

ylabel('$h$ (m)')

%%plot seasonal thickness cycle at pole (Fig 2f)

subplot(4,4,12)

plot(tfin,hfin(end,:),'k')

xlabel('$t$ (final year)')

ylabel('$h_p$ (m)')

% ylim([2 1.1*max(hfin(end,:))])

%%plot ice edge seasonal cycle (Fig 2g)

subplot(4,4,16)

xideg = rad2deg(asin(xi));

plot(tfin,xideg,'k-')

ylim([0 90])

xlabel('$t$ (final year)')

ylabel('$\theta_i$ (deg)')



