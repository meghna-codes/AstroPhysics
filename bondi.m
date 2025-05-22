%%
format long
E = 0.01;
gamma = 4.0/3.0;
n=1;
syms rr dudrc
% FINDING CRITICAL POINTS
% rc1 = solve(rr==(4*(E-phipw(rr))*(gamma-1))/(dphipw(rr)*(gamma+1)));
rc = 0.75/E;
% FINDING VELOCITY GRADIENTS AT CRITICAL POINTS
% acs = rc*dphipw(rc)/2; % square of acoustic velocity at critical point
acs = 1/(2*rc);
ac = sqrt(acs); % acoustic velocity at critical point
% dudr_rc1 = solve(dudrc == (ddphipw(rc) + (1/n)*sqrt(2*dphipw(rc)/rc)*dudrc + (2*dphipw(rc)/(rc*n)) + dphipw(rc))/((1/n)*dudrc + (1/n)*sqrt(2*dphipw(rc)/rc) - dudrc));
% dudr_rc2 = vpa(dudr_rc1);
% dudr_rc = double(dudr_rc2);
dudr1 = (4*(ac^3)/(2*n + 1))*(1+sqrt(n*(n-1.5))); % flow velocity at critical point
dudr2 = (4*(ac^3)/(2*n + 1))*(1-sqrt(n*(n-1.5)));
% DEFINING INITIAL CONDITIONS
delr1 = 1.0e-1;
i1 = ac + (dudr1*delr1);
i2= ac + (dudr2*delr1);
i3 = 1.5*ac;
i4 = 0.5*ac;
% init1= ac + (dudr_rc(1)*delr);
% init2= ac + (dudr_rc(2)*delr);
% DEFINING RANGES OF INTEGRATION
r1 = rc+delr1:0.1:1000.0;
r2 = rc-delr1:-0.1:1.0;
r3 = 127.3:0.1:1000;
r4 = 30:-0.1:1;
%r5 = 127.0:-0.1:100;
% SOLVING FOR VELOCITY GRADIENT AND MACH NUMBER
[r1,u1] = ode45(@dudr_bondi,r1,i1);
[r2,u2] = ode45(@dudr_bondi,r1,i2);
[r3,u3] = ode45(@dudr_bondi,r2,i1);
[r4,u4] = ode45(@dudr_bondi,r2,i2);
% Parabolas
[r5,u5] = ode45(@dudr_bondi,r1,i3);
[r6,u6] = ode45(@dudr_bondi,r2,i3);
[r7,u7] = ode45(@dudr_bondi,r1,i4);
[r8,u8] = ode45(@dudr_bondi,r2,i4);
[r9,u9] = ode45(@dudr_bondi,r3,0.874*ac);
[r10,u10] = ode45(@dudr_bondi,r3,0.88*ac);
%[r101,u101]= ode45(@ dudr_bondi,rspan5,1.2*ac);
[r11,u11]= ode45(@dudr_bondi,r4,1.365*ac);
[r12,u12]= ode45(@dudr_bondi,r4,1.355*ac);
a1 = sqrt((E + 1./r1 - (u1.^2)./2).*(1/n));
a2 = sqrt((E + 1./r2 - (u2.^2)./2).*(1/n));
a3 = sqrt((E + 1./r3 - (u3.^2)./2).*(1/n));
a4 = sqrt((E + 1./r4 - (u4.^2)./2).*(1/n));
% Parabolas
a5 = sqrt((E + 1./r5 - (u5.^2)./2).*(1/n));
a6 = sqrt((E + 1./r6 - (u6.^2)./2).*(1/n));
a7 = sqrt((E + 1./r7 - (u7.^2)./2).*(1/n));
a8 = sqrt((E + 1./r8 - (u8.^2)./2).*(1/n));
a9 = sqrt((E + 1./r9 - (u9.^2)./2).*(1/n));
a10 = sqrt((E + 1./r10 - (u10.^2)./2).*(1/n));
%a101 = sqrt((E + 1./r101 - (u101.^2)./2).*(1/n));
a11 = sqrt((E + 1./r11 - (u11.^2)./2).*(1/n));
a12 = sqrt((E + 1./r12 - (u12.^2)./2).*(1/n));
% a1 = sqrt((E - phipw(r1) - (u1.^2)./2).*(1/n));
% a2 = sqrt((E - phipw(r2) - (u2.^2)./2).*(1/n));
% a3 = sqrt((E - phipw(r3) - (u3.^2)./2).*(1/n));
% a4 = sqrt((E - phipw(r4) - (u4.^2)./2).*(1/n));
M1 = u1./a1;
M2 = u2./a2;
M3 = u3./a3;
M4 = u4./a4;
% Parabolas
M5 = u5./a5;
M6 = u6./a6;
M7 = u7./a7;
M8 = u8./a8;
M9 = u9./a9;
M10 = u10./a10;
M11 = u11./a11;
M12 = u12./a12;
%M101 = u101./a101;
% PLOTS
semilogx(r1,M1,'red','linewidth',3)
%seimlogx(r1(10),M1(10),’<’,’markersize’,3)
hold on
semilogx(r2,M2,'blue','linewidth',3)
hold on
semilogx(r3,M3,'blue','linewidth',3)
hold on
semilogx(r4,M4,'red','linewidth',3)
hold on
semilogx(r5,M5,'black','linewidth',2)
hold on
semilogx(r6,M6,'black','linewidth',2)
hold on
semilogx(r7,M7,'black','linewidth',2)
hold on
semilogx(r8,M8,'black','linewidth',2)
hold on
semilogx(r9,M9,'black','linewidth',2)
hold on
semilogx(r10,M10,'black','linewidth',2)
hold on
% semilogx(r101,M101,’black’,’linewidth’,2)
% hold on
semilogx(r11,M11,'black','linewidth',2)
hold on
semilogx(r12,M12,'black','linewidth',2)
%grid on
xlabel('r (scaled with r_{sch})')
ylabel('Mach number M')
%FUNCTION FOR INTEGRATION USING ODE45
function dudr = dudr_bondi(r,u,varargin)
E = 0.01;
n = 3;
% PW Plot
% dudr = (dphipw(r) - (2.*((E - u.^2/2 - phipw(r))./n))./r)./(((E - u.^2/2 - phipw(r))./n)./u - u);
% Newtonian Plot
dudr = (1./(r.^2)-(2./(r*n)).*(E - u.^2/2 + 1./r))./((E - u.^2/2 + 1./r)./(u.*n)- u);
end
