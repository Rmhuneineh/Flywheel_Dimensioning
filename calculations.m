%% Clear ALL
clear all
close all
clc

%% Load Data
flywheel_data

%% P1
Vd = pi*d^2*c/4; % [m^3] displacement volume
Vc = Vd/(epsilon-1); % [m^3] clearance volume
Vt = Vd + Vc; % [m^3] total cylinder volume
Pr = Pa*Pr_Pa; % [Pa] % residual gas pressure
alpha = alpha_st*(R_alpha); % [-] air-to-fuel ratio
alphap = Pr*Vc*R*Ta/Pa/Vt/Rp/Tr/lambdav; % [-] residual-to-fuel ratio

T1 = ((alpha*Cp + Cb)*Ta + alphap*Cpp*Tr - x*r)/(alpha*Cp + Cb + alphap*Cpp) + delta_T; % [K]
P1 = Pa*R*T1/epsilon*(lambdav*((1+alpha)/alpha)*(epsilon-1)/R/Ta + Pr/Pa/Rp/Tr); % [Pa]

%% P2
P2 = P1*epsilon^m; % [Pa]
T2 = T1*epsilon^(m-1); % [K]

%% P3
Hiv_T2 = Hiv_T0 + (1+alpha)*(Cv-Cvp)*(T2-T0); % [J/Kg]
A = dq;
B = Cvp-2*dq*Td;
C = dq*Td^2 - Cvp*T2 - (1-delta_A)*R_alpha*Hiv_T2/(1+alpha+alphap);
ROOTS = roots([A B C]);
T3 = ROOTS(ROOTS>T2); % [K]
P3 = P2*Rp*T3/R1/T2; % [Pa]

%% P4
P4 = P3*epsilon^-mp; % [Pa]
T4 = T3*epsilon^(1-mp); % [K]

%% P5
P5 = Pr; % [Pa]

%% P6
P6 = Pr; % [Pa]

%% P7
P7 = P1; % [Pa]
%% Crank-slider mechanism
theta = [0:180 180:360 360:540 540:720]; % [deg]
c_r = c/2; % [m] crank radius
area = pi*d^2/4; % [m^2] cylinder cross-sectional area
x_p = c_r*((1-cosd(theta)) + 1/lambda*(1-sqrt(1-lambda^2*sind(theta).^2))) + Vc/area; % [m] piston position
V = area*x_p; % [m^3] cylinder volume vs theta
%% Pressure vs Volume
% Expansion (3 -> 4)
V34 = V(1:181);
P(1:181) = P3*(Vc./V34).^mp;

% Blowdown (4 -> 5)
P(182) = P5;

% Exhaust (5 -> 6)
P(183:362) = P6;

% Intake (6 -> 7 -> 1)
P(363) = P7;
P(364:543) = P1;

% Compression (1 -> 2)
V12 = V(544:723);
P(544:723) = P1*(Vt./V12).^m;

% Combustion (2 -> 3)
P(724) = P3;

figure(1)
hold on
plot(V*1e6, P*1e-5)
title('Pressure vs. Volume')
xlabel('Volume [cm^3]')
ylabel('Pressure [bar]')
grid on
%% Inertia Pressure
beta = asind(lambda*sind(theta)); % [deg]
omega = n*2*pi/60; % [rad/s]
Fin = -AM*Vt*omega^2*c_r*(cosd(theta) + lambda*cosd(2*theta)./cosd(beta)); % [N] Inertia Force
Pin = Fin/area; % [Pa] Inertia Pressure

plot(V*1e6, -Pin*1e-5, 'r')
legend('Indicated Cycle', 'Inertia Pressure')
%% Effective Pressure
Peff = P - Pa + Pin; % [Pa] effective pressure
figure(2)
plot(theta(1:end-1), Peff(1:end-1)*1e-5)
title('Effective Pressure vs Crank Angle')
xlabel('\theta [deg]')
ylabel('P_e_f_f [bar]')
grid on

%% Tangential Tension
gamma = theta + beta; % [deg]
F = Peff*area./cosd(beta); % [N] force along connecting rod
Ft = F.*sind(gamma); % [N] tangential force
Ms = Ft*c_r; % [Nm] shaft moment
t = Ms/(Vd/2); % [Pa] tangential tension
figure(3)
plot(theta, t*1e-5)
title('Tangential Tension vs Crank Angle')
xlabel('\theta [deg]')
ylabel('t [bar]')
grid on

%% Indicated Mean Effective Pressure
imep = P3*(epsilon^(1-mp)-1)/(epsilon-1)/(1-mp) + P1 - Pr + ...
    epsilon*P1*(epsilon^(m-1)-1)/(epsilon-1)/(1-m);
%% Resistive Moment and Tension
Mr = imep*Vd/4/pi*ones(1, 724); % [Nm] resistive moment
tr = Mr/(Vd/2); % [Pa] resitive tension
figure(4)
plot(theta(1:181), t(1:181)*1e-5, theta(1:181), tr(1:181)*1e-5)
title('EXPANSION STROKE: Tension vs Crank Angle')
xlabel('\theta [deg]')
ylabel('[bar]')
legend('Tangential Tension', 'Resistant Tension')
grid on

%% Dynamic Irregularity
omega_max = omega*(2+delta)/2; % [rad/s] maximum angular velocity
omega_min = omega_max - delta*omega; % [rad/s] minimum angular velocity
THETAS = InterX([theta(1:181); t(1:181)], [theta(1:181); tr(1:181)]);
theta_min = THETAS(1, 1); % [deg] theta @ omega_min
theta_max = THETAS(1, 2); % [deg] theta @ omega_max
Ls_min = trapz(theta(theta<=theta_min)*pi/180, Ms(theta<=theta_min)); % [J] 
Lr_min = trapz(theta(theta<=theta_min)*pi/180, Mr(theta<=theta_min)); % [J]
deltaL_min = Ls_min - Lr_min; % [J]
Ls_max = trapz(theta(theta<=theta_max)*pi/180, Ms(theta<=theta_max)); % [J]
Lr_max = trapz(theta(theta<=theta_max)*pi/180, Mr(theta<=theta_max)); % [J]
deltaL_max = Ls_max - Lr_max; % [J]
deltaL_t = deltaL_max + abs(deltaL_min); % [J]
zeta = deltaL_t/imep/Vd; % [-] dynamic irregularity

%% Flywheel
J = zeta*imep*Vd/delta/omega^2; % [kg.m^2] moment of inertia
J_eng = ERM*ic*Vd*c_r^2; % [Kg.m^2] engine moment of inertia
J_fly = J - J_eng; % [Kg.m^2] flywheel moment of inertia
D = (J_fly*320/pi/Gamma_v)^(1/5); % [m] flywheel diameter

%% Shaft Instantaneous Angular Velocity
Ls(1) = 0;
Lr(1) = 0;
for i = 2:length(theta)
    Ls(i) = trapz([theta(i-1) theta(i)]*pi/180, [Ms(i-1) Ms(i)])+Ls(i-1);
    Lr(i) = trapz([theta(1) theta(2)]*pi/180, [Mr(1) Mr(2)])+Lr(i-1);
end
figure(5)
plot(theta, Ls)
hold on
plot(theta, Lr)
title('Work vs Crank Angle')
xlabel('\theta [deg]')
ylabel('L [J]')
legend('Engine Work', 'Resistant Work')
grid on

omega_i = sqrt(omega^2 + 2/J*(Ls-Lr));
shift = mean(omega_i) - omega;
omega_i = omega_i - shift;

figure(6)
plot(theta, omega_i)
title('Shaft Angular Speed vs Crank Angle')
xlabel('\theta [deg]')
ylabel('\omega_s [rad/s]')
grid on

%% Multi-cylinder Engine
delta_phi = 720/ic; % [deg] Phase Shift
% Cylinder: 1
Ms_c1 = Ms; % [Nm]
t_c1 = t; % [Pa]
Mr_c1 = Mr; % [Nm]
tr_c1 = Mr_c1/(Vd/2); % [Pa]

% Cylinder 2
Ms_c2(1:delta_phi+1) = Ms(end-delta_phi:end); % [Nm]
Ms_c2(delta_phi+2:724) = Ms(1:end-delta_phi-1); % [Nm]
t_c2(1:delta_phi+1) = t(end-delta_phi:end);% [Pa]
t_c2(delta_phi+2:724) = t(1:end-delta_phi-1); % [Pa]
Mr_c2(1:delta_phi+1) = Mr(end-delta_phi:end); % [Nm]
Mr_c2(delta_phi+2:724) = Mr(1:end-delta_phi-1); % [Nm]
tr_c2 = Mr_c2/(Vd/2); % [Pa]

% Cylinder 3
Ms_c3(1:2*(delta_phi+1)) = Ms(end-2*(delta_phi+1)+1:end); % [Nm]
Ms_c3(2*(delta_phi+1)+1:724) = Ms(1:end-2*(delta_phi+1)); % [Nm]
t_c3(1:2*(delta_phi+1)) = t(end-2*(delta_phi+1)+1:end);% [Pa]
t_c3(2*(delta_phi+1)+1:724) = t(1:end-2*(delta_phi+1)); % [Pa]
Mr_c3(1:2*(delta_phi+1)) = Mr(end-2*(delta_phi+1)+1:end); % [Nm]
Mr_c3(2*(delta_phi+1)+1:724) = Fin(1:end-2*(delta_phi+1)); % [Nm]
tr_c3 = Mr_c3/(Vd/2); % [Pa]

% Cylinder 4
Ms_c4(1:3*(delta_phi+1)) = Ms(end-3*(delta_phi+1)+1:end); % [Nm]
Ms_c4(3*(delta_phi+1)+1:724) = Ms(1:end-3*(delta_phi+1)); % [Nm]
t_c4(1:3*(delta_phi+1)) = t(end-3*(delta_phi+1)+1:end);% [Pa]
t_c4(3*(delta_phi+1)+1:724) = t(1:end-3*(delta_phi+1)); % [Pa]
Mr_c4(1:3*(delta_phi+1)) = Mr(end-3*(delta_phi+1)+1:end); % [Nm]
Mr_c4(3*(delta_phi+1)+1:724) = Mr(1:end-3*(delta_phi+1)); % [Nm]
tr_c4 = Mr_c4/(Vd/2); % [Pa]

% Cylinder 5
Ms_c5(1:4*(delta_phi+1)) = Ms(end-4*(delta_phi+1)+1:end); % [Nm]
Ms_c5(4*(delta_phi+1)+1:724) = Ms(1:end-4*(delta_phi+1)); % [Nm]
t_c5(1:4*(delta_phi+1)) = t(end-4*(delta_phi+1)+1:end);% [Pa]
t_c5(4*(delta_phi+1)+1:724) = t(1:end-4*(delta_phi+1)); % [deg]
Mr_c5(1:4*(delta_phi+1)) = Mr(end-4*(delta_phi+1)+1:end); % [Nm]
Mr_c5(4*(delta_phi+1)+1:724) = Mr(1:end-4*(delta_phi+1)); % [Nm]
tr_c5 = Mr_c5/(Vd/2); % [Pa]

% Cylinder 6
Ms_c6(1:5*(delta_phi+1)) = Ms(end-5*(delta_phi+1)+1:end); % [Nm]
Ms_c6(5*(delta_phi+1)+1:724) = Ms(1:end-5*(delta_phi+1)); % [Nm]
t_c6(1:5*(delta_phi+1)) = t(end-5*(delta_phi+1)+1:end);% [Pa]
t_c6(5*(delta_phi+1)+1:724) = t(1:end-5*(delta_phi+1)); % [Pa]
Mr_c6(1:5*(delta_phi+1)) = Mr(end-5*(delta_phi+1)+1:end); % [Nm]
Mr_c6(5*(delta_phi+1)+1:724) = Mr(1:end-5*(delta_phi+1)); % [Nm]
tr_c6 = Mr_c6/(Vd/2); % [Pa]

% Engine
Ms_eng = Ms_c1 + Ms_c2 + Ms_c3 + Ms_c4 + Ms_c5 + Ms_c6;
t_eng = Ms_eng/(Vd/2);
Mr_eng = imep*ic*Vd/4/pi*ones(1, 724);
tr_eng = Mr_eng/(Vd/2);
figure(7)
plot(theta, t_c1*1e-5, theta, t_c2*1e-5, theta, t_c3*1e-5,...
    theta, t_c4*1e-5, theta, t_c5*1e-5, theta, t_c6*1e-5,...
    theta, t_eng*1e-5)
title('Tangential Tension')
xlabel('\theta [deg]')
ylabel('[bar]')
grid on
legend('Cyl 1', 'Cyl 2', 'Cyl 3', 'Cyl 4', 'Cyl 5', 'Cyl 6', 'multi-cyl')

%% Dynamic Irregularity
THETAS_eng = InterX([theta(1:121); t_eng(1:121)], [theta(1:121); tr_eng(1:121)]);
theta_min_eng = THETAS_eng(1, 1); % [deg]
theta_max_eng = THETAS_eng(1, end); % [deg]
Ls_min_eng = trapz(theta(theta<=theta_min_eng)*pi/180, Ms_eng(theta<=theta_min_eng)); % [J] 
Lr_min_eng = trapz(theta(theta<=theta_min_eng)*pi/180, Mr_eng(theta<=theta_min_eng)); % [J]
deltaL_min_eng = Ls_min_eng - Lr_min_eng; % [J]
Ls_max_eng = trapz(theta(theta<=theta_max_eng)*pi/180, Ms(theta<=theta_max_eng)); % [J]
Lr_max_eng = trapz(theta(theta<=theta_max_eng)*pi/180, Mr(theta<=theta_max_eng)); % [J]
deltaL_max_eng = Ls_max_eng - Lr_max_eng; % [J]
deltaL_t_eng = deltaL_max_eng + abs(deltaL_min_eng); % [J]
zeta_eng = deltaL_t_eng/imep/Vd/ic; % [-] dynamic irregularity

%% Engine Flywheel
J1 = zeta_eng*imep*Vd/delta/omega^2; % [kg.m^2] moment of inertia
J_eng1 = ERM*Vd*ic*c_r^2; % [Kg.m^2] engine moment of inertia
J_fly_eng = J1 - J_eng1; % [Kg.m^2] flywheel moment of inertia
D1 = (J_fly_eng*320/pi/Gamma_v)^(1/5); % [m] flywheel diameter

%% Work
Ls_eng(1) = 0;
Lr_eng(1) = 0;
for i = 2:length(theta)
    Ls_eng(i) = trapz([theta(i-1) theta(i)]*pi/180, [Ms_eng(i-1) Ms_eng(i)])+Ls_eng(i-1);
    Lr_eng(i) = trapz([theta(i-1) theta(i)]*pi/180, [Mr_eng(i-1) Mr_eng(i)])+Lr_eng(i-1);
end
figure(8)
plot(theta(1:121), Ls_eng(1:121)/Vd*1e-5)
hold on
plot(theta(1:121), Lr_eng(1:121)/Vd*1e-5)
plot(theta(1:121), t_eng(1:121)*1e-5)
plot(theta(1:121), tr_eng(1:121)*1e-5)
title('Work vs Crank Angle')
xlabel('\theta [deg]')
ylabel('L [J]')
legend('Engine Work', 'Resistant Work', 'Tangential Tension', 'Resistant Tension')
grid on

omega_i_eng = sqrt(omega^2 + 2/J*(Ls_eng-Lr_eng));
shift_eng = mean(omega_i_eng) - omega;
omega_i_eng = omega_i_eng - shift_eng;

figure(9)
plot(theta(1:121), omega_i_eng(1:121))
title('Shaft Angular Speed vs Crank Angle')
xlabel('\theta [deg]')
ylabel('\omega_s [rad/s]')
grid on