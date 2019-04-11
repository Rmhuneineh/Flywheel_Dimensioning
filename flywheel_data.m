%% Engine Characteristic
ic = 6; % [-] number of cylinders
d = 72.5*1e-3; % [m] bore
c = 66*1e-3; % [m] stroke
epsilon = 9.4; % [-] compression ratio
lambda = 0.25; % [-] crank slider parameter
n = 5400; % [rpm] engine speed
AM = 1.2*1e3; % [kg/m^3] alternative masses
ERM = 0.85*1e3; % [kg/m^3] equivalent rotating masses

%% Data for point 1 calculation
Pa = 100658; % [Pa] environment pressure
Ta = 297; % [K] environment temperature
lambdav = 0.82; % [-] volumetric efficiency
Pr_Pa = 1.15; % [-] residual to environment pressure ratio
Tr = 903; % [K] residual temperature
alpha_st = 14.7; % [-] stoichiometric air-to-fuel ratio
R_alpha = 0.95; % [-] relative air-to-fuel ratio
Cb = 2500; % [J/(Kg.K)] fuel specific heat
x = 1; % [-] fuel vaporized fraction
r = 320000; % [J/Kg] fuel vaporization heat
delta_T = 30; % [degC] temperature increment for the air-fuel-residual mixture during the intake stroke
R = 287.2; % [J/(Kg.K)] air elastic constant
Rp = 288; % [J/(Kg.K)] burnt gas elastic constant
R1 = 271; % [J/[Kg.K)] air-fuel-residual mixture elastic constant
Cp = 1009; % [J/(Kg.K)] air specific heat
Cpp = 1150; % [J/Kg.K)] burnt gas specific heat

%% Data for point 2 calculation
m = 1.35; % [-] comression index
% To calculate lower heating value at T2
Cv = 744; % [J/(Kg.K)] air specific heat
Cvp = 824; % [J/(Kg.K)] burnt gas specific heat
Hiv_T0 = 44000000; % [J/Kg] lower heating value at T0 = 288 K
T0 = 288; % [K]

%% Data for point 3 calculation
% Dissciation heat
dq = 0.54; % [J/(Kg.K^2)]
Td = 1850; % [K]
Cpp3 = 1328; % [J/(Kg.K)] burnt gas specific heat
delta_A = 0.06; % [-] heat loss coefficient

%% Data for point 4 calculation
mp = 1.27; % [-] expansion index

%% Data for flywheel design
Gamma_v = 7.7*1e3; % [Kg/m^3] density
delta = 0.01; %  [-] kinematic irregularity
