%Quanser Aero 2DoF Model

%DC Motor Model
Rm=8.4; % Terminal resistance  (?)
kt=0.042;   % Torque constant (N ? m/A)
km=0.042;   % Motor back-emf constant (V/(rad/s))
Lm=1.16e-3; % Rotor inductance (H)


%Propeller
kd=1e-5;    % Drag/Air resistance coefficient (N ? m/(rad/s))
Jh=3.04e-9; % Propeller hub inertia (kg ? m2)
Jprop=7.2e-6;  % Propeller inertia (kg ? m2)
Jeq=Jprop+Jh;

%DC Motor Static Model
G=kt/(Rm*kd+kt*km); %Static gain of transfer function ((rad/s)/V)

% Moment of Inertia of helicopter body (kg-m^2)
L_body = 6.5*0.0254; % length of horizontal body (metal tube)
m_body = 0.094; % mass of horizontal body (metal tube) 
J_body = m_body * L_body^2 / 12; % horizontal cylinder rotating about CM

% Moment of Inertia of yoke fork that rotates about yaw axis (kg-m^2)
m_yoke = 0.526; % mass of entire yoke assembly (kg)
% h_yoke = 9*0.0254; % height of yoke assembly (m)
r_fork = 0.04/2; % radius of each fork (approximated as cylinder)
J_yoke = 0.5*m_yoke*r_fork^2;

% Moment of Inertia from motor + guard assembly about pivot (kg-m^2)
m_prop = 0.43; % mass of dc motor + shield + propeller shield
% m_motor = 0.203; % mass of dc motor
r_prop = 6.25*0.0254; % distance from CM to center of pitch axis
J_prop = m_prop * r_prop^2; % using parallel axis theorem

% Equivalent Moment of Inertia about Pitch and Yaw Axis (kg-m^2)
Jp = J_body + 2*J_prop; % pitch: body and 2 props
Jy = J_body + 2*J_prop + J_yoke; % yaw: body, 2 props, and yoke

% Voltage-torque constant (N-m/V) [found experimentally]
Kpp = 0.0011; % (pre-production unit: 0.0015)
Kyy = 0.0022; % (pre-production unit: 0.0040)
Kpy = 0.0021; % thrust acting on pitch from yaw (pre-production unit: 0.0020)
Kyp = -0.0027; % thrust acting on yaw from pitch (pre-production unit: -0.0017)

%Rotation-torque constant (Nm/(rad/s))
Kppw = Kpp/G; 
Kyyw = Kyy/G; 
Kpyw = Kpy/G; 
Kypw = -Kyp/G; 

% Stiffness (N-m/rad)[found experimentally] 
Ksp = 0.037463;

% Viscous damping (N-m-s/rad) [found experimentally]
Dp = 0.0071116; % pitch axis (pre-production unit: Dp = 0.0226) 
Dy = 0.0220; % yaw axis (pre-production unit: Dy = 0.0211)
% 