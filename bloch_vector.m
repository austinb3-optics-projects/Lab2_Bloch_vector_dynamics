clear; close all;

%% Initialization 
hbar = 1.054571817e-34;             % reduced plancks constant [J*s]
uB = 9.274e-24;                     % Bohr magnetron [J/T]
uN = 5.051e-27;                     % Nuclear magneton [J/T]
ge = -2.002; gp = 5.586; gn = -3.826;       % g factors 
yE =ge*uB/hbar; yP = gp*uN/hbar; yN = gn*uN/hbar;   % gyromagnetic ratios 
B0 = 1.4;           % magnetic field strength of a neodynium magnet [T]
wE_L = -yE*B0; wP_L = -yP*B0; wN_L = -yN*B0;   % larmor frequency

theta = 0; phi = 0;     % spherical coordinate system parameters
% u = [sin(theta).*cos(phi); sin(theta)*sin(phi); cos(theta)];    % arbitrary spin direction vector
sigmax = [0,1;1,0]; sigmay = [0,-1i;1i,0]; sigmaz = [1,0;0,-1];     % pauli matrices these are in the z basis
% Sx = hbar/2*sigmax; Sy = hbar/2*sigmay; Sz = hbar/2*sigmaz;         % spin operators
% S = [Sx;Sy;Sz];             % spin operator
% sigma = 2/hbar.*S;

psi0 = [1;0];           % initial state: Spin up along z

Omega = wE_L;           % generalized rabi frequency.

% setup of time evolution operator?
N = 1000;       % total time steps?
dt = 1e-16;
t_end = N*dt;      % let final time be 1 second.
t = 0:dt:t_end;
% looking at values of dt to where the code breaks down based off of |psi|
% dt = 1e-4 => no good.
% so long as dt is 1e-15 or smaller, the assumptions are valid.
% I think that this time step is predicated


% magnetic field vector direction 
A = 0.5;
ux = A.*cos(Omega.*t);
uy = A.*sin(Omega.*t);
uz = sqrt(1 - A.^2);

% setup for loop over time to determine U(t,t0), H(t), psi(t) based on u(t)
U = ones(2,2,length(t));  % time evolution operator
U_temp = ones(2,2,length(t));
PSI = ones(length(psi0),length(t));
expect_ox = zeros(1,length(t));
expect_oy = zeros(1,length(t));
expect_oz = zeros(1,length(t));
sigma = zeros(1,length(t));
psi = zeros(1,length(t));
E = hbar*wE_L./2;       % prefactor for Hamiltonian

%% Computation

for i = 1:length(t)

% perform dot product 
    dotproduct = sigmax.*ux(i) + sigmay.*uy(i) + sigmaz.*uz;
    H = E.*dotproduct;
% use prod() to perform product operation for U(t) at each time step
% make sure to use expm for e^A where A is a matrix

    U_temp(:,:,i) = expm(-1i./hbar.*dt.*H);
% U_temp(:,:,i)
    U(:,:,i) = prod(U_temp,3);      % perform the product over the third dimension 
                                % which ends up as elementwise multiplication
                                % between the elements of the 2x2 matrices
                                % stored in U_temp
    PSI(:,i) = U(:,:,i)*psi0;
    expect_ox(i) = PSI(:,i)'*sigmax*PSI(:,i);
    expect_oy(i) = PSI(:,i)'*sigmay*PSI(:,i);
    expect_oz(i) = PSI(:,i)'*sigmaz*PSI(:,i);
    
    sigma(i) = sqrt(expect_oz(i)^2 + expect_oy(i)^2 + expect_ox(i)^2);
    psi(i) = sqrt(PSI(1,i)^2 + PSI(2,i)^2);
end

% A = [1,2;3,4];
% B = [2,3;4,5];
% C(:,:,1) = A;
% C(:,:,2) = B;
% C
% D = prod(C,3)

%% Plotting

figure;
plot(t,expect_ox,'LineWidth',1);
xlabel('$t$','Interpreter','latex');
ylabel('$\langle \sigma_x \rangle(t)$','Interpreter','latex');
title('Expectation value along $x$','Interpreter','latex');
set(gca,'FontSize',15);

figure;
plot(t,expect_oy,'LineWidth',1);
xlabel('$t$','Interpreter','latex');
ylabel('$\langle \sigma_y \rangle(t)$','Interpreter','latex');
title('Expectation value along $y$','Interpreter','latex');
set(gca,'FontSize',15);

figure;
plot(t,expect_oz,'LineWidth',1);
xlabel('$t$','Interpreter','latex');
ylabel('$\langle \sigma_z \rangle(t)$','Interpreter','latex');
title('Expectation value along $z$','Interpreter','latex');
set(gca,'FontSize',15);

figure;
plot(t,abs(PSI(1,:)).^2,'LineWidth',1.5); hold on;
plot(t,abs(PSI(2,:)).^2,'LineWidth',1.5); hold off;
xlabel('$t$','Interpreter','latex');
ylabel('$|c_{\pm,z}|^2$','Interpreter','latex');
title('Expansion coefficients','Interpreter','latex');
legend('|c_{+,z}|','|c_{-,z}|');
set(gca,'FontSize',15);