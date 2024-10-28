clear; close all;

%% Physical Constants and Parameters
hbar = 1.054571817e-34;             % reduced Planck's constant [J*s]
uB = 9.274e-24;                     % Bohr magneton [J/T]
ge = -2.002;                        % electron g-factor
yE = ge*uB/hbar;                    % electron gyromagnetic ratio
B0 = 1.4;                          % magnetic field strength [T]
w0 = -yE*B0;                       % Larmor frequency

%% Simulation Parameters
A = 0.5;                           % field rotation amplitude (0 ≤ A ≤ 1)
Omega = 0.1*abs(w0);              % rotation frequency (test Omega << w0)
N = 1000;                         % number of time steps
dt = 2*pi/(100*abs(w0));          % time step (100 points per Larmor period)
t = (0:N-1)*dt;                   % time array

%% Initial State and Operators
psi0 = [1;0];                     % initial state: |+⟩z
sigmax = [0,1;1,0];              % Pauli matrices
sigmay = [0,-1i;1i,0];
sigmaz = [1,0;0,-1];

%% Magnetic Field Direction
ux = A*cos(Omega*t);
uy = A*sin(Omega*t);
uz = sqrt(1 - A^2)*ones(size(t));

%% Initialize Arrays
U = zeros(2,2,N);
PSI = zeros(2,N);
expect_ox = zeros(1,N);
expect_oy = zeros(1,N);
expect_oz = zeros(1,N);
bloch_norm = zeros(1,N);
psi_norm = zeros(1,N);

%% Time Evolution
U(:,:,1) = eye(2);                % Initialize evolution operator
PSI(:,1) = psi0;
E = hbar*w0/2;                    % Energy scale

% Calculate initial expectation values
expect_ox(1) = real(psi0'*sigmax*psi0);
expect_oy(1) = real(psi0'*sigmay*psi0);
expect_oz(1) = real(psi0'*sigmaz*psi0);
bloch_norm(1) = sqrt(expect_ox(1)^2 + expect_oy(1)^2 + expect_oz(1)^2);
psi_norm(1) = sqrt(abs(psi0'*psi0));

% Time evolution loop
for i = 2:N
    % Calculate Hamiltonian
    H = E*(sigmax*ux(i) + sigmay*uy(i) + sigmaz*uz(i));
    
    % Calculate incremental evolution
    U_temp = expm(-1i/hbar*dt*H);
    
    % Update total evolution
    U(:,:,i) = U_temp * U(:,:,i-1);
    
    % Calculate state
    PSI(:,i) = U(:,:,i)*psi0;
    
    % Calculate expectation values
    expect_ox(i) = real(PSI(:,i)'*sigmax*PSI(:,i));
    expect_oy(i) = real(PSI(:,i)'*sigmay*PSI(:,i));
    expect_oz(i) = real(PSI(:,i)'*sigmaz*PSI(:,i));
    
    % Calculate norms
    bloch_norm(i) = sqrt(expect_ox(i)^2 + expect_oy(i)^2 + expect_oz(i)^2);
    psi_norm(i) = sqrt(abs(PSI(:,i)'*PSI(:,i)));
end

%% Plotting
% figure('Position', [100 100 1200 800]);
figure;
% Expectation values
subplot(2,2,1);
plot(t, expect_ox, 'r-', t, expect_oy, 'g-', t, expect_oz, 'b-', 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('Expectation Value');
title('Pauli Matrix Expectation Values');
legend('\langle\sigma_x\rangle', '\langle\sigma_y\rangle', '\langle\sigma_z\rangle');
grid on;

% Conservation checks
subplot(2,2,2);
plot(t, bloch_norm, 'k-', t, psi_norm, 'r--', 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('Norm');
title('Conservation Checks');
legend('||⟨\vec{\sigma}⟩||', '||\psi||');
ylim([0.99 1.01]);
grid on;

% State probabilities
subplot(2,2,3);
plot(t, abs(PSI(1,:)).^2, 'b-', t, abs(PSI(2,:)).^2, 'r-', 'LineWidth', 1.5);
xlabel('Time (s)');
ylabel('Probability');
title('State Probabilities');
legend('|c_{+,z}|^2', '|c_{-,z}|^2');
grid on;

% 3D Bloch vector trajectory
subplot(2,2,4);

% Plot coordinate axes
axis_length = 1.2;
arrow_size = 0.1;

% Plot arrows
quiver3(0, 0, 0, axis_length, 0, 0, 0, 'r', 'LineWidth', 2, 'MaxHeadSize', arrow_size, 'DisplayName', 'x');
hold on;
quiver3(0, 0, 0, 0, axis_length, 0, 0, 'g', 'LineWidth', 2, 'MaxHeadSize', arrow_size, 'DisplayName', 'y');
quiver3(0, 0, 0, 0, 0, axis_length, 0, 'b', 'LineWidth', 2, 'MaxHeadSize', arrow_size, 'DisplayName', 'z');

% Add axis labels at the tips
text(axis_length+0.1, 0, 0, '|+⟩_x', 'Color', 'r', 'FontSize', 10);
text(-axis_length-0.1, 0, 0, '|-⟩_x', 'Color', 'r', 'FontSize', 10);
text(0, axis_length+0.1, 0, '|+⟩_y', 'Color', 'g', 'FontSize', 10);
text(0, -axis_length-0.1, 0, '|-⟩_y', 'Color', 'g', 'FontSize', 10);
text(0, 0, axis_length+0.1, '|+⟩_z', 'Color', 'b', 'FontSize', 10);
text(0, 0, -axis_length-0.1, '|-⟩_z', 'Color', 'b', 'FontSize', 10);

% Plot Bloch vector trajectory
plot3(expect_ox, expect_oy, expect_oz, 'k-', 'LineWidth', 1, 'DisplayName', 'Trajectory');
plot3(expect_ox(1), expect_oy(1), expect_oz(1), 'go', 'MarkerSize', 10, 'DisplayName', 'Start');
plot3(expect_ox(end), expect_oy(end), expect_oz(end), 'ro', 'MarkerSize', 10, 'DisplayName', 'End');

% Plot Bloch sphere
[X,Y,Z] = sphere(50);
surf(X,Y,Z, 'FaceColor', [0.7 0.7 0.7], 'FaceAlpha', 0.1, ...
    'EdgeColor', [0.7 0.7 0.7], 'EdgeAlpha', 0.2, 'DisplayName', 'Bloch Sphere');

xlabel('\langle\sigma_x\rangle');
ylabel('\langle\sigma_y\rangle');
zlabel('\langle\sigma_z\rangle');
title('Bloch Vector Trajectory');
axis equal;
grid on;
view(45,30);
xlim([-1.2 1.2]);
ylim([-1.2 1.2]);
zlim([-1.2 1.2]);
legend('Location', 'northeastoutside');