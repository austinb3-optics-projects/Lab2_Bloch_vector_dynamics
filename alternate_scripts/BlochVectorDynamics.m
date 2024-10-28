% Main script: BlochVectorDynamics.m
function BlochVectorDynamics
    % Define simulation parameters
    omega0 = 1.0;  % Base frequency
    Omega = 0.1;   % Rotation frequency
    A = 0.5;       % Field amplitude parameter
    n_steps = 1000; % Number of time steps
    total_time = 50.0;
    dt = total_time/n_steps;
    times = linspace(0, total_time, n_steps);
    
    % Pauli matrices
    sigma_x = [0 1; 1 0];
    sigma_y = [0 -1i; 1i 0];
    sigma_z = [1 0; 0 -1];
    
    % Initialize state vectors and measurement arrays
    psi = [1; 0];  % Initial state |+⟩z
    sigma_x_exp = zeros(1, n_steps);
    sigma_y_exp = zeros(1, n_steps);
    sigma_z_exp = zeros(1, n_steps);
    norm_vals = zeros(1, n_steps);
    bloch_vector_norm = zeros(1, n_steps);
    
    % Store initial values
    sigma_x_exp(1) = real(psi' * sigma_x * psi);
    sigma_y_exp(1) = real(psi' * sigma_y * psi);
    sigma_z_exp(1) = real(psi' * sigma_z * psi);
    norm_vals(1) = abs(psi' * psi);
    bloch_vector_norm(1) = sqrt(sigma_x_exp(1)^2 + sigma_y_exp(1)^2 + sigma_z_exp(1)^2);
    
    % Time evolution
    for i = 2:n_steps
        % Calculate Hamiltonian
        H = calculateHamiltonian(times(i-1), omega0, Omega, A, sigma_x, sigma_y, sigma_z);
        
        % Evolution operator
        U = evolutionOperator(H, dt);
        
        % Evolve state
        psi = U * psi;
        
        % Calculate expectation values
        sigma_x_exp(i) = real(psi' * sigma_x * psi);
        sigma_y_exp(i) = real(psi' * sigma_y * psi);
        sigma_z_exp(i) = real(psi' * sigma_z * psi);
        norm_vals(i) = abs(psi' * psi);
        bloch_vector_norm(i) = sqrt(sigma_x_exp(i)^2 + sigma_y_exp(i)^2 + sigma_z_exp(i)^2);
    end
    
    % Plot results
    plotResults(times, sigma_x_exp, sigma_y_exp, sigma_z_exp, norm_vals, bloch_vector_norm);
end

% Function to calculate Hamiltonian
function H = calculateHamiltonian(t, omega0, Omega, A, sigma_x, sigma_y, sigma_z)
    ux = A * cos(Omega * t);
    uy = A * sin(Omega * t);
    uz = sqrt(1 - A^2);
    
    H = 0.5 * omega0 * (ux * sigma_x + uy * sigma_y + uz * sigma_z);
end

% Function to calculate evolution operator
function U = evolutionOperator(H, dt)
    % Second-order expansion of evolution operator
    U = eye(2) - 1i*dt/pi * H - (dt^2/(2*pi^2)) * (H * H);
end

% Function to plot results
function plotResults(times, sigma_x_exp, sigma_y_exp, sigma_z_exp, norm_vals, bloch_vector_norm)
    figure('Position', [100 100 1200 800]);
    
    % Plot expectation values
    subplot(2,2,1);
    plot(times, sigma_x_exp, 'r-', 'LineWidth', 1.5);
    hold on;
    plot(times, sigma_y_exp, 'g-', 'LineWidth', 1.5);
    plot(times, sigma_z_exp, 'b-', 'LineWidth', 1.5);
    xlabel('Time');
    ylabel('Expectation Value');
    title('Spin Components');
    legend('⟨σx⟩', '⟨σy⟩', '⟨σz⟩');
    grid on;
    
    % Plot state norm
    subplot(2,2,2);
    plot(times, norm_vals, 'k-', 'LineWidth', 1.5);
    xlabel('Time');
    ylabel('State Norm');
    title('State Normalization');
    grid on;
    
    % Plot Bloch vector norm
    subplot(2,2,3);
    plot(times, bloch_vector_norm, 'b-', 'LineWidth', 1.5);
    xlabel('Time');
    ylabel('Bloch Vector Norm');
    title('Bloch Vector Normalization');
    grid on;
    
    % 3D trajectory plot
    subplot(2,2,4);
    plot3(sigma_x_exp, sigma_y_exp, sigma_z_exp, 'b-', 'LineWidth', 1.5);
    grid on;
    xlabel('X');
    ylabel('Y');
    zlabel('Z');
    title('Bloch Vector Trajectory');
    view(45, 30);
    
    % Add Bloch sphere wireframe
    hold on;
    [X,Y,Z] = sphere(50);
    surf(X,Y,Z, 'FaceAlpha', 0.1, 'EdgeAlpha', 0.1);
    axis equal;
end

% Script to run multiple test cases
function runTestCases()
    % Case 1: Static field (A = 0)
    disp('Running simulation with static field (A = 0)...');
    omega0 = 1.0;
    Omega = 0.1;
    A = 0.0;
    BlochVectorDynamics(omega0, Omega, A);
    
    % Case 2: Rotating field with Ω << ω0
    disp('Running simulation with Ω << ω0...');
    omega0 = 10.0;
    Omega = 0.1;
    A = 0.5;
    BlochVectorDynamics(omega0, Omega, A);
    
    % Case 3: Rotating field with Ω >> ω0
    disp('Running simulation with Ω >> ω0...');
    omega0 = 0.1;
    Omega = 10.0;
    A = 0.5;
    BlochVectorDynamics(omega0, Omega, A);
end