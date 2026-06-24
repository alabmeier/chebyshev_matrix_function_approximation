clear
clc
close all

%% Part 1 - Relative Error Analysis

% Parameters for Chebyshev approximation
a = 0;
ng = 1000;
t = 0.001;

% Test for different values of N and C
N_values = 100:100:1000;
C_values = [0.01 0.1 0.25 0.5];

relative_errors = zeros(length(C_values), length(N_values));

for i = 1:length(C_values)

    C = C_values(i);

    for j = 1:length(N_values)

        N = N_values(j);

        % Set right endpoint
        b = 4*(N^2);

        % Create the sparse matrix A
        main_diag = 2*ones(N,1);
        upper_lower_diag = -1*ones(N,1);

        A = spdiags([upper_lower_diag main_diag upper_lower_diag], ...
                    [-1 0 1], N, N);

        A(1,N) = -1;
        A(N,1) = -1;

        % Create function handle
        f = @(x) exp(-t*x);

        % Create corresponding initial condition
        v = rand(N,1);
        v = v/norm(v);

        % Compute ground truth using expm
        w_true = expm(-t*A) * v;

        % Compute approximation using Chebyshev
        n = round(C * N);
        c_approx = cheb_coeff(f, a, b, n, ng);
        w_approx = cheb_matvec(A, v, a, b, c_approx);

        % Calculate relative error
        relative_errors(i,j) = norm(w_approx - w_true) / norm(w_true);

    end

end

% Plot results
figure;
semilogy(N_values, relative_errors, 'o-', 'LineWidth', 1.5);
hold on;
yline(1e-8, 'r--', 'Threshold 10^{-8}', 'LineWidth', 1.5);
xlabel('Problem size N');
ylabel('Relative Error');
grid on;
legend([arrayfun(@(x) sprintf('C=%.1f', x), C_values, ...
       'UniformOutput', false), {'Threshold'}]);
title('Relative Error vs Problem Size for Different Values of C');


%% Part 2 - Computational Runtime Analysis

% Parameters for Chebyshev approximation
a = 0;
t = 0.001;
ng = 1000;

% Test for different values of N
N_values = 500:500:3000;
times = zeros(size(N_values));

for i = 1:length(N_values)

    N = N_values(i);
    n = round(0.5*N);

    % Set right endpoint
    b = 4*(N^2);

    % Create the sparse matrix A
    main_diag = 2*ones(N,1);
    upper_lower_diag = -1*ones(N,1);

    A = spdiags([upper_lower_diag main_diag upper_lower_diag], ...
                [-1 0 1], N, N);

    A(1,N) = -1;
    A(N,1) = -1;

    % Create function handle
    f = @(x) exp(-t*x);

    % Create corresponding initial condition
    v = rand(N,1);
    v = v/norm(v);

    % Measure time for Chebyshev method
    tic;
    c_approx = cheb_coeff(f, a, b, n, ng);
    w_approx = cheb_matvec(A, v, a, b, c_approx);
    times(i) = toc;

end

% Plot results
figure;
plot(N_values, times, 'LineWidth', 1.5);
hold on;

% Validate with quadratic curve
p = polyfit(N_values, times, 2);
N_fit = linspace(min(N_values), max(N_values), 100);
times_fit = polyval(p, N_fit);

plot(N_fit, times_fit, 'LineWidth', 1.5);

xlabel('N');
ylabel('Computational time in seconds');
title('Computational Time vs N for Chebyshev Method');
grid on;
legend('Actual times', 'Quadratic fit');
