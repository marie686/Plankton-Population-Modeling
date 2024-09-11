%for 3 species 3 resources
%max growth rate is 1 
max_growth_rate= repmat(1, 1, 3);

% Methods section stated m=D=.25 per day
mortality_rate=repmat(.25, 3, 1);

%experimentally found K
K =[1, 0.75, 0.25;
 0.25, 1, 0.75;
 0.75, 0.25, 1;];

% cij is resource consumption of resource i by species j stored in matrix Content 
Content=[0.1, 0.2, 0.15;
 0.15, 0.1, 0.2;
 0.2, 0.15, 0.1;];
tspan=[0,100];

%resource supply
S=[10, 10, 10];

growth_rate = @(mu, R, col) min(mu .* (R ./ (R + K(:, col))), [], 'all');


sumfunc = @(row, y) Content(row,1) * growth_rate(max_growth_rate(1), y(4:6), 1) * y(1) + ...
                  Content(row,2) * growth_rate(max_growth_rate(2), y(4:6), 2) * y(2) + ...
                  Content(row,3) * growth_rate(max_growth_rate(3), y(4:6), 3) * y(3);


% Define the ODE system
f = @(t, y) [
    y(1) * (growth_rate(max_growth_rate(1), y(4:6), 1) - mortality_rate(1)); 
    y(2) * (growth_rate(max_growth_rate(2), y(4:6), 2) - mortality_rate(2)); 
    y(3) * (growth_rate(max_growth_rate(3), y(4:6), 3) - mortality_rate(3)); 
    .25 * (S(1) - y(4)) - sumfunc(1, y);
    .25 * (S(2) - y(5)) - sumfunc(2, y);
    .25 * (S(3) - y(6)) - sumfunc(3, y);
];


% Plot each trajectory
figure;
hold on;
for i = 1:size(K, 1)
    n = 3;
    R_0 = S;
    N_0 = 0.1 + (1:n) / 100;
    y0 = [N_0, R_0];
    [t, y] = ode45(f, tspan, y0);
    plot(t, y(:, 1:3));
end

xlabel('Time');
ylabel('Population');
title('Trajectories of the ODE System');
grid on;
legend('Species 1', 'Species 2', 'Species 3');
hold off;

