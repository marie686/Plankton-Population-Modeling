%for 9 species 3 resources
%max growth rate is 1 
max_growth_rate2= repmat(1, 1, 9);

% Methods section stated m=D=.25 per day
mortality_rate=repmat(.25, 9, 1);


%experimentally found K
Content9 =[ 0.1, 0.2, 0.15, 0.05, 0.01, 0.4, 0.3, 0.2, 0.25;
 0.15, 0.1, 0.2, 0.15, 0.3, 0.35, 0.25, 0.02, 0.35,
 0.2, 0.15, 0.1, 0.25, 0.05, 0.2, 0.4, 0.15, 0.1];

% cij is resource consumption of resource i by species j stored in matrix Content 
Ke=[1, 0.75, 0.25, 0.7, 0.2, 0.65, 0.68, 0.38, 0.46;
 0.25, 1, 0.75, 0.2, 1.01, 0.55, 0.83, 1.1, 0.85;
 0.75, 0.25, 1, 1.1, 0.7, 0.95, 0.6, 0.5, 0.77;];
tspan=[0,100];

%resource supply
%What if Supply is varied?
S=[50, 10, 10];
growth_rate = @(mu, R, col) min(mu .* (R ./ (R + Ke(:, col))), [], 'all');


sumfuncT = @(row, y) Content9(row,1) * growth_rate(max_growth_rate2(1), y(10:12), 1) * y(1) + ...
                  Content9(row,2) * growth_rate(max_growth_rate2(2), y(10:12), 2) * y(2) + ...
                  Content9(row,3) * growth_rate(max_growth_rate2(3), y(10:12), 3) * y(3)+ ...
                  Content9(row,4) * growth_rate(max_growth_rate2(4), y(10:12), 4) * y(4) + ...
                  Content9(row,5) * growth_rate(max_growth_rate2(5), y(10:12), 5) * y(5) + ...
                  Content9(row,6) * growth_rate(max_growth_rate2(6), y(10:12), 6) * y(6)+ ...
                  Content9(row,7) * growth_rate(max_growth_rate2(7), y(10:12), 7) * y(7) + ...
                  Content9(row,8) * growth_rate(max_growth_rate2(8), y(10:12), 8) * y(8) + ...
                  Content9(row,9) * growth_rate(max_growth_rate2(9), y(10:12), 9) * y(9);


% Define the ODE system
f = @(t, y) [
    y(1) * (growth_rate(max_growth_rate2(1), y(10:12), 1) - mortality_rate(1)); 
    y(2) * (growth_rate(max_growth_rate2(2), y(10:12), 2) - mortality_rate(2)); 
    y(3) * (growth_rate(max_growth_rate2(3), y(10:12), 3) - mortality_rate(3));
    y(4) * (growth_rate(max_growth_rate2(1), y(10:12), 4) - mortality_rate(4)); 
    y(5) * (growth_rate(max_growth_rate2(2), y(10:12), 5) - mortality_rate(5)); 
    y(6) * (growth_rate(max_growth_rate2(3), y(10:12), 6) - mortality_rate(6));
    y(7) * (growth_rate(max_growth_rate2(1), y(10:12), 7) - mortality_rate(7)); 
    y(8) * (growth_rate(max_growth_rate2(2), y(10:12), 8) - mortality_rate(8)); 
    y(9) * (growth_rate(max_growth_rate2(3), y(10:12), 9) - mortality_rate(9));
    .25 * (S(1) - y(4)) - sumfuncT(1, y);
    .25 * (S(2) - y(5)) - sumfuncT(2, y);
    .25 * (S(3) - y(6)) - sumfuncT(3, y);
];
% Plot each trajectory
figure;
hold on;
for i = 1:size(Ke, 1)
    n = 9;
    R_0 = S;
    %what is N_0 is varied?
    N_0 = 100 + (1:n) / 100;
    y0 = [N_0, R_0];
    [t, y] = ode45(f, tspan, y0);
    plot(t, y(:, 1:9));
end

xlabel('Time');
ylabel('Population');
title('Trajectories of the ODE System 9 Species 3 Resources');
grid on;
legend('Species 1', 'Species 2', 'Species 3', 'Species 4', 'Species 5', 'Species 6', 'Species 7', 'Species 8', 'Species 9');
hold off;
