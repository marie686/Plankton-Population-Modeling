%same thing with 5 species
max_5=repmat(1, 1, 5);
%mortality rate and D set at .25
mort5=repmat(.25, 5, 1);
%Ks from the article all are intermediate competitors
K5=[ 0.39, 0.34, 0.3, 0.24, 0.23;
 0.22, 0.39, 0.34, 0.3, 0.27;
 0.27, 0.22, 0.39, 0.34, 0.3;
0.3, 0.24, 0.22, 0.39, 0.34;
 0.34, 0.3, 0.22, 0.2, 0.39;];
%What if species 1 is a strong competitor for resource 4?
K6=[ 0.39, 0.34, 0.3, 0.24, 0.23;
 0.22, 0.39, 0.34, 0.3, 0.27;
 0.27, 0.22, 0.39, 0.34, 0.3;
0.2, 0.24, 0.22, 0.39, 0.34;
 0.34, 0.3, 0.22, 0.2, 0.39;];
%What if all species are weak competitors?
K7=[ 0.59, 0.64, 0.73, 0.54, 0.63;
 0.62, 0.79, 0.54, 0.43, 0.47;
 0.47, 0.62, 0.59, 0.84, 0.53;
0.52, 0.74, 0.42, 0.59, 0.54;
 0.54, 0.63, 0.62, 0.52, 0.49;];
%What if there are two strong competitors for the same resource?
K8=[ 0.39, 0.34, 0.3, 0.24, 0.23;
 0.22, 0.15, 0.15, 0.3, 0.27;
 0.27, 0.22, 0.39, 0.34, 0.3;
0.3, 0.24, 0.22, 0.39, 0.34;
 0.34, 0.3, 0.22, 0.2, 0.39;];
%What if they are all strong competitors for one resource?
K9=[ 0.2, 0.12, 0.15, 0.17, 0.13;
 0.22, 0.15, 0.15, 0.3, 0.27;
 0.27, 0.22, 0.39, 0.34, 0.3;
0.3, 0.24, 0.22, 0.39, 0.34;
 0.34, 0.3, 0.22, 0.2, 0.39;];
%Cs also from the article
Content5=[0.04, 0.04, 0.07, 0.04, 0.04;
 0.08, 0.08, 0.08, 0.1, 0.08;
 0.1, 0.1, 0.1, 0.1, 0.14;
 0.05, 0.03, 0.03, 0.03, 0.03;
 0.07, 0.09, 0.07, 0.07, 0.07;];
%over 100 days
tspan=[0,100];
%Supply concentrations
S5=[6, 10, 14, 4, 9];
growth_rate = @(mu, R, col) min(mu .* (R ./ (R + K9(:, col))), [], 'all');
sumfunction = @(row, y) Content5(row,1) * growth_rate(max_5(1), y(6:10), 1) * y(1) + ...
                  Content5(row,2) * growth_rate(max_5(2), y(6:10), 2) * y(2) + ...
                  Content5(row,3) * growth_rate(max_5(3), y(6:10), 3) * y(3) + ...
                    Content5(row,4) * growth_rate(max_5(4), y(6:10), 4) * y(4) + ...
                  Content5(row,5) * growth_rate(max_5(5), y(6:10), 5) * y(5);
f5 = @(t, y) [
    y(1) * (growth_rate(max_5(1), y(6:10), 1) - mort5(1)); 
    y(2) * (growth_rate(max_5(2), y(6:10), 2) - mort5(2)); 
    y(3) * (growth_rate(max_5(3), y(6:10), 3) - mort5(3));
    y(4) * (growth_rate(max_5(4), y(6:10), 4) - mort5(4)); 
    y(5) * (growth_rate(max_5(5), y(6:10), 5) - mort5(5)); 
    .25*(S5(1) - y(6)) - sumfunction(1,y); 
    .25*(S5(2) - y(7)) - sumfunction(2,y); 
    .25*(S5(3) - y(8)) - sumfunction(3,y);
    .25*(S5(4) - y(9)) - sumfunction(4,y); 
    .25*(S5(5) - y(10)) - sumfunction(5,y);
];
figure;
hold on;
for i = 1:size(K5, 1)
    n = 5;
    R_0 = S5;
    N_0 = 0.1 + (1:n) / 100;
    y0 = [N_0, R_0];
    [t, y] = ode45(f5, tspan, y0);
    plot(t, y(:, 1:5));
end
xlabel('Time');
ylabel('Population');
title('Trajectories of the ODE System 5 Species, 5 resources');
grid on;
legend('Species 1', 'Species 2', 'Species 3', 'Species 4', 'Species 5');
hold off;
