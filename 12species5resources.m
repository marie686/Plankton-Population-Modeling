%same thing with 12 species
max_50=repmat(1, 1, 12);
%mortality rate and D set at .25
mort12=repmat(.25, 12, 1);
%Ks from the article 
K12=[ 0.39, 0.34, 0.3, 0.24, 0.23, 0.41, 0.2, 0.45, 0.14, 0.15, 0.38, 0.28;
0.22, 0.39, 0.34, 0.3, 0.27, 0.16, 0.15, 0.05, 0.38, 0.29, 0.37, 0.31;
 0.27, 0.22, 0.39, 0.34, 0.3, 0.07, 0.11, 0.05, 0.38, 0.41, 0.24, 0.25;
 0.3, 0.24, 0.22, 0.39, 0.34, 0.28, 0.12, 0.13, 0.27, 0.33, 0.04, 0.41;
 0.34, 0.3, 0.22, 0.2, 0.39, 0.4, 0.5, 0.26, 0.12, 0.29, 0.09, 0.16;];
%over 100 days
tspan=[0,100];
%Supply concentrations
S50=[6, 10, 14, 4, 9];
Content90=[0.04, 0.04, 0.07, 0.04, 0.04, 0.22, 0.1, 0.08 ,0.02, 0.17, 0.25, 0.03;
 0.08, 0.08, 0.08, 0.1, 0.08, 0.14, 0.22, 0.04, 0.18, 0.06, 0.2, 0.04;
 0.1, 0.1, 0.1, 0.1, 0.14, 0.22, 0.24, 0.12, 0.03, 0.24, 0.17, 0.01;
 0.05, 0.03, 0.03 ,0.03, 0.03, 0.09 ,0.07, 0.06, 0.03, 0.03, 0.11, 0.05;
 0.07 ,0.09, 0.07, 0.07, 0.07, 0.05, 0.24, 0.05, 0.08, 0.1, 0.02, 0.04;];
%over 100 days
tspan=[0,100];


growth_rate = @(mu, R, col) min(mu .* (R ./ (R + K12(:, col))), [], 'all');
sumfunction2 = @(row, y) Content90(row,1) * growth_rate(max_50(1), y(13:17), 1) * y(1) + ...
                  Content90(row,2) * growth_rate(max_50(2), y(13:17), 2) * y(2) + ...
                  Content90(row,3) * growth_rate(max_50(3), y(13:17), 3) * y(3)+ ...
                  Content90(row,4) * growth_rate(max_50(4), y(13:17), 4) * y(4) + ...
                  Content90(row,5) * growth_rate(max_50(5), y(13:17), 5) * y(5) + ...
                  Content90(row,6) * growth_rate(max_50(6), y(13:17), 6) * y(6)+ ...
                  Content90(row,7) * growth_rate(max_50(7), y(13:17), 7) * y(7) + ...
                  Content90(row,8) * growth_rate(max_50(8), y(13:17), 8) * y(8) + ...
                  Content90(row,9) * growth_rate(max_50(9), y(13:17), 9) * y(9) +...
                  Content90(row,10) * growth_rate(max_50(7), y(13:17), 10) * y(10) + ...
                  Content90(row,11) * growth_rate(max_50(8), y(13:17), 11) * y(11) + ...
                  Content90(row,12) * growth_rate(max_50(9), y(13:17), 12) * y(12);
f5 = @(t, y) [
    y(1) * (growth_rate(max_50(1), y(13:17), 1) - mort12(1)); 
    y(2) * (growth_rate(max_50(2), y(13:17), 2) - mort12(2)); 
    y(3) * (growth_rate(max_50(3), y(13:17), 3) - mort12(3));
    y(4) * (growth_rate(max_50(4), y(13:17), 4) - mort12(4)); 
    y(5) * (growth_rate(max_50(5), y(13:17), 5) - mort12(5)); 
    y(6) * (growth_rate(max_50(6), y(13:17), 6) - mort12(6)); 
    y(7) * (growth_rate(max_50(7), y(13:17), 7) - mort12(7)); 
    y(8) * (growth_rate(max_50(8), y(13:17), 8) - mort12(8));
    y(9) * (growth_rate(max_50(9), y(13:17), 9) - mort12(9)); 
    y(10) * (growth_rate(max_50(10), y(13:17), 10) - mort12(10));
    y(11) * (growth_rate(max_50(11), y(13:17), 11) - mort12(11)); 
    y(12) * (growth_rate(max_50(12), y(13:17), 12) - mort12(12)); 
   
    .25*(S5(1) - y(13)) - sumfunction2(1,y); 
    .25*(S5(2) - y(14)) - sumfunction2(2,y); 
    .25*(S5(3) - y(15)) - sumfunction2(3,y);
    .25*(S5(4) - y(16)) - sumfunction2(4,y); 
    .25*(S5(5) - y(17)) - sumfunction2(5,y);
];

figure;
hold on;
for i = 1:size(K12, 1)
    n = 12;
    R_0 = S50;
    N_0 = 0.1 + (1:12) / 100;
    y0 = [N_0, R_0];
    [t, y] = ode45(f5, tspan, y0);
    plot(t, y(:, 1:12));
end
xlabel('Time');
ylabel('Population');
title('Trajectories of the ODE System 12 Species, 5 resources');
grid on;
legend('Species 1', 'Species 2', 'Species 3', 'Species 4', 'Species 5','Species 6', 'Species 7', 'Species 8', 'Species 9', 'Species 10', 'Species 11', 'Species 12');
hold off;