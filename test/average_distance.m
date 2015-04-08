close all
clear all

% We have an area of membrane of 
% A = pi*R^2 = 0.03 units^2
% So
R = sqrt(0.03/pi)

N_sims = 1000000;
N_channels = 5;

% Sampling in a circle
figure
theta=linspace(0, 2*pi, 100000)'; % Just for drawing a circle
radius = R + zeros(length(theta),1);
plot(radius.*cos(theta), radius.*sin(theta), 'k-')
hold on

mean_distance_all = zeros(N_sims, 1);
mean_distance_nearest = zeros(N_sims, 1);
pairings = combnk(1:N_channels,2); % This lists channel pairings

for simulations = 1:N_sims

    for i=1:N_channels
        % Sampling uniformly in a circle is slightly subtle!
        t(i) = 2*pi*rand;
        r(i) = sqrt(rand).*R;        
    end
    
    x = r.*cos(t);
    y = r.*sin(t);
    
    plot(x, y,'b.')
        
    % For the average distance between each two
    for pair = 1:length(pairings)
        x_dist = x(pairings(pair,1))-x(pairings(pair,2));
        y_dist = y(pairings(pair,1))-y(pairings(pair,2));
        dist_each_pair(pair) = sqrt(x_dist.^2 + y_dist.^2);
    end    
    mean_distance_all(simulations) = mean(dist_each_pair);    
    
    % For the average distance of nearest neighbour
    % Loop over each channel
    for i=1:N_channels
        distance_to_others = [];
        for j=1:N_channels
           % Don't get distance to ourselves
           if (i==j)
               continue
           end
           distance_to_others = [distance_to_others; sqrt((x(i)-x(j)).^2 + (y(i)-y(j)).^2)];
        end        
        dist_nearest_neighbour(i) = min(distance_to_others);
    end
    assert(length(dist_nearest_neighbour)==N_channels)
    mean_distance_nearest(simulations) = mean(dist_nearest_neighbour);
end

figure
subplot(2,1,1)
hist(mean_distance_all)
subplot(2,1,2)
hist(mean_distance_nearest)

% We are currently in um. Need to convert to nm.
fprintf('Mean distance between any two channels picked at random is %6.4g nm.\n',mean(mean_distance_all)*1000)
fprintf('Mean distance between any channel picked at random and nearest neighbour is %6.4g nm.\n',mean(mean_distance_nearest)*1000)




