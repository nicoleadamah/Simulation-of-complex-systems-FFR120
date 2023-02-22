%% Standard Vicsek model.
close all
clear all
clc 
tic
% PARAMETERS
t_steps = 10000;
N = 100; 
L = 100; 
R = 1; 
dt = 1; 
v = 1; 
eta = 0.1;
limit = 45; % sight limit in degrees
k = 8;
%% EXERCISE 8.4 - 8.7 & 8.10 - 8.11
for m = 1:1
    [x, theta] = Initialize(N, L);
    gA = zeros(1, t_steps);
    gC = zeros(1, t_steps);
    for i = 2:t_steps
    neighbours = zeros(N) ~= 0;
    for j = 1:N
        neighbours(j, :) = Find(x, j, R, L);
    end
    velocities = UpdateV(v, theta, N);
    x = UpdatePos(x, velocities, dt, L);
%     theta = UpdateOrient(theta, neighbours, eta, dt,N);
    theta = UpdateOrientKNN(theta, x, k, eta, dt, limit, N); %unmute to use KNN, dont forget to mute row above
    gA(i) = GlobalAlignment(N, velocities, v);
    gC(i) = GlobalClustering(x, L, R);
    if (rem(i, 1000) == 0)
        disp(i)
    end  
    end
end
% 
%% PLOT 
% 8.4 a
voronoi(x(:,1), x(:,2))
axis equal
% 8.4 b
f2 = figure;
plot(gA)
hold on
plot(gC)
legend('?', 'C')
ylim([0,1.1])
xlabel('t')
ylabel('? & C')
toc
%% FUNCTIONS
%Initialize positions and theta
function [pos,the] = Initialize(N, L)
     pos = zeros(N, 2);
     the = zeros(N, 1);
     for i = 1:N
        r1 = rand(); 
        r2 = rand();
        x = L*r1 - L/2;
        y = L*r2 - L/2;
        pos(i, :) = [x, y];
        the(i) = rand()*2*pi;
     end 
end

% Update positions with periodic boundary conditions
function pos = UpdatePos(pos, velocities,dt, L)
    for i = 1:size(pos, 1)
        dR = velocities.*dt;
        pos(i, :) = pos(i, :) + dR(i, :);
    end
    for i = 1:size(pos, 1)
        for j = 1:size(pos, 2)
            if (pos(i, j) > L/2)
                pos(i, j) = pos(i, j) - L;
            elseif (pos(i ,j) < -L/2)
                pos(i, j) = pos(i, j) + L;
            end
        end
    end
end

% Update velocities
function velocities = UpdateV(v, theta,N)
    velocities = zeros(N, 2);
    for i = 1:N
        velocities(i, 1) = v*(cos(theta(i)));
        velocities(i, 2) = v*(sin(theta(i)));
    end
end

% Find neighbors within Rf to update theta
function Neighbours = Find(p, index, Rf, L)
    Neighbours = zeros(1, size(p, 1)) ~= 0;
    current_p = p(index, :);
    % For cases when the absolute(dist) > L/2 = wrapped around 
    for i = 1:size(p, 1)
        xdist = current_p(1) - p(i, 1);
        if (abs(xdist) > L/2)
            xdist = L - xdist;
        end
        ydist = current_p(2) - p(i, 2);
        if (abs(ydist) > L/2)
            ydist = L - ydist;
        end
        dist = sqrt(xdist^2 + ydist^2); %Euclidean
        if (dist < Rf)
            Neighbours(i) = true;
        end
    end
end
% Update theta
function theta = UpdateOrient(old_theta, neighbour, eta, dt,N)
    theta = zeros(N, 1);
    for i = 1:N
        r = rand - 0.5;
        S = find(neighbour(i,:));
        if (length(S) == 1)
            mean_theta = old_theta(i);
        else
            mean_theta = atan2(mean(sin(old_theta(S))), mean(cos(old_theta(S))));
        end
        theta(i) = mean_theta + eta*r*dt;
    end
end
% Find neighbors with KNN
function theta = UpdateOrientKNN(old_theta, x, k, eta, dt, limit, N)
    theta = zeros(N, 1);
    for i = 1:N
        dist = zeros(N, 1);
        for m = 1:N
            dist(m) = Inf;
        end
        for j = 1:N
            angle = acosd(dot(x(i, :), x(j, :))/(norm(x(i, :))*norm(x(j, :))));%cosine similarity
            if (angle < limit)
                dist(j) = norm(x(i, :) - x(j, :));
            end
        end
        [~, idx] = sort(dist);
        nearest = idx(1:k);
        mean_theta = atan2(mean(sin(old_theta(nearest))), mean(cos(old_theta(nearest))));
        theta(i) = mean_theta + eta*(rand - 0.5)*dt;
    end
end
% Calculate global alignment coefficient
function gA = GlobalAlignment(N, velocities, v)
    gA = norm(sum(velocities)) /(v* N);
end


%Calculate global clustering coefficient
function gA = GlobalClustering(x, L, Rf)
pos = [[x(:, 1) - L, x(:, 2) + L]; [x(:, 1), x(:, 2) + L]; [x(:, 1) + L, x(:, 2) + L];
     [x(:, 1) - L, x(:, 2)];     [x(:, 1), x(:, 2)];     [x(:, 1) + L, x(:, 2)];
     [x(:, 1) - L, x(:, 2) - L]; [x(:, 1), x(:, 2) - L]; [x(:, 1) + L, x(:, 2) - L]]; 
[vertices, cells] = voronoin(pos);
area_idx = zeros(1, size(pos, 1));
area = zeros(sum(area_idx), 1);

for i = 1:size(pos, 1)
    if (max(abs(pos(i, 1))) < 1/2*L)
        if (max(abs(pos(i, 2))) < 1/2*L)
            area_idx(i) = 1;
        end
    end
end

j = 1;
for i = find(area_idx)
    V1 = vertices(cells{i}, 1); 
    V2 = vertices(cells{i}, 2);
    area(j) = polyarea(V1, V2);
    j = j + 1;
end
gA = sum(area < pi*Rf^2)/length(area);
end