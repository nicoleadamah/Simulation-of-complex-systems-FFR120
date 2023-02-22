%Prisoner's dilemma 13.2 & 13.3
%Nicole Adamah 2022
close all 
clear all
clc 
%% PARAMETERS
R = 0.07; 
S = 0.08; 
P = 1; 
T=0;
timesteps = 500;
L=30;
N = 7;
mu = 0.01;
d=7;%number of defectors, nr 5 is for 11.2c, nr 6 is for 11.3 and nr 7 for 11.4
lattice=Initialize(L, N, d);
initialLattice = lattice;
data = zeros(timesteps, N+1);%for 11.4
%% SIMULATION
for t = 1:timesteps
    lattice_years = zeros(L);
    for i = 1:L
        for j = 1:L
%       Interaction with the four nearest von Neumann neighbors (top, bottom, left, right). Wrap around if at edges        
            if (i == 1)             
                topNeighbor=lattice(L,j);
            else
               topNeighbor = lattice(i-1,j); 
            end
            if (i == L)             
                bottomNeighbor=lattice(1,j);
            else
                bottomNeighbor = lattice(i+1,j);
            end 
            if (j == 1)             
               leftNeighbor=lattice(i,L);  
            else
               leftNeighbor = lattice(i,j-1);  
            end
            if (j == L)             
               rightNeighbor=lattice(i,1);
            else
               rightNeighbor = lattice(i,j+1);
            end
            currentAgent=lattice(i,j);
            myNeighbors=[topNeighbor bottomNeighbor leftNeighbor rightNeighbor]; 
            years=[];
                for k=1:4
                    years(k)=PrisonersModel1(currentAgent,T,S,R,P,myNeighbors(k),N);

                end  
                
            lattice_years(i, j) = sum(years);      
        end
    end
    
    % Update strategy
    newLattice = lattice;
    for i = 1:L
        for j = 1:L
            if (i == 1)             
                topNeighbor1=lattice_years(L,j);
            else
               topNeighbor1 = lattice_years(i-1,j); 
            end
            if (i == L)             
                bottomNeighbor1=lattice_years(1,j);
            else
                bottomNeighbor1 = lattice_years(i+1,j);
            end 
            if (j == 1)             
               leftNeighbor1=lattice_years(i,L);  
            else
               leftNeighbor1 = lattice_years(i,j-1);  
            end
            if (j == L)             
               rightNeighbor1=lattice_years(i,1);
            else
               rightNeighbor1 = lattice_years(i,j+1);
            end
            currentAgent1 = lattice_years(i,j);
            all=[topNeighbor1, bottomNeighbor1, leftNeighbor1, rightNeighbor1,currentAgent1 ]; 

            optima=min(all);%change to max for 11.2c
            idx=find(all==optima);               
            if length(idx)>1   % If two scores tie randomly select one of them  
                idx1=idx(randi([1 length(idx)]));
            else 
                idx1=idx;
            end
            % Update original lattice with best strategies    
            if (idx1 == 1)
                if (i-1== 0)
                    newLattice(i, j) = lattice(L, j);
                else
                    newLattice(i, j) = lattice(i - 1, j);
                end

            elseif (idx1 == 2)
                if (i+1== L+1)
                    newLattice(i, j) = lattice(1, j);
                else
                    newLattice(i, j) = lattice(i + 1, j);
                end 
                
            elseif (idx1 == 3)
                     
                if (j-1== 0)
                    newLattice(i, j) = lattice(i, L);
                else
                    newLattice(i, j) = lattice(i, j - 1);
                end  
                
            elseif (idx1 == 4)
                if (j+1 == L+1)
                    newLattice(i, j) = lattice(i, 1);
                else
                    newLattice(i, j) = lattice(i, j + 1);
                end
            elseif(idx1 == 5)
                    newLattice(i, j) = lattice(i, j);
            end
        end
    end
    lattice = newLattice;
%   Mutate with probability mu 
    lattice = mutate(lattice, 2, L, N, mu); % change between 2 and 3 depding on problem 11.2 or 11.3 
    idxData = zeros(1, N+1);
    for i = 1:N+1
        idxData(i) = sum(lattice(:) == i-1);
    end
    data(t, :) = idxData;
end 
%% PLOT 11.2 
% figure; 
% imagesc(initialLattice); 
% colormap([1 0 0;0 0 1]);  
% title(['\bf{$Initial defectors=$' num2str(d) ', $R=$' num2str(R)  '}'],'FontSize',12,'Interpreter','Latex')
% ylabel('t=0'); 
% 
% figure; 
% imagesc(lattice); 
% colormap([1 0 0;0 0 1]);  
% title(['\bf{$Initial defectors=$' num2str(d) ', $R=$' num2str(R)  '}'],'FontSize',12,'Interpreter','Latex')
% ylabel('t=20'); 
%% PLOT 11.3 
% figure; 
% imagesc(lattice); 
% colormap([1 0 0;0 0 1]);  
% title(['\bf{$S=$' num2str(S) ', $R=$' num2str(R)  '}'],'FontSize',12,'Interpreter','Latex')
% ylabel('t=100');
%% PLOT 11.4 
% figure; 
% imagesc(initialLattice); 
% colormap(flipud(jet(N+1)))
% colorbar
% title(['\bf{$S=$' num2str(S) ', $R=$' num2str(R)  '}'],'FontSize',12,'Interpreter','Latex')
% ylabel('t=0'); 
% 
% figure; 
% imagesc(lattice); 
% colormap(flipud(jet(10)))
% colorbar
% title(['\bf{$S=$' num2str(S) ', $R=$' num2str(R)  '}'],'FontSize',12,'Interpreter','Latex')
% ylabel('t=200'); 
% figure
% for i = 1:N+1
%     plot(1:1:timesteps, data(:,i), 'color', colors(i,:))
%     hold on
% end
% legend('n = 0', 'n = 1', 'n = 2', 'n = 3', 'n = 4','n = 5', 'n = 6', 'n = 7')
% xlabel('Time')
% ylabel('Population fraction')
% title('Number of agents with a particular strategy,  n = [0, N]')
%% PLOT 11.5
data = data(101:500, :);
v=var(data);
figure; 
imagesc(v);
set(gca, 'YDir', 'normal')
colormap(flipud(jet(10)))
colorbar
title(['\bf{$S=$' num2str(S) ', $R=$' num2str(R)  '}'],'FontSize',12,'Interpreter','Latex')
% ylabel('t=500'); 
%% FUNCTIONS
%Prison dilemma
function years_in_prison1 = PrisonersModel1(n, T,S,R,P,m,N) 
years_in_prison1 = 0;
for i=1:N
    if n < m
        years_in_prison1 = (n)*R + (N-1-n)*P + (1)*T; %player 1 defects, player 2 cooperates

    elseif n > m
        years_in_prison1 = (m)*R + (N-1-m)*P + (1)*S; %Player 2 defects, Player 1 cooperates.

    elseif n==m 
        years_in_prison1 = (m)*R + (N-m)*P;
    end
end
end
%% Mutation 
function latticeM = mutate(lattice, mut, L, N, mu)
if (mut==2)
     for i = 1:L
        for j = 1:L
            if (rand() < mu)
                lattice(i, j) = randi(N+1) - 1;
            end
        end
     end
     latticeM=lattice;
end
if (mut==3)
    for i = 1:L
        for j = 1:L
            if (rand() < mu)
                if (rand() < 0.5)
                    lattice(i, j) = 0;
                else
                    lattice(i, j) = N;
                end
            end
        end
     latticeM=lattice;        
    end 
end
end 
    
    
%% Initialize the lattice with defectors
function InitialLattice = Initialize(L,N, defector) 

% Initialization. 
if (defector == 1)
    lattice = N*ones(L);
    lattice(ceil(L/2), ceil(L/2)) = 0;
    InitialLattice = lattice;
end
if (defector == 2)
    lattice = N*ones(L);
    lattice(ceil(L/3), ceil(2*L/3)) = 0;
    lattice(ceil(2*L/3), ceil(L/3)) = 0;
    InitialLattice = lattice;
end
if (defector == 3)
    lattice = N*ones(L);
    lattice(ceil(L/2), ceil(L/2)) = 0;
    lattice(ceil(L/4), ceil(3*L/4)) = 0;
    lattice(ceil(3*L/4), ceil(L/4)) = 0;
    InitialLattice = lattice;
end
if (defector == 4)
    lattice = N*ones(L);
    lattice(ceil(L/5), ceil(4*L/5)) = 0;
    lattice(ceil(2*L/5), ceil(3*L/5)) = 0;
    lattice(ceil(3*L/5), ceil(2*L/5)) = 0;
    lattice(ceil(4*L/5), ceil(L/5)) = 0;
    InitialLattice = lattice;
end
if (defector==5)
    lattice = zeros(L);
    lattice(ceil(L/2), ceil(L/2)) = N;
    InitialLattice = lattice;
end
if (defector == 6)
    strategy = [0, N];
    r = randi([1, 2], L); 
    lattice = strategy(r);
    InitialLattice = lattice;
end 
if (defector == 7)
    lattice = randi([0, N], L); 
    InitialLattice = lattice;
end 

end
