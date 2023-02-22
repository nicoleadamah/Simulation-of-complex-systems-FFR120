%Prisoner's dilemma 13.2 & 13.3
%Nicole Adamah 2022
close all
clear all
clc
%% Prisoner's dilemma with multiple rounds. 

N = 10;
T = 0;%the player that defects while the other cooperates will be punshied with T years
S = 1.5;%The other player that cooperates while the other will be punished with S years 
R = 0.5;% both players cooperates
P = 1; % both players betray
m = 6;% other player's strategy
%T < R < P < S
n_range = 0:N;
years_in_prison = zeros(length(n_range),1);

for n=0:length(n_range)-1
    m = m;% other player's strategy, amount of ones(1=cooperate)
    if n < m
        years_in_prison(n+1) = (n)*R + (N-1-n)*P + (1)*T; %player 1 defects, player 2 cooperates

    elseif n > m
        years_in_prison(n+1) = (m)*R + (N-1-m)*P + (1)*S; %Player 2 defects, Player 1 cooperates.

    elseif n==m 
        years_in_prison(n+1) = (m)*R + (N-m)*P;
    end
end
%% 11.1a
f1 = figure;
scatter(n_range,years_in_prison,100,'filled');
x = [6 6]; y = [0 12];
line(x, y,'Color', 'black', 'LineStyle', '--')
ylim([6 10])
xlabel('n')
ylabel('Years in prison')
title(['\bf{S=$' num2str(S) ', R=$' num2str(R) ', P=$' num2str(P) ', T=$' num2str(T)  '}'],'FontSize',12,'Interpreter','Latex')
%% 11.1b
f2 = figure;
imagesc(0:N, 0:N, years_in_prison)
set(gca,'YDir','normal')
hold on;
colorbar
x = [-1 12]; y = [-2 11];
line(x, y, 'Color', 'black', 'LineStyle', '--')
xlabel('m')
ylabel('n')
title('Years in prison')