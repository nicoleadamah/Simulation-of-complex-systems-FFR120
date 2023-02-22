% Homework 3: Disease Spreading, Simulation of Complex Systems FFR120
% Nicole Adamah 2022
% 11.2 
clear all; 
close all; 
clc;
%% Parameters:
N=1000;     
n=100;      
infect=10;  
d=0.8;      
b_list=linspace(0.05,1,12);
g_list=[0.01,0.02];
trials = 4;
data1 = zeros(trials,length(b_list));
data2 = zeros(trials,length(b_list));

for m = 1:trials
disp(m);
a = 1;
for k =1:length(g_list)
    g = g_list(1,k);
for j = 1:length(b_list)
    B = b_list(1,j);
    t = 0;
    x=randi(n,1,N)-1;                       % Random Location along x
    y=randi(n,1,N)-1;                       % Random Location along y
    [~,I]=sort((x-n/2).^2+(y-n/2).^2);      % Find closest ones to the center
    preI=zeros(1,N);
    preI(I(1:infect))=1;
    I=logical(preI);                        % Infection status array
    R=false(1,N);                           % Recovered status array
    S=logical(1-I);                         % Susceptible status array
    for h = 1:1000
%% SIMULATION
      dx=2*(round(rand(1,N))-0.5).*(rand(1,N)<d);     % Random Walks along x
      dy=2*(round(rand(1,N))-0.5).*(rand(1,N)<d);     % Random Walks along x
      x=mod(x+dx,n); 
      y=mod(y+dy,n);                   % Performing walks, Periodic Boundary Conditions, Alternative -solid walls- would be x(x>n)=n;
      for i=1:N
        if (I(i)==true)&&(rand<B)       % Only infected agents can infect others, and we also roll the dice for infection
          infection=(x==x(i))&(y==y(i));  % Determine the indices for those who sit at the site of infection
          S(infection)=false;             % There are no longer any susceptibles at infection area
          I(infection)=not(R(infection)); % All non-recovered agents will turn infected at infection site
        end
      end
      recovery=(rand(1,N)<g);         % Recovery array
      R = R | (I&recovery);           % Recovery operation
      I = I & not(recovery);          % The ones recovered are no longer infected
      nrI(t+1)=sum(I); 
      nrR(t+1)=sum(R); 
      nrS(t+1)=sum(S);
      finalR = nrR(end);
      t=t+1;
    end
    if k == 1
        data1(a,j) = data1(a,j) +  finalR;
    else
        data2(a,j) = data2(a,j) +  finalR;
    end
end
end 
if (m == 1)
    avgData1 = data1;
    avgData2 = data2;
else
    avgData1 = data1/trials;
    avgData2 = data2/trials;
end
end
%% PLOT 11.2 a & b
h1=figure; set(h1,'Color','w','Units','Pixels')
set(h1,'Color','w','Units','Pixels');
a2=axes('Units','Pixels'); 
box on; 
hold on;
scatter(b_list, avgData1(1,:), 'b')
scatter(b_list, avgData2(1,:),'g')
xlabel('Beta','FontSize',18,'Interpreter','Latex');
ylabel('Rinf','FontSize',18,'Interpreter','Latex');
legend('g = 0.01', 'g=0.02');



%% PLOT 11.2 c

