% Homework 3: Disease Spreading, Simulation of Complex Systems FFR120
% Nicole Adamah 2022
%11.1
clear all; 
close all; 
clc;
%% Parameters:
N=1000;     
n=100;      
infect=10;  
d=0.8;      % Probability of random walk
g=0.01;     % Recovery rate
B=0.6;      % Infection rate
trials = 1;

for m = 1:trials
disp(m);
t = 0;
x=randi(n,1,N)-1;                       % Random location along x
y=randi(n,1,N)-1;                       % Random location along y
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
      y=mod(y+dy,n);                   % Performing walks, Periodic Boundary Conditions
      for i=1:N
        if (I(i)==true)&&(rand<B)       % Only infected agents can infect others,roll the dice for infection
          infection=(x==x(i))&(y==y(i));  % Determine the indices for the agents who sit at the site of infection
          S(infection)=false;             % There are no longer any susceptibles at infection area
          I(infection)=not(R(infection)); % the not recovered agents will turn infected at infection site
        end
      end
      recovery=(rand(1,N)<g);         % Recovery array
      R = R | (I&recovery);           % Recovery operation, the output is true if either or both of the inputs are true
      I = I & not(recovery);          % The recovered are no longer infected, the output is true when both inputs are true
      nrI(t+1)=sum(I); 
      nrR(t+1)=sum(R); 
      nrS(t+1)=sum(S);
      finalR = nrR(end);
      t=t+1;
    end
end
%% PLOT 11.1
h=figure; 
set(h,'Color','w','Units','Pixels');
a1=axes('Units','Pixels'); 
box on;  
hold on;
xlabel('$x$','FontSize',20,'Interpreter','Latex');
ylabel('$y$','FontSize',20,'Interpreter','Latex');
plot(x(I),y(I),'r.','MarkerSize',7);  
plot(x(R),y(R),'g.','MarkerSize',7);
plot(x(S),y(S),'b.','MarkerSize',7);
title(['t= ' num2str(t)]);             
xlim([-1 n]); 
ylim([-1 n]);

h1=figure; set(h,'Color','w','Units','Pixels')
set(h1,'Color','w','Units','Pixels');
a2=axes('Units','Pixels'); 
box on; 
hold on;
xlabel('Time Steps','FontSize',18,'Interpreter','Latex');
ylabel('Number of agents','FontSize',18,'Interpreter','Latex');
title(['\bf{$d=$' num2str(d) ', $\beta=$' num2str(B) ', $\gamma=$' num2str(g) '}'],'FontSize',12,'Interpreter','Latex')
plot(nrS,'b','LineWidth',2);
plot(nrI,'r','LineWidth',2); 
plot(nrR,'g','LineWidth',2); 
y1 = linspace(0,1000,6);
ylim([0 1000]);
legend('susceptible', 'infected','recovered');

