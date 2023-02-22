% Homework 3: Disease Spreading, Simulation of Complex Systems FFR120
% Nicole Adamah 2022
clear all; 
close all; 
clc;
%% Parameters:
N=1000;     
n=100;      
infect=10;  
d=0.8;      
g=0.01;     
B=0.8;      
mu_list=linspace(0.05,0.1,12);%11.2
mu=0.001;
g_list=[0.01,0.02];
trials = 1;
data1 = zeros(trials,length(mu_list));
data2 = zeros(trials,length(mu_list));
imm = 0.2;
for m = 1:trials
disp(m);
a = 1;

% for j = 1:length(mu_list)
% mu = mu_list(1,j);
t = 0;
x=randi(n,1,N)-1;                      
y=randi(n,1,N)-1;                      
[~,I]=sort((x-n/2).^2+(y-n/2).^2);      
preI=zeros(1,N);
preI(I(1:infect))=1;
I=logical(preI);                        
R=false(1,N);                           
S=logical(1-I);                         
IM=false(1,N); 
for h = 1:1000
%% SIMULATION
      dx=2*(round(rand(1,N))-0.5).*(rand(1,N)<d);     
      dy=2*(round(rand(1,N))-0.5).*(rand(1,N)<d);   
      x=mod(x+dx,n); 
      y=mod(y+dy,n);                    
      for i=1:N
        if (I(i)==true)&&(rand<B)       
          infection=(x==x(i))&(y==y(i));  
          S(infection)=false;             
          I(infection)=not(R(infection)); 
        end
      end
      recovery=(rand(1,N)<g);    
      R = R | (I&recovery);                                        
      I = I & not(recovery);          
      immun = (rand(1,N)<imm);
      IM = IM | (R&immun);
      S = S & (IM);
      nrI(t+1)=sum(I);               
      nrR(t+1)=sum(R); 
      nrS(t+1)=sum(S);
      nrIM(t+1)=sum(IM);
      finalR = nrR(end);
      finalIM = nrIM(end);
      t=t+1;
      end
%       data1(a,j) = data1(a,j) +  finalD;
    
if (m == 1)
    avgData1 = data1;
    avgData2 = data2;
else
    avgData1 = data1/trials;
    avgData2 = data2/trials;
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
title(['\bf{$d=$' num2str(d) ', $\beta=$' num2str(B) ', $\gamma=$' num2str(g) ', immun=' num2str(imm) '}'],'FontSize',12,'Interpreter','Latex')
plot(nrS,'b','LineWidth',2);
plot(nrI,'r','LineWidth',2); 
plot(nrR,'g','LineWidth',2); 
y1 = linspace(0,1000,6);
ylim([0 1000]);
legend('susceptible', 'infected','recovered');
h2=figure; 
set(h2,'Color','w','Units','Pixels')
a2=axes('Units','Pixels'); 
box on; 
hold on;
scatter(mu_list, avgData1(1,:), 'b')
xlabel('mu','FontSize',18,'Interpreter','Latex');
ylabel('Dinf','FontSize',18,'Interpreter','Latex');