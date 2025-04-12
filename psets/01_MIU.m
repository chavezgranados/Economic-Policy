 % Topic 1: MIU model
 % Fuente: Karel Mertens
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 clear all;
 close all;
 
 % Utility function: u(C,m)=log(C)+phi*log (m)
 % Production function Y = AK^(1-alpha)* N^alpha
 
 % Parameter values
 alfa     = 0.58;
 beta     = 0.988;
 delta    = 0.025;
 mu       = 0.015;
 sigma    = 1;
 chi      = 1/0.39;
 yss      = 1;    % Output is normalized to 1
 rho_a    = 0.9;
 rho_theta   = 0.481; 
 sigma_a     = 0.01;
 sigma_theta = 0.009;
 
 % Steady state capital, investment and consumption:
 kss    = (1-alfa)*yss/(1/beta-1+delta);
 iss    = delta*kss; 
 css    = yss-iss;
 
 % Determine phi;
 phi      = (0.16)^chi*css^(-sigma)*(1+mu-beta)/(1+mu)*yss^chi;
 
 % Note: the implied values of A is
 A      = yss/kss^(1-alfa);

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 % Model Solution
 % The Linearized System is of the form
 % 
 % Fyp*Ey(t+1)+Fxp x(t+1)+Fy*y(t)+Fx*x(t)=0
 % 
 % where y are the controls and x are the states
K=1; mlag=2; a=3; theta=4; 
Y=1; C=2; I=3; m=4; R=5; pi=6;

% Note: R is 1+i where i is the nominal interest rate

Fy=zeros(10,6);  Fx=zeros(10,4); 
Fyp=zeros(10,6); Fxp=zeros(10,4);

%1. Investment Euler Equation
eqn = 1;
Fy(eqn,C)=sigma; 
Fyp(eqn,C)=-sigma;
Fxp(eqn,a)=(1-beta*(1-delta));
Fxp(eqn,K) = -alfa*(1-beta*(1-delta));

%2. Bond Euler Equation
eqn = 2;
Fy(eqn,C)=sigma; 
Fyp(eqn,C)=-sigma;
Fy(eqn,R)=1;
Fyp(eqn,pi) = -1;

%3. Money Demand Equation
eqn = 3;
Fy(eqn,C)=sigma; 
Fy(eqn,m)=-chi;
Fy(eqn,R)=-beta/(1+mu-beta);

%4. Resource constraint
eqn = 4;
Fy(eqn,C)  = css/yss;
Fx(eqn,K)  = - 1+alfa-iss/yss*(1-delta)/delta;
Fx(eqn,a)  = - 1;
Fxp(eqn,K) = iss/yss/delta;

%5. Production function
eqn = 5;
Fy(eqn,Y)  = -1;
Fx(eqn,K)  = (1-alfa);
Fx(eqn,a)  = 1;

%6. Capital Accumulation
eqn = 6;
Fy(eqn,I)  = -1;
Fx(eqn,K)  = - (1-delta)/delta;
Fxp(eqn,K) = 1/delta;

%7. Inflation
eqn = 7;
Fy(eqn,m)     = -1;
Fx(eqn,mlag)  =  1;
Fy(eqn,pi)    = -1;
Fx(eqn,theta) = 1;

%8. Money lag
eqn=8;
Fxp(eqn,mlag) = -1;
Fy(eqn,m) = 1;

%9. Technology process
eqn = 9;
Fxp(eqn,a)  = -1;
Fx(eqn,a)  = rho_a; 

%10. Money Growth process
eqn = 10;
Fxp(eqn,theta)  = -1;
Fx(eqn,theta)  = rho_theta;

At = [-Fxp -Fyp]; Bt = [Fx Fy];
[H,G]=solab(At,Bt,size(Fx,2));

%% Impulse responses
%% 1.Technology Shock
 s = [];
 s(:,1) = [0;0;0.01;0];
 for i=1:60
 s(:,i+1) = G*s(:,i);
 end
 
 z = (H*s)';
 s = s';
 % Plot Impulse Responses
 subplot(2,1,1)
 plot(100*s(:,K),'k-s','MarkerSize',3)
 hold on 
 plot(100*z(:,C),'r-+','MarkerSize',3)
 hold on 
 plot(100*z(:,Y),'b-d','MarkerSize',3)
 hold on 
 plot(100*z(:,I),'g-^','MarkerSize',3)
 hold on
 plot(100*s(:,a),'y-v','MarkerSize',3)
 hold on 
 ylabel('Percent Deviations')
 xlabel('Quarters')
 hold off
 legend('k','c','y','i','a')
 axis([0 25 -1 3.5])
 
 subplot(2,1,2)
 plot(100*(((1+mu)/beta*exp(z(:,R)))-((1+mu)/beta)),'k-s','MarkerSize',3)
 hold on 
 plot(100*((1/beta*exp(z(1:length(z)-1,R)-z(2:length(z),pi)))-(1/beta)),'r-+','MarkerSize',3)
 hold on 
 plot(100*(((1+mu)*exp(z(:,pi)))-(1+mu)),'b-d','MarkerSize',3)
 hold on 
 plot(100*z(:,m),'g-^','MarkerSize',3)
%  hold on 
%  plot(100*(z(:,Y)-z(:,m)),'m-^','MarkerSize',3)
 hold on
 plot(100*s(:,theta),'y-v','MarkerSize',3)
 hold on 
 ylabel('Percent Deviations')
 xlabel('Quarters')
 hold off
 legend('i','r','\pi','m','\theta')
 axis([0 25 -0.2 0.2])
 saveas(gcf,'miumodel1tech','psc2')
 
 
%% 2.Money growth Shock
 s=[];
 s(:,1) = [0;0;0;0.01];
 for i=1:60
 s(:,i+1) = G*s(:,i);
 end
 
 z = (H*s)';
 s = s';
 
 % Plot Impulse Responses
 subplot(2,1,1)
 plot(100*s(:,K),'k-s','MarkerSize',3)
 hold on 
 plot(100*z(:,C),'r-+','MarkerSize',3)
 hold on 
 plot(100*z(:,Y),'b-d','MarkerSize',3)
 hold on 
 plot(100*z(:,I),'g-^','MarkerSize',3)
 hold on
 plot(100*s(:,a),'y-v','MarkerSize',3)
 hold on 
 ylabel('Percent Deviations')
 xlabel('Quarters')
 hold off
 legend('k','c','y','i','a')
  axis([0 10 -0.005 0.015])
  
 subplot(2,1,2)
 plot(100*(((1+mu)/beta*exp(z(:,R)))-((1+mu)/beta)),'k-s','MarkerSize',3)
 hold on 
 plot(100*((1/beta*exp(z(1:length(z)-1,R)-z(2:length(z),pi)))-(1/beta)),'r-+','MarkerSize',3)
 hold on 
 plot(100*(((1+mu)*exp(z(:,pi)))-(1+mu)),'b-d','MarkerSize',3)
 hold on 
 plot(100*z(:,m),'g-^','MarkerSize',3)
 hold on 
 plot(100*(z(:,Y)-z(:,m)),'m-^','MarkerSize',3)
 hold on
 plot(100*s(:,theta),'y-v','MarkerSize',3)
 hold on 
 ylabel('Percent Deviations')
 xlabel('Quarters')
 hold off
 legend('i','r','\pi','m','v','\theta')
  axis([0 10 -1 2])
 saveas(gcf,'miumodel1money','psc2')