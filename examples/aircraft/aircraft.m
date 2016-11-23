%
% aircraft.m
%
% created on: 20.01.2016
%     author: rungger
%
% you need to run ./aircraft binary first 
%


function aircraft
clear set
close all

dbstop if error
%% simulation

% target set
lb=[63 -3*pi/180 0];
ub=[75 0 2.5];

cf=1;
eta=[cf*25.0/362 cf*3*pi/180/66 cf*56.0/334]; 
%% initial state

x0=[81 -pi/180 55];

% load controller from file
controller=Controller('reach.scs');

% simulate closed loop system
y=x0;
v=[];

while(1) 
    
    if ( lb(1) <= y(end,1)-eta(1)/2.0 && y(end,1)+eta(1)/2.0 <= ub(1) &&...
          lb(2) <= y(end,2)-eta(2)/2.0 && y(end,2)+eta(2)/2.0 <= ub(2) &&...
         lb(3) <= y(end,3)-eta(3)/2.0 && y(end,3)+eta(3)/2.0 <= ub(3) &&...
          ( -0.91 <= ((y(end,1)+eta(1)/2.0)*sin(y(end,2)-eta(2)/2.0))  ))
   %  if (   lb(1) <= y(end,1) && y(end,1)<= ub(1) &&...
      %      lb(2) <= y(end,2) && y(end,2) <= ub(2) &&...
       %     lb(3) <= y(end,3) && y(end,3) <= ub(3) &&...
        %    -0.91 <= y(end,1)*sin(y(end,2)) )
   
    break;
  end 
  u=controller.input(y(end,:));
  
 [t x]=ode45(@(t,x) aircraft_ode(t,x,u),[0 .25], y(end,:));
 y=[y; x(end,:)];
 v=[v; u(end,:)];
 
end


%% plot the domain
% colors
colors=get(groot,'DefaultAxesColorOrder');

figure
plot(y(:,1),'k.-','color',colors(2,:),'markersize',5);
ylabel('x_1');
hold on

figure
plot(y(:,2),'-o','color',colors(5,:),'markersize',5);
ylabel('x_2');
hold on
 
figure
plot(y(:,3),'.-','color',colors(1,:),'markersize',5);
ylabel('x_3');
hold on

%figure
%plot3(y(:,1),y(:,2),y(:,3),'.-','color',colors(5,:),'markersize',10);
%hold on 

figure
plot(v(:,1),'k.-','color',colors(2,:),'markersize',5);
ylabel('u_1');
hold on

figure
plot(v(:,2),'k.-','color',colors(2,:),'markersize',5);
ylabel('u_2');
hold on



function dxdt = aircraft_ode(t,x,u)

  dxdt = zeros(3,1);
  c=(1.25+4.2*u(2));

  dxdt(1)= 1/60e3 *(u(1)*cos(u(2))-(2.7+3.08*c^2)*x(1)^2-60e3*9.81*sin(x(2)));
  dxdt(2)= (1/(60e3*x(1)))*(u(1)*sin(u(2))+68.6*c*x(1)^2-60e3*9.81*cos(x(2)));
  dxdt(3)=x(1)*sin(x(2));

end

 end
