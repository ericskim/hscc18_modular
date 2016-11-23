%
% dcdc.m
%
% created on: 09.10.2015
%     author: rungger
%
% see readme file for more information on the safety example
%
% you need to run ./dcdc binary first 
%
% so that the files: safe.scs
%
% is created
%

function dcdc
clear set
close all

%% simulation

% initial state
x0=[1.35 5.755];
tau_s=0.5;

% load controller from file
controller=Controller('safe.scs');
domainsize=size(controller.domain);
% simulate closed loop system
y=x0;
v=[];
loop=100;
while(loop>0)
	loop=loop-1;

  u=controller.input(y(end,:));
  
  %-------------here choose your controller input-------------%
  in=u(1,:);
  %-----------------------------------------------------------%
  
  v=[v; in];
  [t x]=ode45(@unicycle_ode,[0 tau_s], y(end,:), odeset('abstol',1e-4,'reltol',1e-4),in');

  y=[y; x(end,:)];
end


%% plot the vehicle domain
% colors
colors=get(groot,'DefaultAxesColorOrder');


% load the symbolic set containig the abstract state space
set=UniformGrid('safedomain.scs');
plot(set,'.')
hold on

% plot initial state  and trajectory
plot(y(:,1),y(:,2),'k.-')
hold on
plot(y(1,1),y(1,2),'.','color',colors(5,:),'markersize',20)

% plot safe set
v=[1.15 5.45; 1.55 5.45; 1.15 5.85; 1.55 5.85 ];
patch('vertices',v,'faces',[1 2 4 3],'facecolor','none','edgec',colors(2,:),'linew',1)
hold on


box on
axis([1.1 1.6 5.4 5.9])


%set(gcf,'paperunits','centimeters','paperposition',[0 0 16 10],'papersize',[16 10])


end

function dxdt = unicycle_ode(t,x,u)
    % parameter initialization
    xc=70;
    xl=3;
    rc=0.005;
    rl=0.05;
    ro=1;
    vs=1;

    dxdt = zeros(2,1);
    switch u
        case 1
            dxdt(1)=-rl/xl*x(1)+vs/xl;
			dxdt(2)=-1/(xc*(ro+rc))*x(2);
        case 2
            dxdt(1)=-(1/xl)*(rl+ro*rc/(ro+rc))*x(1)-(1/xl)*ro/(5*(ro+rc))*x(2)+vs/xl;
            dxdt(2)=(1/xc)*5*ro/(ro+rc)*x(1)-(1/xc)*(1/(ro+rc))*x(2);
    end

end

