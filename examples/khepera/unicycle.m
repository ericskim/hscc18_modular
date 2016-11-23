%
% unicycle.m
%
% created on: 21.01.2016
%     author: rungger
%
% see readme file for more information on the unicycle example
%
% you need to run ./unicycle binary first 
%
% so that the files: unicycle_ss.bdd 
%                    unicycle_obst.bdd
%                    unicycle_target.bdd
%                    unicycle_controller.bdd 
% are created
%

function unicycle
clear set
close all



%% simulation

%% target set
L=[3 0 0; 0 3 0; 0 0 .1];
c=[5.25;5.25; 0];

% initial state
x0=[12.490160 17.810000 2.974979];

% load controller from file
con=Controller('reach.scs');

%keyboard



y=x0;
v=[];
while(1)

  
  if ( (y(end,:)-c')*L'*L*(y(end,:)'-c)<=1 )
    break;
  end 

  u=con.input(y(end,:));
  v=[v; u(1,:)];
  [t x]=ode45(@unicycle_ode,[0 .3], y(end,:),[],u(1,:));
  x(end,:)
  y=[y; x(end,:)];
end

%keyboard


% load the symbolic set containig the abstract state space
%set=UniformGrid('problemdomain.scs');
%plotCells(set,'projection',[1 2],'facecolor','none','edgec',[0.8 0.8 0.8],'linew',.1)
set=UniformGrid('problemdomain.scs','projection',[1,2]);
plotCells(set,'facecolor','none','edgec',[0.8 0.8 0.8],'linew',.1);
hold on


%% plot the unicycle domain
% colors
colors=get(groot,'DefaultAxesColorOrder');

%% load the symbolic set containig obstacles
set=UniformGrid('obstacles.scs','projection',[1 2]);
plotCells(set,'facecolor',colors(1,:)*0.5+0.5,'edgec',colors(4,:),'linew',.1)

% plot the real obstacles and target set
plot_domain
hold on

%% load the symbolic set containig target set
set=UniformGrid('target.scs','projection',[1 2]);
plotCells(set,'facecolor',colors(2,:)*0.5+0.5,'edgec',colors(2,:),'linew',.1)

% plot initial state  and trajectory
plot(y(:,1),y(:,2),'k.-')
plot(y(1,1),y(1,2),'.','color',colors(5,:),'markersize',20)


box on
axis([0 24.5 0 23.6])


end

function dxdt = unicycle_ode(t,x,u)

  dxdt = zeros(3,1);

  dxdt(1)=u(1)*cos(x(3));
  dxdt(2)=u(1)*sin(x(3));
  dxdt(3)=u(2);


end

function plot_domain

colors=get(groot,'DefaultAxesColorOrder');


% v=[0.2   12.6200;0.7938   12.6200;0.2   17.2200;0.7938   17.2200];
% patch('vertices',v,'faces',[1 2 4 3],'facec',colors(1,:),'edgec',colors(1,:));
% v=[2.4587   12.6200;3.8139   12.6200;2.4587   17.2200;3.8139   17.2200];
% patch('vertices',v,'faces',[1 2 4 3],'facec',colors(1,:),'edgec',colors(1,:));
% v=[0.2   17.2000;3.8139   17.2000;0.2   23.5000;3.8139   23.5000];
% patch('vertices',v,'faces',[1 2 4 3],'facec',colors(1,:),'edgec',colors(1,:));
% v=[7.3568    7.0800;10.9190    7.0800;7.3568   10.6400;10.9190   10.6400];
% patch('vertices',v,'faces',[1 2 4 3],'facec',colors(1,:),'edgec',colors(1,:));
% v=[7.3374   10.6400;8.3442   10.6400;7.3374   14.2600;8.3442   14.2600];
% patch('vertices',v,'faces',[1 2 4 3],'facec',colors(1,:),'edgec',colors(1,:));
% v=[7.3181   13.9600;10.8997   13.9600;7.3181   14.8000;10.8997   14.8000 ];
% patch('vertices',v,'faces',[1 2 4 3],'facec',colors(1,:),'edgec',colors(1,:));
% v=[10.2221   19.8600;10.8029   19.8600;10.2221   23.5000;10.8029   23.5000];
% patch('vertices',v,'faces',[1 2 4 3],'facec',colors(1,:),'edgec',colors(1,:));
% v=[7.3762   18.5400;10.8610   18.5400;7.3762   20.1000;10.8610   20.1000];
% patch('vertices',v,'faces',[1 2 4 3],'facec',colors(1,:),'edgec',colors(1,:));
% v=[14.3651   11.6000;17.6757   11.6000;14.3651   14.7800;17.6757   14.7800];
% patch('vertices',v,'faces',[1 2 4 3],'facec',colors(1,:),'edgec',colors(1,:));
% v=[14.3845    7.0600;17.7144    7.0600;14.3845    8.1600;17.7144    8.1600];
% patch('vertices',v,'faces',[1 2 4 3],'facec',colors(1,:),'edgec',colors(1,:));
% v=[14.1134         0;24.0451         0;14.1134    3.5000;24.0451    3.5000];
% patch('vertices',v,'faces',[1 2 4 3],'facec',colors(1,:),'edgec',colors(1,:));
% v=[21.1798   12.4600;23.9870   12.4600;21.1798   19.4400;23.9870   19.4400];
% patch('vertices',v,'faces',[1 2 4 3],'facec',colors(1,:),'edgec',colors(1,:));
% v=[ 7.3181   0.0200;10.8416   0.0200;7.3181    3.3600;10.8416    3.3600];
% patch('vertices',v,'faces',[1 2 4 3],'facec',colors(1,:),'edgec',colors(1,:));
% v=[17.0562   21.6400;17.6950   21.6400;17.0562   23.5000;17.6950   23.5000];
% patch('vertices',v,'faces',[1 2 4 3],'facec',colors(1,:),'edgec',colors(1,:));
% v=[21.1411    6.9200;23.9677    6.9200;21.1411    9.0000;23.9677    9.0000];
% patch('vertices',v,'faces',[1 2 4 3],'facec',colors(1,:),'edgec',colors(1,:));
% v=[14.3070   18.3600;17.6563   18.3600;14.3070   21.9000;17.6563   21.9000];
% patch('vertices',v,'faces',[1 2 4 3],'facec',colors(1,:),'edgec',colors(1,:));

v=[0.2517    4.6800;0.7938    4.6800;0.2517    9.1000;0.7938    9.1000];
patch('vertices',v,'faces',[1 2 4 3],'facec',colors(1,:),'edgec',colors(1,:));
v=[2.4394    4.7000;3.8914    4.7000;2.4394    9.0600;3.8914    9.0600];
patch('vertices',v,'faces',[1 2 4 3],'facec',colors(1,:),'edgec',colors(1,:));
v=[ 0.2517    2.5800;3.9301    2.5800;0.2517    5.0000;3.9301    5.0000];
patch('vertices',v,'faces',[1 2 4 3],'facec',colors(1,:),'edgec',colors(1,:));


v=[5.3434    0.3000;5.7499    0.3000;5.3434    3.0200;5.7499    3.0200];
patch('vertices',v,'faces',[1 2 4 3],'facec',colors(1,:),'edgec',colors(1,:));
v=[5.3434    7.0400;5.7499    7.0400;5.3434    8.9800;5.7499   8.9800];
patch('vertices',v,'faces',[1 2 4 3],'facec',colors(1,:),'edgec',colors(1,:));
v=[5.3434    11.7800;5.7499    11.7800;5.3434    14.6400;5.7499   14.6400];
patch('vertices',v,'faces',[1 2 4 3],'facec',colors(1,:),'edgec',colors(1,:));
v=[5.3434    18.6400;5.7499    18.6400;5.3434    21.6400;5.7499   21.6400];
patch('vertices',v,'faces',[1 2 4 3],'facec',colors(1,:),'edgec',colors(1,:));
v=[2.2070   10.7000;3.9301   10.7000;2.2070   11.0600;3.9301   11.0600];
patch('vertices',v,'faces',[1 2 4 3],'facec',colors(1,:),'edgec',colors(1,:));


v=[12.2323    7.7200;12.6195    7.7200;12.2323   14.2800;12.6195   14.2800];
patch('vertices',v,'faces',[1 2 4 3],'facec',colors(1,:),'edgec',colors(1,:));
v=[12.2323    18.7200;12.6195    18.7200;12.2323   21.8400;12.6195   21.8400];
patch('vertices',v,'faces',[1 2 4 3],'facec',colors(1,:),'edgec',colors(1,:));
v=[12.2323    0.6000;12.6195    0.6000;12.2323   2.9200;12.6195   2.9200];
patch('vertices',v,'faces',[1 2 4 3],'facec',colors(1,:),'edgec',colors(1,:));


v=[ 8.0150    4.9400;10.6093    4.9400;8.0150    5.3000;10.6093    5.3000];
patch('vertices',v,'faces',[1 2 4 3],'facec',colors(1,:),'edgec',colors(1,:));
v=[ 14.4426    4.9400;17.1336    4.9400;14.4426    5.3000;17.1336    5.3000];
patch('vertices',v,'faces',[1 2 4 3],'facec',colors(1,:),'edgec',colors(1,:));


v=[7.6085   16.4600;10.5706   16.4600;7.6085   16.8000;10.5706   16.8000];
patch('vertices',v,'faces',[1 2 4 3],'facec',colors(1,:),'edgec',colors(1,:));
v=[14.4232   16.4600;17.6757   16.4600;14.4232   16.8000;17.6757   16.8000];
patch('vertices',v,'faces',[1 2 4 3],'facec',colors(1,:),'edgec',colors(1,:));


v=[19.1245    7.4600;19.4536    7.4600;19.1245    9.7200;19.4536    9.7200];
patch('vertices',v,'faces',[1 2 4 3],'facec',colors(1,:),'edgec',colors(1,:));
v=[19.1245    12.6600;19.4536    12.6600;19.1245    15.0600;19.4536    15.0600];
patch('vertices',v,'faces',[1 2 4 3],'facec',colors(1,:),'edgec',colors(1,:));
v=[19.1245    18.5600;19.4536    18.5600;19.1245    21.4400;19.4536    21.4400];
patch('vertices',v,'faces',[1 2 4 3],'facec',colors(1,:),'edgec',colors(1,:));

v=[21.1992   10.5800;22.7674   10.5800;21.1992   10.8600;22.7674   10.8600];
patch('vertices',v,'faces',[1 2 4 3],'facec',colors(1,:),'edgec',colors(1,:));


end