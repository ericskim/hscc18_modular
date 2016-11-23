function testProjection
% aims to test if the domain projection is correct

% test Controller.domain
controller=Controller('safe.scs');
controller1=Controller('safe.scs','projection',[1]);
controller2=Controller('safe.scs','projection',[2]);

do=controller.domain;
do1=controller1.domain;
do2=controller2.domain;

d1cmp=sort(unique(do(:,1),'rows'),1);
d2cmp=sort(unique(do(:,2),'rows'),1);

do1=sort(do1,1);
do2=sort(do2,1);

r1=d1cmp-do1;
r2=d2cmp-do2;

disp(sum(sum(r1)));
disp(sum(r2));

% test UniformGrid
set=UniformGrid('safedomain.scs');
set1=UniformGrid('safedomain.scs','projection',[1]);
set2=UniformGrid('safedomain.scs','projection',[2]);

pt=set.points;
pt1=set1.points;
pt2=set2.points;

pt1cmp=sort(unique(pt(:,1),'rows'),1);
pt2cmp=sort(unique(pt(:,2),'rows'),1);

pt1=sort(pt1,1);
pt2=sort(pt2,1);

r3=pt1cmp-pt1;
r4=pt2cmp-pt2;

disp(sum(r3));
disp(sum(r4));
end