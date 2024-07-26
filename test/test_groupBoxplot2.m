x1=[1,2,2,3,5,3];
x2=[2,5,4,5,8,6];
g1=[x1,x2];
v1=[repmat('A',1,numel(x1)),repmat('B',1,numel(x2))];
%group2
x3=[2,8,9,2,1,6];
x4=[5,4,3,22,11,6];
g2=[x3,x4];
v2=[repmat('A',1,numel(x3)),repmat('B',1,numel(x4))];
%group3
x5=[10,12,22,4];
x6=[12,15,4,25];
g3=[x5,x6];
v3=[repmat('A',1,numel(x5)),repmat('B',1,numel(x6))];
%%

G=[g1,g2,g3]
vg1 = [repmat('1',1,numel(v1)),repmat('2',1,numel(v2)),repmat('3',1,numel(v3))];
vg2=[v1,v2,v3] ;
%%
clc
figure;
boxplot(G', {vg1';vg2'}, 'factorseparator', 1, 'factorgap', 1,...
    'colorgroup',vg2','labelverbosity','majorminor');