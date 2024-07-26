% test out 3d plot idea:

testx = 0:10;
testy = [0 1 2 3 4 5 5 4 4 3 0];

% figure; plot(testx, testy)

y1 = ones(11, 1);
z1 = testy;

x2 = testx; 
y2 = ones(11, 1) + 2; 
z2 = testy;


figure; plot3(x1, y1, z1, x2, y2, z2)



% turn regular plot data into step-function data:

x = 0:10:
y = [1 2 3 4 5 4 6 4 3 2 0];

xS = 0