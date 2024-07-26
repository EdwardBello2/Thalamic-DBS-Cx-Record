% script for fitting a sigmoid to the frequency sweep data



% first run the analysis script below:
script_runEntropyDirectISIord1_All_userVarH;
% outputs 2 variables Entropy_AllNeu and Entropy_grpLabel


% convert Entropy_grpLabel to numerical vector:
grp_Freq = cellfun(@str2num, Entropy_grpLabel);
grp_FreqLog = log10(grp_Freq);

% Now run cftool for sigmoid-fitting!
% cftool

% figure; scatter(grp_Freq, Entropy_AllNeu)
% hold on;
% 
% 
% x = -100:0.01:130;
% 
% % custom sigmoid function:
% a = -20;
% b = 10;
% c = 0.1;
% sh = 50;
% 
% y = (a ./ (exp(-c*(x - sh)) + 1)) + b;
% 
% 
% figure; plot(x,y);
% 


%%  test the look of a fit:
figure; scatter(grp_Freq, Entropy_AllNeu)
hold on

x = -10:0.01:130;

a = -20;
b = 1.774;
c = 0.96491;
sh = 50;

y = (a ./ (exp(-c*(x - sh)) + 1)) + b;


plot(x,y);
% 
% 
% function y = mySigmoid(x)
% a = -1;
% b = 3;
% c = 1;
% sh = 5;
% 
% y = (a ./ (exp(-c*(x - sh)) + 1)) + b;
% 
% end
