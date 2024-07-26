data = rand(20,24)
month = repmat({'jan' 'feb' 'mar' 'apr' 'may' 'jun' 'jul' 'aug' 'sep' 'oct' 'nov' 'dec'},1,2);
simobs = [repmat({'sim'},1,12),repmat({'obs'},1,12)];
boxplot(data,{month,simobs},'colors',repmat('rb',1,12),'factorgap',[5 2],'labelverbosity','minor');


figure; boxplot(y, {g2, g1}, 'colors', 'rbrbrbrbrbrb', )

figure; boxplot(y, {g2, g1}, 'colors', 'rbrbrbrbrbrb', 'factorgap', [5, 2])



figure; boxplot(y, {g2, g1}, 'ColorGroup', g1, 'FactorGap', 10)

figure; boxplot(y, {g2, g1}, 'ColorGroup', g1, 'FactorSeparator', [2,1])


%% Tutorial on grouping boxplots of different sample sizes


x = [1,2,3,4,5,1,2,3,4,6]; % change this
group = [1,1,2,2,2,3,3,3,4,4]; % change this
positions = [1 1.25 2 2.25]; % change this
boxplot(x,group, 'positions', positions);

set(gca,'xtick',[mean(positions(1:2)) mean(positions(3:4)) ])
set(gca,'xticklabel',{'Direct care','Housekeeping'})

color = ['c', 'y', 'c', 'y'];
h = findobj(gca,'Tag','Box');
for j=1:length(h)
   patch(get(h(j),'XData'),get(h(j),'YData'),color(j),'FaceAlpha',.5);
end

c = get(gca, 'Children');

hleg1 = legend(c(1:2), 'Feature1', 'Feature2' );


%% My own version of the above



% x = [1,2,3,4,5,1,2,3,4,6]; % change this
group = [1,1,2,2,2,3,3,3,4,4]; % change this
positions = [1 1.25 2 2.25]; % change this
boxplot(y,{g2,g1}, 'positions', positions);

set(gca,'xtick',[mean(positions(1:2)) mean(positions(3:4)) ])
set(gca,'xticklabel',{'Direct care','Housekeeping'})

color = ['c', 'y', 'c', 'y'];
h = findobj(gca,'Tag','Box');
for j=1:length(h)
   patch(get(h(j),'XData'),get(h(j),'YData'),color(j),'FaceAlpha',.5);
end

c = get(gca, 'Children');

hleg1 = legend(c(1:2), 'Feature1', 'Feature2' );
