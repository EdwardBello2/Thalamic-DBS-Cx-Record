% script for plotting out microelectrode thalamic mapping track

%% Display vertical ticks for each depth-recording of interest for microstim

stimTrack = []; % manually fill with the Depth and stim uA columns

stimTicks = stimTrack((stimTrack(:,2) > 0),1)
nTicks = length(stimTicks);
figure;
for iTick = 1:nTicks
    plot([0 1], [stimTicks(iTick), stimTicks(iTick)], ...
        'Color', 'k')
    hold on;
    
    
end
set(gca, 'XLim', [0 10], 'YLim', [20 40]);



%% Display vertical ticks for each depth-recording of interest for motor exam

motExTrack = [];

motExTicks = motExTrack((motExTrack(:,2) > 0),1)
nTicks = length(motExTicks);
figure;
for iTick = 1:nTicks
    plot([0 1], [motExTicks(iTick), motExTicks(iTick)], ...
        'Color', 'k')
    hold on;
    
    
end
set(gca, 'XLim', [0 10], 'YLim', [20 40]);

