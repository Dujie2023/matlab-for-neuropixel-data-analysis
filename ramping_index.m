%% calculate AUC changes when switch trial changes
neuroactivity = neuroactivity_bftone;

num_trials = size(neuroactivity,2);
numneuron = size(neuroactivity,1);
switch_trial = find(diff(blocktype_use));
repeat_times = 20;% sampling times

%switch twice, but few trials in 3rd block
if length(switch_trial) == 2
    temp = find(diff(BlockType));
    inds_use_temp = inds_use;
    inds_use_temp(temp(2)+1:end) = 0;
    blocktype_use_temp = BlockType(inds_use_temp);
    
    switch_trial = switch_trial(1);
    num_trials = sum(inds_use_temp);
end

medium_trials = 50:5:num_trials-50;
trials_bwblock = switch_trial-51:switch_trial + 50;
[setswitch_trials, ~] = setdiff(medium_trials, trials_bwblock);
AUC_setswitch_repeat = zeros(length(setswitch_trials),numneuron, repeat_times);

type_sampling = zeros(100,1);
type_sampling(51:100) = 1;

for i = 1:repeat_times
    
    for iswitch = 1:length(setswitch_trials)
        sampling_block1_trials = randperm(setswitch_trials(iswitch),50);
        sampling_block2_trials = randperm(num_trials-setswitch_trials(iswitch),50)+setswitch_trials(iswitch);
        sampling_trials = [sampling_block1_trials sampling_block2_trials];
       
        for jneuron = 1:numneuron
            [AUC_setswitch_repeat(iswitch, jneuron, i), ~] = cal_auc2(neuroactivity(jneuron,sampling_trials), type_sampling);
        end
        
    end
    
end

AUC_setswitch = mean(AUC_setswitch_repeat,3);
histogram(AUC_setswitch_repeat(:));
%% calculate distance(AUC) to identify ramping neurons
trials_bwblock = (switch_trial-51):(switch_trial + 50);
[~, idx_AUCinblock] = setdiff(setswitch_trials, trials_bwblock);
AUC_inblock = AUC_setswitch(idx_AUCinblock,:);

distance_AUC = mean(AUC_inblock - AUC_block(1,:),1);

figure;
scatter(AUC_block(1,:), distance_AUC);
xlabel('AUC');
ylabel('distance(AUC)');

hold on;
%plot AUC sig threshold;
yl1 = ylim;
y_aucsig = linspace(yl1(1), yl1(2), 20);
x_aucsig = block_sig_thre01*ones(1,20);
plot(x_aucsig, y_aucsig,'b--');
hold off;
saveas(gcf,fullfile(np_data_pathway,'example_neurons', 'distance(AUC).png'));

%% plot 
[~, unit_sig] = maxk(AUC_block(1,:),20);

for iunit = 1:length(unit_sig)
    
    iinterest = unit_sig(iunit);
    
    figure;
    set(gcf, 'position', [200 200 1000 600]);
    subplot(1,2,1);
    x = 1:sum(inds_use);
    y1 = smooth(neuroactivity(iinterest,:),5);
    plot(x,y1,'k','LineWidth',1);
    hold on;

    %plot block change
    yl = ylim; 
    y2 = linspace(yl(1),yl(2), 20);
    plot(switch_trial*ones(20,1),y2,'k--','LineWidth',1);% block change
    hold off;

    xlabel('# Trial','FontSize',12);
    ylabel('Firing rate (Hz)','FontSize',12);
    txtauc = ['auc = ', num2str(AUC_block(1,iinterest))];
    txtdisauc = ['distance(auc) = ', num2str(distance_AUC(iinterest))];
    text(180,460, txtauc, 'Units', 'pixels');    
    text(180,480, txtdisauc, 'Units', 'pixels');    

    subplot(1,2,2);
    y3 =AUC_setswitch(:,iinterest);
    distance2switch = setswitch_trials - switch_trial;
    plot(distance2switch, y3);

    xlabel('Distance to switch trial (# Trial)','FontSize',12);
    ylabel('AUC','FontSize',12);
    
    pngname = ['neuron ', num2str(iinterest), '.png'];
    saveas(gcf,fullfile(np_data_pathway, 'example_neurons',  pngname));
    close;
end