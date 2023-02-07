% use linear regression to identify neurons with ramping activity
% neurons with ramping acitivity: neural acitivity decrease or increase consistently
% across time. Context information can be decoded from neural activity of these
% neurons.
neuroactivity = neuroactivity_bftone;
x = 1:1:size(neuroactivity,2);% #trial
is_ramping = nan(size(neuroactivity,1),1);
Stats = nan(size(neuroactivity,1),4);

for unit = 1:size(neuroactivity,1)
    x1 = ones(size(x',1),1);%trial number, miss trial not included
    X = [x1 x'];
    y = neuroactivity(unit,:);%mean neural activity
    [~,~,~,~,stats] = regress(y',X);
    Stats(unit,:) = stats;
    is_ramping(unit) = stats(3) < 0.01;
end


%% randomly assign switch trials
repeat_times = 20;
num_trials = size(neuroactivity,2);
numneuron = size(neuroactivity,1);
switch_trials = randperm(num_trials-50, repeat_times) + 25;

type = zeros(1, num_trials);
randomAUC = zeros(repeat_times, numneuron);

for i = 1: repeat_times
    switch_trial = switch_trials(i);
    type_i = type;
    type_i(1:switch_trial) = 1;
    for jneuron = 1:numneuron
        [randomAUC(i, jneuron), ~] = cal_auc2(neuroactivity(jneuron,:), type_i);
    end
end

VarAUC = var(randomAUC);
dAUC = sum(abs(randomAUC - AUC_block(1,:)))/repeat_times;

%% the ratio of variance
%condition of 1 switch 
switch_trial = find(diff(BlockType(inds_use)));
bw_trials = 100;
bwblock = (switch_trial-bw_trials/2+1):(switch_trial+bw_trials/2);
block1 = 1:(switch_trial-bw_trials/2);
block2 = (switch_trial+bw_trials/2+1):num_trials;

var_inblock = zeros(1,numneuron);
var_bwblock = zeros(1,numneuron);

for jneuron = 1:numneuron
    %variance in blocks
    var_inblock(jneuron) = (var(neuroactivity(jneuron, block1))+var(neuroactivity(jneuron, block2)))/2;%in block1

    %variance between blocks
    var_bwblock(jneuron) = var(neuroactivity(jneuron, bwblock));
end

ratio_of_variance = var_inblock./var_bwblock;
%% overlap with auc significance?
is_blocksig = AUC_block(1,:)>block_sig_thre01;
num_ramping = sum(is_ramping);
num_blocksig = sum(is_blocksig);
num_overlap = sum(is_ramping' & is_blocksig);

figure;
scatter(AUC_block(1,:), Stats(:,1));
xlabel('AUC');
ylabel('R^2');
hold on;
%plot AUC sig threshold;
yl1 = ylim;
y_aucsig = linspace(yl1(1), yl1(2), 20);
x_aucsig = block_sig_thre01*ones(1,20);
plot(x_aucsig, y_aucsig,'b--');
hold off;
saveas(gcf,fullfile(np_data_pathway, 'example_neurons', 'R2.png'));

figure;
scatter(AUC_block(1,:), VarAUC);
xlabel('AUC');
ylabel('var(AUC)');
hold on;
%plot AUC sig threshold;
yl1 = ylim;
y_aucsig = linspace(yl1(1), yl1(2), 20);
x_aucsig = block_sig_thre01*ones(1,20);
plot(x_aucsig, y_aucsig,'b--');
hold off;
saveas(gcf,fullfile(np_data_pathway, 'example_neurons', 'var(AUC).png'));

figure;
scatter(AUC_block(1,:), dAUC);
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

figure;
scatter(AUC_block(1,:), ratio_of_variance);
xlabel('AUC');
ylabel('ratio of variance');
hold on;
%plot AUC sig threshold;
yl1 = ylim;
y_aucsig = linspace(yl1(1), yl1(2), 20);
x_aucsig = block_sig_thre01*ones(1,20);
plot(x_aucsig, y_aucsig,'b--');
hold off;
saveas(gcf,fullfile(np_data_pathway,'example_neurons', 'ratio of variance.png'));

%% plot example neuron
[~, unit_sig] = maxk(AUC_block(1,:),40);
%unit_sig = find(~is_ramping' & is_blocksig);
x_switch = find(diff(BlockType(inds_use)));

for iunit = 1:length(unit_sig)
    iinterest = unit_sig(iunit);

    figure;
    x = 1:sum(inds_use);
    y1 = smooth(neuroactivity(iinterest,:),5);
    plot(x,y1,'k','LineWidth',1);
    hold on;

    yl = ylim;
    y2 = linspace(yl(1),yl(2), 20);
    for iblock_change = 1:length(x_switch)
        plot(x_switch(iblock_change)*ones(20,1),y2,'k--','LineWidth',1);% block change
    end
    hold off;
    
    txtauc = ['auc = ', num2str(AUC_block(1,iinterest))];
    txtvar = ['var(auc) = ', num2str(VarAUC(iinterest))];
    txtdauc = ['distance(auc) = ', num2str(dAUC(iinterest))];
    txtr2 = ['R^2 = ', num2str(Stats(iinterest,1))];
    txtratio = ['ratio of variance = ', num2str(ratio_of_variance(iinterest))];
    
    text(300,320, txtr2, 'Units', 'pixels');
    text(300,300, txtvar, 'Units', 'pixels');
    text(300,280, txtdauc, 'Units', 'pixels');    
    text(300,260, txtauc, 'Units', 'pixels');
    text(300,240, txtratio, 'Units','pixels');
    
    xlabel('# Trial','FontSize',12);
    ylabel('Firing rate (Hz)','FontSize',12);
    
    pngname = ['neuron ', num2str(iinterest), '.png'];
    saveas(gcf,fullfile(np_data_pathway, 'example_neurons',  pngname));
    close;
end
