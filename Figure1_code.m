%% Figure 1 Code

clear
clc
load ../mat' files'/Asynch_Lambda.mat
load ../mat' files'/pos.mat
load ../mat' files'/active_origin_list

%% General statistical and stuctural analysis for Figure 1

All_rep_profiles = async_reads(:,cell2mat(Sorted_Names_all_strains(:,3)));

All_rep_profiles_cell = cell(size(Sorted_Names_all_strains,1),1);
for i=1:size(Sorted_Names_all_strains,1)
    All_rep_profiles_cell{i} = async_reads(:,Sorted_Names_all_strains{i,3});
end

Chromosomes = zeros(16,2);
Chr_profiles = cell(16,1);
Chr_profiles_single = cell(16,1);
Chr_pos = cell(16,1);
Chr_profiles_cell = cell(16,1);
position_vector = cell(16,1);
width_vector = cell(16,1);
All_profiles_per_chr = cell(length(All_rep_profiles_cell),1);
for c = 1:16
    current = find(pos(:,1)==c);
    Chromosomes(c,1) = current(1);
    Chromosomes(c,2) = current(end);
    Chr_profiles{c,:} = Mean_profiles(current(1):current(end),:); 
    Chr_pos{c} = pos(current(1):current(end),2:3);
    position_vector{c} = Chr_pos{c}*[0.5;0.5];
    for k=2:length(position_vector{c})
        width_vector{c}(1,1) = position_vector{c}(1);
        width_vector{c}(k,1)= position_vector{c}(k)-position_vector{c}(k-1);
    end
    for i=1:length(All_rep_profiles_cell)
        All_profiles_per_chr{i} = All_rep_profiles_cell{i,1}(current(1):current(end),:);  
    end
    Chr_profiles_cell{c} = All_profiles_per_chr;
end
oris = zeros(16,2);
oris_pos = cell(16,1);
for c=1:16
    current = find(origins(:,1)==c);
    oris(c,1) = current(1);
    oris(c,2) = current(end);
    oris_pos{c} = origins(current(1):current(end),2);
end

Lambdas = (Asynch_Lambda(indices,2));

Single_profiles = async_reads(:,indices);
n=floor(12000000/size(pos,1)); %No of bp per position
Autocorrelation = nan(567,length(indices));
for i = 1:length(indices)
Autocorrelation(:,i) = autocorr(Single_profiles(:,i),n);
end

low_auto_ind = NaN(length(indices),1);
for i=1:length(indices)
low_auto_ind(i) = find((Autocorrelation(:,i)<=0.5),1);
end

Mean_WT = nanmean(async_reads(:,ref_indices),2);
Mean_Rtt = nanmean(async_reads(:,Rtt_indices),2);
Mean_profiles = [Mean_WT,Mean_Rtt];
Mean_WT_Lambda = nanmean(cell2mat(Lambdas(1:9)));
Mean_Rtt_Lambda = nanmean(cell2mat(Lambdas(10:18)));
Mean_Lambdas = [Mean_WT_Lambda,Mean_Rtt_Lambda];
Lambdas = [Lambdas(1:9),Lambdas(10:18)];
Lambda_Std = nanstd(cell2mat(Lambdas));
Lambda_SE = nan(length(Lambda_Std),1);
for j=1:length(Lambda_Std)
Lambda_SE(j) = Lambda_Std(j)/sqrt(length(Lambdas(:,j)));
end

n=floor(12000000/size(pos,1)); %No of bp per position
Autocorrelation = nan(567,size(Mean_profiles,2));
for i = 1:size(Mean_profiles,2)
Autocorrelation(:,i) = autocorr(Mean_profiles(:,i),n);
end

low_auto_ind_av = NaN(size(Mean_profiles,2),1);
for i=1:size(Mean_profiles,2)
low_auto_ind_av(i) = find((Autocorrelation(:,i)<=0.5),1);
end

low_auto_ind_WT = low_auto_ind(1:9);
low_auto_ind_Rtt = low_auto_ind(10:18);
low_auto_ind = [low_auto_ind_WT,low_auto_ind_Rtt];

low_auto_ind_Std = nanstd(low_auto_ind);
low_auto_ind_SE = nan(size(low_auto_ind,2),1);
for i= 1:size(low_auto_ind,2)
    low_auto_ind_SE(i) = low_auto_ind_Std(i)/sqrt(length(low_auto_ind(:,i)));
end

%% Panel A - Repication profiles 

for i=16
    chr_len=max(Chr_pos{i,1}(:,2));
    Profiles_cell_Rtt = Chr_profiles{i,1}(:,Rtt); %for Rtt109
    Profiles_cell_Wt_1 = Chr_profiles{i,1}(:,WT); %for WT
    figure('units','normalized','outerposition',[0 0 1 1])
    n = length(oris_pos{i,1});
    for j = [1,3:6,8,11,19,20]%1:n %These are the late origin lines
        line([oris_pos{i,1}(j),oris_pos{i,1}(j)],[-1.5,1.5],'color',[0.5,0.5,0.5, .2],'linewidth',3,'HandleVisibility','off');
    end
    for j = [2,7,9,10,12:18]%1:n %These are the early origin lines
        line([oris_pos{i,1}(j),oris_pos{i,1}(j)],[-1.5,1.5],'color',[0.396, 0.364, 0.364],'linewidth',3,'HandleVisibility','off');
    end
    hold on
    dummyh = line(nan, nan, 'Linestyle', 'none', 'Marker', 'none', 'Color', 'none');
    dummyh1 = line(nan, nan, 'Linestyle', 'none', 'Marker', 'none', 'Color', 'none');
    plot(position_vector{i,1},(Profiles_cell_Wt_1)','color',[0.670, 0.654, 0.627],'linewidth',4);
    plot(position_vector{i,1},(Profiles_cell_Rtt)','color',[0.050, 0.392, 0.972],'linewidth',4);
    
    set(gca,'XTick',0:100000:chr_len,'fontsize',24,'linewidth',2);
    axis([-2,chr_len,-1.5,1.5]);
    set(gcf,'color','w');
    ylim([-0.8 0.7]);
    xlabel('Chromosome 16','fontsize',28);
    ylabel('Relative DNA content','fontsize',28);
    legend({'\color[rgb]{0.670, 0.654, 0.627} WT',...
        '\color[rgb]{0.050, 0.392, 0.972} \Delta Rtt109'},...
        'fontsize',20,'Location','northwest');
    legend boxoff
    hold off
end

%%  Panel B - Autocorrelation

figure
dummyh2 = line(nan, nan, 'Linestyle', 'none', 'Marker', 'none', 'Color', 'none');
dummyh3 = line(nan, nan, 'Linestyle', 'none', 'Marker', 'none', 'Color', 'none');
dummyh4 = line(nan, nan, 'Linestyle', 'none', 'Marker', 'none', 'Color', 'none');
hold on
plot((0:n)*0.5,(Autocorrelation(:,1)'),'color',[0.670, 0.654, 0.627],'linewidth',4);
plot((0:n)*0.5,(Autocorrelation(:,2)'),'color',[0.050, 0.392, 0.972],'linewidth',4);

set(gca,'fontsize',24,'linewidth',2); 
set(gcf,'color','w');
xlim([0 120]);
ylim([0 1]);
legend({'\color[rgb]{0.670, 0.654, 0.627} WT',...
    '\color[rgb]{0.050, 0.392, 0.972} Rtt109'},'fontsize',20); 
legend boxoff 
xlabel('Distance (kb)','fontsize',28);
ylabel('Autocorrelation','fontsize',28);
hold off

%% Panel C - Replicon length

figure
bar_data = flipud(Lambda([WT,Rtt],1));
b = barh(bar_data,'edgealpha',0);
b.FaceColor = 'flat';
b.CData(1,:) = [0.670, 0.654, 0.627];
b.CData(2,:) = [0.050, 0.392, 0.972];
hold on
herrorbar(bar_data,1:2,flipud(Lambda_SE([WT,Rtt],1)),'.');
set(gca, 'YTick', 1:2, 'YTickLabel',{'Rtt109','WT'},'FontSize',16);
xlabel('\lambda (kb)','FontSize',28);
set(gcf,'color','w');
box off

%% Panel D - Replicon length and autocorrelation

figure
colors = [0.670, 0.654, 0.627;0.050, 0.392, 0.972]; 
scatter(low_auto_ind_av.*0.5,Mean_Lambdas,100,colors,'O','filled'); 
hold on
text(low_auto_ind_av.*0.5,Mean_Lambdas,{'WT','Rtt109'},'FontSize',20);
errorbar(low_auto_ind_av.*0.5,Mean_Lambdas,Lambda_SE,'.');
herrorbar(low_auto_ind_av.*0.5,Mean_Lambdas,low_auto_ind_SE.*0.5,'.');
set(gca,'fontsize',24,'linewidth',2);
set(gca,'YScale','log') 
ylim([40 110]);
set(gcf,'color','w');
ylabel('Replicon length (kb)','fontsize',20);
xlabel('Autocorrelation (kb)','fontsize',20);