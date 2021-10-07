% blap_ODE_system_proposal(endtimes, timestep, cell_density, protein_abund, blap_conc, p_change, x_change)
endtimes = 3600*2; % 2 hours in seconds
timestep = 1; % in seconds
cell_density = 1e9; % cells per Liter
protein_abund = [0.0520366765832331;0.0545386805308612;0.0857070042654188;0.0881078120842245;0.000440634256626082;0.00114502724318699;0.0127500286488247;0.0139948821104316;7.80621199066546;3.57979266612128;0.183810580839153;0.568659709462479;0.0963871431955497;1.14579927273508;0.0394839712424266;0.168540247916159;0.0124121404861996;0.00115046489521783;0.196125362492972;0.0120011757358293;0.239021944276831;0.514762104888601;0.130493812185125;2.62897341896851;0.221232056471700;0.276994204117248;1.68699919087255;0.150011683588273;0.0592678392422175;0.0310702399754362;0.000642710762887895;0.00144353090653194];
blap_conc = 3;
p_change = ones(35,1);
x_change = ones(31,1);

%% Establish genes of interest
genes = ["GPX1", "GPX2", "GPX3", "GPX4", "GPX5", "GPX6", "GPX7", "GPX8", "PRDX1", "PRDX2", ...
    "PRDX3", "PRDX6", "CAT", "TXN", "TXN2", "TXNRD1", "TXNRD2", "TXNRD3", "GLRX", "GLRX2", ...
    "G6PD", "GLUD1", "GSR", "GSTP1", "POR", "NQO1", "SOD1", "SOD2", "SOD3", "AQP3", ...
    "AQP8", "AQP9"];
    
%% Get files
metadata_file = pwd() + "\\Data\\2021-06-08_scrnaseqdata.csv";
files = dir(pwd() + "\\Data\\HNSCC_Broad_scRNAseq");
%%
metadata = readtable(metadata_file);
%%
gene_index = [];
protein = readtable(files(3).folder + "\\" + metadata.cell(1) + ".csv");
for i = 1:length(genes)
    gene_index = [gene_index, find(protein.Var1 == genes(i))];
end
%%
protein_vals = zeros(32,height(metadata));
%%
for i = 1:height(metadata)
    protein = readtable(files(3).folder + "\\" + metadata.cell(i) + ".csv");
    protein_vals(:,i) = table2array(protein(gene_index,2));
    
end
% Convert from PPM to micromolar
protein_vals = protein_vals.*5.29e-4;

%% Run control sim
protein_averages = mean(protein_vals, 2);

%%

num_x = 31;
% Sensitivities of all parameters
percentchange = 10.;
h2o2e_sensitivities = ones(35,2);
nadph_sensitivities = ones(35,2);
trx_sensitivities = ones(35,2);
gsh_sensitivities = ones(35,2);
%out_sensitivities = ones(34,2);
y_ctrl = blap_ODE_system_proposal(endtimes, timestep, cell_density, protein_abund, blap_conc, p_change, x_change);
%y_up = zeros(35,endtimes+1,num_x+1);
%y_down = zeros(35,endtimes+1,num_x+1);
h2o2e_ctrl = y_ctrl(end,2);
nadph_ctrl = y_ctrl(end,24)/y_ctrl(end,25);
trx_ctrl = y_ctrl(end,15)/y_ctrl(end,16);
gsh_ctrl = y_ctrl(end,7)/y_ctrl(end,8);
%out_ctrl = y_ctrl(7201,2);
%i = 1;
%disp(i)
num_params = 35;

%%
for param1 = 1:num_params
    %for kd_rate = 1:length(kd_rates)
    %disp(i)
    p_up = p_change;
    p_down = p_change;
    p_up(param1,1) = 1.1;
    p_down(param1,1) = 0.9;
    
    y_up =  blap_ODE_system_proposal(endtimes, timestep, cell_density, protein_abund, blap_conc, p_up, x_change);
    y_down =  blap_ODE_system_proposal(endtimes, timestep, cell_density, protein_abund, blap_conc, p_down, x_change);
    
    h2o2e_up = y_up(end,3);
    h2o2e_down = y_down(end,3);
    h2o2e_sensitivities(param1,1) = ((h2o2e_up-h2o2e_ctrl)/h2o2e_ctrl)/(percentchange/100);
    h2o2e_sensitivities(param1,2) = ((h2o2e_down-h2o2e_ctrl)/h2o2e_ctrl)/(percentchange/100);
    
    nadph_up = y_up(end,24)/y_up(end,25);
    nadph_down = y_down(end,24)/y_down(end,25);
    nadph_sensitivities(param1,1) = ((nadph_up-nadph_ctrl)/nadph_ctrl)/(percentchange/100);
    nadph_sensitivities(param1,2) = ((nadph_down-nadph_ctrl)/nadph_ctrl)/(percentchange/100);
    
    trx_up = y_up(end,15)/y_up(end,16);
    trx_down = y_down(end,15)/y_down(end,16);
    trx_sensitivities(param1,1) = ((trx_up-trx_ctrl)/trx_ctrl)/(percentchange/100);
    trx_sensitivities(param1,2) = ((trx_down-trx_ctrl)/trx_ctrl)/(percentchange/100);
    
    gsh_up = y_up(end,7)/y_up(end,8);
    gsh_down = y_down(end,7)/y_down(end,8);
    gsh_sensitivities(param1,1) = ((gsh_up-gsh_ctrl)/gsh_ctrl)/(percentchange/100);
    gsh_sensitivities(param1,2) = ((gsh_down-gsh_ctrl)/gsh_ctrl)/(percentchange/100);
end

%{
% Adjusting b-lap cycling rate
blap_ooms = -10:10;
sim = 1;
figure
set(gca,'YTickLabel',[]);
for blap_oom = blap_ooms
    y = blap_ODE_system_proposal(7200, timestep, 29, 10^blap_oom);
    subplot(2,1,1)
    semilogy(y(:,1)./60/timestep,y(:,2),...
    'LineWidth',4,...
    'Color',[.3 .3 .3])
    hold on
    subplot(2,1,2)
    semilogy(y(:,1)./60/timestep,y(:,3),...
    'LineWidth',4,...
    'Color','r')
    hold on
end
%}
%%
save = 1;
t_sens = 3600*2;
filename = 'Oct4_fig2.png';

% Use below for double parameter sensitivity
%{
names = cell(1,528);
name_ind = 1;

for i = 1:33
    new_start = i + 1;
    for j = new_start:33
        names{name_ind} = strcat('k',string(i),'+k',string(j));
        name_ind = name_ind+1;
    end
end
%}

% Use below for single parameter sensitivity
%%{
names = {'k1','k2','k3','k4','k5','k6','k7','k8','k9','k10','k11','k12', ...
    'k13','k14','k15','k16','k17','k18','k19','k20','k21','k22','k23', ...
    'k24','k25','k26','k27','k28','k29','k30','k31','k32','k33','k34', 'k35'};
%%}


for sens = 1:4
    % H2O2 Sensitivity
    if sens == 1
        Objective_low_sum = h2o2e_sensitivities(:,2);
        Objective_high_sum = h2o2e_sensitivities(:,1);
        filename = 'h2o2e_sens.png';
        plot_title = 'Intracellular H2O2 Sensitivities';
        
        %{
        figure
        for i = 1:33
            plot(1:7201,y_up(i,:,3))
            hold on
        end
        title('Intracellular H2O2')
        saveas(gcf,"H2O2TimeSeries.png")
        %}
    end
    
    % NADPH sens
    if sens == 2
        Objective_low_sum = nadph_sensitivities(:,2);
        Objective_high_sum = nadph_sensitivities(:,1);
        filename = 'nadph_sens.png';
        plot_title = 'Intracellular NADPH/NADP+ Sensitivities';
        
        %{
        figure
        for i = 1:33
            plot(1:7201,y_up(i,:,24)./y_up(i,:,25))
            hold on
        end
        title('NADPH/NADP+ Ratio')
        saveas(gcf,"NADPH-NADP+TimeSeries.png")
        %}
    end
    
    % GSH sens
    if sens == 3
        Objective_low_sum = gsh_sensitivities(:,2);
        Objective_high_sum = gsh_sensitivities(:,1);
        filename = 'gsh_sens.png';
        plot_title = 'Intracellular GSH/GSSG Sensitivities';
        
        %{
        figure
        for i = 1:33
            plot(1:7201,y_up(i,:,7)./y_up(i,:,8))
            hold on
        end
        title('GSH/GSSG Ratio')
        saveas(gcf,"GSH-GSSGTimeSeries.png")
        %}
    end
    
    % Trx sens
    if sens == 4
        Objective_low_sum = trx_sensitivities(:,2);
        Objective_high_sum = trx_sensitivities(:,1);
        filename = 'Trx_sens.png';
        plot_title = 'Intracellular rTrx/oTrx Sensitivities';
        
        %{
        figure
        for i = 1:33
            plot(1:7201,y_up(i,:,15)./y_up(i,:,16))
            hold on
        end
        title('Trx-SH/Trx-SS Ratio')
        saveas(gcf,"rTrx-oTrxTimeSeries.png")
        %}
    end
    % Sort the values based on the lower change
    % Sort the higher values and the names arrays
    %    using the same indices
    [~,ind]=sort(Objective_high_sum,'ascend');
    Objective_low_sum = Objective_low_sum(ind);
    Objective_high_sum = Objective_high_sum(ind);
    names_Objective=names(ind);
   
    
    % Create a figure and plot the low and high horizontally
    figure('Position', [10 10 800 1500])
    h = barh(Objective_high_sum);
    hold on
    xmin=max([max(Objective_low_sum),max(Objective_high_sum),-min(Objective_low_sum),-min(Objective_high_sum)]);
    xmax=max([max(Objective_low_sum),max(Objective_high_sum),-min(Objective_low_sum),-min(Objective_high_sum)]);
    xlim([-1.025*xmax 1.025*xmax])
    barh(Objective_low_sum,'r')
    bh = get(h,'BaseLine');
    %set(bh,'BaseValue',in_ctrl);
    title(plot_title);


    set(gca,'yticklabel',names_Objective)
    set(gca,'Ytick',[1:length(names_Objective)],'YTickLabel',[1:length(names_Objective)])
    set(gca,'yticklabel',names_Objective)

    % get the current tick labels
    ticklabels = get(gca,'YTickLabel');
    % prepend a color for each tick label
    ticklabels_new = cell(size(ticklabels));

    
    for i = [3, 4, 5, 13, 20, 24, 25, 26]
        ordered_i = find(ind==i);
        ticklabels_new{ordered_i} = ['\color{red} ' ticklabels{ordered_i}];
    end
    for i = [8, 9, 10, 11, 12, 21, 27, 28]
        ordered_i = find(ind==i);
        ticklabels_new{ordered_i} = ['\color{blue} ' ticklabels{ordered_i}];
    end
    for i = [7, 22]
        ordered_i = find(ind==i);
        ticklabels_new{ordered_i} = ['\color[rgb]{0, 0.8, 0} ' ticklabels{ordered_i}];
    end
    for i = [14, 15, 16, 17, 18, 19]
        ordered_i = find(ind==i);
        ticklabels_new{ordered_i} = ['\color[rgb]{1,0.6,0.2} ' ticklabels{ordered_i}];
    end
    for i = [29, 30, 31, 32, 33, 34, 35]
        ordered_i = find(ind==i);
        ticklabels_new{ordered_i} = ['\color{black} ' ticklabels{ordered_i}];
    end
    for i = [1, 2, 6, 23]
        ordered_i = find(ind==i);
        ticklabels_new{ordered_i} = ['\color{magenta} ' ticklabels{ordered_i}];
    end
    % set the tick labels
    set(gca, 'YTickLabel', ticklabels_new);
    
    %xlabel('colors for legend')
    l1 = plot([NaN,NaN], 'r');
    l2 = plot([NaN,NaN], 'b');
    l3 = plot([NaN,NaN], 'color', '#00CC00');
    l4 = plot([NaN,NaN], 'color', '#FF9933');
    l5 = plot([NaN,NaN], 'k');
    l6 = plot([NaN,NaN], 'm');
    [~,hObj]=legend([l1, l2, l3, l4, l5, l6], {'GSH', 'Prx', ...
        'NADPH', 'Protein Thiol', 'Drug Metabolism', 'Other'},...
        'Location', 'SouthWest');           % return the handles array
    hL=findobj(hObj,'type','line');  % get the lines, not text
    set(hL,'linewidth',1.5)
    

    saveas(gcf,filename)
    end
%Objective_low_sum = (((y_down(:,7201,3)-in_ctrl)./in_ctrl)./.1)';
%Objective_high_sum = (((y_up(:,7201,3)-in_ctrl)./in_ctrl)./.1)';
%}
