close all
clear all

spacings = [6.3, 9.6, 47.5, 88.5];
linestyles = {'b','k','r','g'};


figure
for spacing_idx = 1:length(spacings)
    spacing = spacings(spacing_idx);

    d = importdata(['Results/Slice_ER_' num2str(spacing) '.csv']);
    ca_er = d.data(:,1);
    fprintf('Max at ER for spacing %gnm is %g\n',spacing,max(ca_er));

    [n_er, ca_er_out] = hist(ca_er,1000);

    d = importdata(['Results/Slice_Membrane_' num2str(spacing) '.csv']);
    ca_mem = d.data(:,1);
    fprintf('Max at membrane for spacing %gnm is %g\n',spacing,max(ca_mem));

    [n_mem, ca_mem_out] = hist(ca_mem,1000);

    cumulative_n_er = zeros(length(n_er),1);
    cumulative_n_er(1) = n_er(1);
    for i=2:length(n_er)
        cumulative_n_er(i) = cumulative_n_er(i-1) + n_er(i);
    end

    cumulative_n_mem = zeros(length(n_mem),1);
    cumulative_n_mem(1) = n_mem(1);
    for i=2:length(n_mem)
        cumulative_n_mem(i) = cumulative_n_mem(i-1) + n_mem(i);
    end
    
    plot(ca_er_out,100 - 100.*cumulative_n_er./cumulative_n_er(end),[linestyles{spacing_idx} '-'],'LineWidth',2.0)
    hold on
    plot(ca_mem_out,100 - 100.*cumulative_n_mem./cumulative_n_mem(end),[linestyles{spacing_idx} '--'],'LineWidth',2.0)

end

legend('6.3nm spacing, at ER','6.3nm spacing, at outer membrane','9.6nm spacing, at ER','9.6nm spacing, at outer membrane',...
       '47.5nm spacing, at ER','47.5nm spacing, at outer membrane','88.5nm spacing, at ER','88.5nm spacing, at outer membrane','Location','NorthEast')
xlabel('[Ca] (uM)','FontSize',14)
ylabel('Proportion of membrane >= concentration (%)','FontSize',14)
ylim([0 100])
xlim([0 4])
set(gca,'FontSize',13)
