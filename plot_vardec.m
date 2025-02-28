function plot_vardec(fevd, varnames, shocknames)

% Function to plot the FEVD

% Inputs:   fevd        = N x N x horizon array of forecast error shares 
%           varnames    = names of variables  
%           shocknames  = names of shocks

% Returns:  Plots of FEVD

if size(shocknames,2) > 1
    [N,~,horizon] = size(fevd);
else
    [N,horizon] = size(fevd);
end

for i=1:N
    
    subplot(1,N,i)
    if size(shocknames,2) > 1
        area(0:horizon-1, squeeze(fevd(i,:,:))','LineWidth',1); hold on; 
    else 
        othershocks = 1 - (fevd(i,:))';
        area(0:horizon-1, [(fevd(i,:))', othershocks],'LineWidth',1); hold on; 
    end
    set(gca,'FontSize',15)
    title(varnames{i})
    ylabel('% FEVD', 'FontSize', 18);
    xlim([0 horizon-1]);
    ylim([0 1])
          
end
    if size(shocknames,2) > 1
        legend(shocknames)
    else
        shocknames = [shocknames, {'Other Shocks'}];
        legend(shocknames)
    end
end