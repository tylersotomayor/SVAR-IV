function plotirf_partial(point,upper,lower,varnames,shockname, prc)


[N,horizon] = size(point);

nrow = ceil(sqrt(N));
ncol = ceil(N / nrow);

figure;
for i=1:N
    
    subplot(nrow, ncol, i)
    
    fill([0:horizon-1 fliplr(0:horizon-1)]' ,[upper(i, :)'; fliplr(lower(i, :))'],...
    [0 0.4470 0.7410],'EdgeColor','None'); hold on;
    plot(0:horizon-1,(point(i, :)),'-','LineWidth',1.5,'Color','k'); hold on;
    line(get(gca,'Xlim'),[0 0],'Color',[1 0 0],'LineStyle','--','LineWidth',1); hold off;


    ylabel(varnames{i}, 'FontSize', 16);   
    title(shockname, 'FontSize', 16);  
    xlim([0 horizon-1]);
    set(gca,'FontSize',16)
    
end

legend({strcat(num2str(prc), '% confidence bands'),'IRF'},'FontSize',16,'Orientation','Horizontal')


end