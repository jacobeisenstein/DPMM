function h = scatterMixture(data, z)
%function h = scatterMixture(data, z)
    colors = ['r','g','b','c','m','k','y'];
    symbols = ['.','o','x','+','*','s','d','v','^','<','>','p','h'];
    for i = min(z):max(z)
        hold on;
        format = strcat(colors(1+rem(i,numel(colors))),symbols(1+rem(3*i,numel(symbols))));
        scatter(data(z==i,1),data(z==i,2),format);
    end
    legend;
end