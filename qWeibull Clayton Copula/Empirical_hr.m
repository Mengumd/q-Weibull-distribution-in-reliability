 count = 8;
    binWid = 12/count;
    H = histogram(lifes,count,'Normalization','count','BinWidth',binWid)
    W = histogram(lifes,count,'Normalization','cumcount','BinWidth',binWid)
    nRisk = 36-W.Values;nRisk = [36 nRisk(1:end-1)]
    H = histogram(lifes,count,'Normalization','count','BinWidth',binWid)
    
xlabel('Time to Failure');
ylabel('N');
title('Nubmer of Failures');

    hr=H.Values./nRisk./binWid;
    figure;bar(binWid*(1:count)-binWid/2,hr,1);
xlabel('Time to Failure');
ylabel('hr');
title('Empirical Hazard Rate');
colormap('default')