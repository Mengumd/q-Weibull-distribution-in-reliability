figure;
hold on;
t=0:0.1:20;
legendTextCell = {};
for q = [0.5 1.5]
    for kse = [0.5 1.5];
        sigma = 5;
        y = (kse/sigma) * (t/sigma).^(-kse-1) .* (exp_q((-1/(2-q))*(t/sigma).^(-kse),q));
        R = 1 - (1-((1-q)/(2-q))*(t./sigma).^(-kse)).^((2-q)/(1-q));
        h = y./R;
        plot(t,h,'-');
        legendTextCell{end+1} = ['q=' num2str(q) ',\xi=' num2str(kse)];
    end
end
legend(legendTextCell);
axis([0 20 0 0.5])



xlabel('t');
ylabel('h(t)');
title('hazard rate for q-Frechet');