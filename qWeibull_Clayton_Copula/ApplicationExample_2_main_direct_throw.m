% global lifes;
% lifes = [0.058 0.070 0.090 0.105 0.113 0.121 0.153 0.159 0.224 0.421 0.570 0.596 0.618 0.834 1.019 1.104 1.497 2.027 2.234 2.372 2.433 2.505 2.690 2.877 2.879 3.166 3.455 3.551 4.378 4.872 5.085 5.272 5.341 8.952 9.188 11.399];

[Obj, sol, output]=ABC_direct_throw_local_search_adaptive(1);
objMean = mean(Obj);
objStd = std(Obj);
solMean = mean(sol);
solStd = std(sol);

