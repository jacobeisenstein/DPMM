% Copyright (C) 2007 Jacob Eisenstein: jacobe at mit dot edu
% distributable under GPL, see README.txt

[Y,z,mu,ss,p] = drawGmm(200);
subplot(1,2,1);
title('generative clusters');
scatterMixture(Y,z);
params = dpmm(Y,100);
subplot(1,2,2);
title('dpmm clustering');
scatterMixture(Y,params(end).classes);
