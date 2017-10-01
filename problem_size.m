cheetah_dense = (9^3)*(12^3)
cheetah_sparse = 9 * (24^3)


[horizon, inputs] = meshgrid(1:20,1:20);
dense = (inputs.*horizon).^3;
sparse_12 = horizon.*(inputs + 12).^3;
sparse_24 = horizon.*(inputs + 24).^3;
ratio = dense./sparse_12;
ratio_24 = dense./sparse_24;
dense_better = ratio < 1;
dense_better_24 = ratio_24 < 1;
contourf(double(dense_better + dense_better_24),2);
title('Dense Better? 12 States and 24 States');
xlabel('horizon length');
ylabel('number of inputs');

set(gcf,'color','w');

export_fig 'problemsize1.eps'