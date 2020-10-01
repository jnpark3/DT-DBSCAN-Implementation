% MATLAB code for cluster visualization
% Author: Jian Park (21PARKJ@sagehillschool.org)

dsource = 'cluster.out';
M = csvread( dsource );

X = M(:,1);
Y = M(:,2);
C = M(:,3);

scatter(X,Y,10, C);

outfile = strcat( dsource, '.png' );
saveas(gcf, outfile );
