clc, clear; format shortG; format compact;

addpath(genpath(strcat(pwd,'\dace')));

load("ak_mcs_results.mat",'xTrain','S','initial_xTrain','krgMdl')
xTrain(1:size(initial_xTrain,1),:) = [];

x1 = linspace(-5,5,200);
x2 = linspace(-5,5,200);
[X,Y] = meshgrid(x1,x2);

Z=zeros(size(X));
Z_bar=zeros(size(X));
for i=1:size(X,1)
    for j=1:size(X,2)
        Z(i,j) = g_func([X(i,j),Y(i,j)]);
        Z_bar(i,j) = predictor([X(i,j),Y(i,j)], krgMdl);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot samples used in the active learning

f1 = figure;
set(f1,'units','inches','position',[1,1,5,5]);
leg_title = {};

s = scatter(S(:,1),S(:,2),'go','filled'); hold on; alpha(s,0.2); leg_title{1} = 'All samples';
scatter(xTrain(:,1),xTrain(:,2),'bo','filled'); hold on; leg_title{2} = 'Added samples';
scatter(initial_xTrain(:,1),initial_xTrain(:,2),100*ones(1,size(initial_xTrain,1)), ...
    'red','filled','pentagram'); hold on; leg_title{3} = 'Initial DOE';

contour(X,Y,Z,[0,0],LineColor="black",LineWidth=1); leg_title{4} = 'Exact g(x)';
contour(X,Y,Z_bar,[0,0],LineColor="magenta",LineWidth=1,LineStyle='--'); leg_title{5} = 'Surrogate g(x)';


xlabel('$x_1$',Interpreter='latex');
ylabel('$x_2$',Interpreter='latex');
legend(leg_title,Location="southeast");
box on


% Save figure content to PDF file
exportgraphics(gcf, 'example2.svg');

rmpath(genpath(strcat(pwd,'\dace')));