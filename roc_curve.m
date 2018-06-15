function [auc,sn,sp] = roc_curve(deci,label_y,colour)
[val,Index] = sort(deci,'descend');
roc_y = label_y(Index);
%%
%roc
stack_x = cumsum(roc_y ==0)/sum(roc_y ==0);
stack_y = cumsum(roc_y == 1)/sum(roc_y == 1);
xlabel({'False Positive Rate'})
ylabel('True Positive Rate');
plot(stack_x,stack_y,'Color',colour,'LineWidth',1);
auc=sum((stack_x(2:length(roc_y),1)-stack_x(1:length(roc_y)-1,1)).*stack_y(2:length(roc_y),1));
sp=1-stack_x;
sn=stack_y;
 title(['AUC = ' num2str(auc) ]);
end