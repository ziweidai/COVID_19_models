function heatmap_cluster(data,rownames,colnames,datarange)
D=pdist(data);
tree=linkage(D,'average');
leafOrder_row=optimalleaforder(tree,D);
D=pdist(data');
tree=linkage(D,'average');
leafOrder_column=optimalleaforder(tree,D);
if length(rownames)>=1
	heatmap(data(leafOrder_row,leafOrder_column),colnames(leafOrder_column),...
		rownames(leafOrder_row),[],'GridLine',':','ShowAllTicks',true,'TickAngle',45,...
		'MinColorValue',datarange(1),'MaxColorValue',datarange(2));
else
	heatmap(data(leafOrder_row,leafOrder_column),[],[],[],'GridLine',':',...
		'ShowAllTicks',true,'TickAngle',45,'MinColorValue',datarange(1),'MaxColorValue',datarange(2));
end
end