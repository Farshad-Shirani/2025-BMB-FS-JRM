function reduceEdgeLines(q1_plotHandle, numEdgeLines_x1, numEdgeLines_x2 )

%%Extract X,Y and Z data from surface plot
x = q1_plotHandle.XData;
y = q1_plotHandle.YData;
z = q1_plotHandle.ZData;

%%Create vectors out of surface's XData and YData
x = x(1,:);
y = y(:,1);
%%Divide the lengths by the number of lines needed
xnumlines = numEdgeLines_x1; 
ynumlines = numEdgeLines_x2; 
xspacing = round(length(x)/xnumlines);
yspacing = round(length(y)/ynumlines);
%%Plot the mesh lines 
% Plotting lines in the X-Z plane
hold on
for i = 1:yspacing:length(y)
    Y1 = y(i)*ones(size(x)); % a constant vector
    Z1 = z(i,:);
    plot3(x,Y1,Z1,'-k');
end
% Plotting lines in the Y-Z plane
for i = 1:xspacing:length(x)
    X2 = x(i)*ones(size(y)); % a constant vector
    Z2 = z(:,i);
    plot3(X2,y,Z2,'-k');
end
hold off