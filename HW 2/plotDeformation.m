function plotDeformation(X,x)
%PLOTDEFORMATION plots the trianglular element's reference and deformed
%configurations

X_plot = [X(:,[1,2,3]) X(:,1)];
plot(X_plot(1,:),X_plot(2,:),'g')
hold on;
for i=1:3
    [xaf,yaf] = ds2nfu(X(1,i),X(2,i));
    str = sprintf('X%d',i);
    annotation('textbox', [xaf,yaf,0.1,0.1],...
        'String',str);
end

x_plot = [x(:,[1,2,3]) x(:,1)];
plot(x_plot(1,:),x_plot(2,:),'r');
axis square;
end