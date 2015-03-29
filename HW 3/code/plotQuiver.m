function fig = plotQuiver(fig,force,scale,IEN,X)
%PLOTQUIVER plots the force vectors on a figure handle 'fig 'at centroids
%of triangular elements given by IEN and X. X is assumed to be
%column vectors that need to be reshaped. Force is a numEle-by-3 matrix

ax = fig.CurrentAxes;

X = reshape(X,3,[]);
X = X';

% Calculate the centroids
X1 = X(IEN(:,1),:);
X2 = X(IEN(:,2),:);
X3 = X(IEN(:,3),:);
Centroids = 1/3*(X1+X2+X3);



% Do the quiver
hold(ax,'on');
quiver3(ax,force(:,1),force(:,2),force(:,3),...
    Centroids(:,1),Centroids(:,2),Centroids(:,3),scale);
hold(ax,'off');

end