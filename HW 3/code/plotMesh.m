function  plotMesh(IEN,x,fileName,isScatterOn)
%plotMesh Plots the mesh and prints it to file in eps format
%   IEN is the connectivity matrix
%   x is a column vector containing co-ordinates of all the nodes in the
%   mesh
%   fileName is the name of the file with or without extension
%   isScatterOn is a boolean flag to decide whether or not to do a scatter
%   plot for nodes. We turn it off for quadratic triangular elements

dim = 3;
x = reshape(x,dim,[]);
f = figure('visible','off');
trisurf(IEN(:,1:3),x(1,:),x(2,:),x(3,:));

if(isScatterOn)
    hold on;
    scatter3(x(1,:),x(2,:),x(3,:));
    hold off;
end

% Save the color eps file
print(f,fileName,'-depsc');

% Save the figure file
idx = strfind(fileName,'.');
if(isempty(idx))
    figFileName = strcat(fileName,'.fig');
else
    figFileName = strcat(fileName(1:idx-1),'.fig');
end

savefig(f,figFileName);

end

