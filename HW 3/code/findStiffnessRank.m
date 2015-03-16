function rankOut = findStiffnessRank(K)
% FINDSTIFFNESSRANK accepts a 4-dimensional stiffness array as input
% and returns the rank of a matrix obtained by unrolling the array to a 2D
% matrix.

%% Unroll K

[n1,n2,n3,n4] = size(K);

K_unrolled = zeros(n1*n2,n3*n4);

for i=1:n1
    for a=1:n2
        row = sub2ind([n1,n2],i,a);
        for k=1:n3
            for b=1:n4
                col = sub2ind([n3,n4],k,b);
                K_unrolled(row,col) = K(i,a,k,b);
            end
        end
        
    end
end

%% Find the rank
%rankOut = rank(K_unrolled);
[~,eigValues] = eig(K_unrolled);

% Find number of non-zero eigen values
tol = norm(K_unrolled)*10^-10;
rankOut = 0;
% fprintf('The eigen values: \n');
for k=1:size(K_unrolled,1)
    %fprintf('%17.16f ',abs(eigValues(k,k)));
    if abs(eigValues(k,k)) > tol
        rankOut = rankOut+1;
    end
end
 %fprintf('%17.16f\n\n');
end