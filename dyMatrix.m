function M = dyMatrix(N,h,xSize)

    % Build d/dy matrix
    e = ones((N-2)^2,1);
    M = spdiags([-1*e e], [-xSize+2 xSize-2], (N-2)^2, (N-2)^2)/(2*h);
end    