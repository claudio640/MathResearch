function M = dxMatrix(N,h,xSize)

   
    % Build d/dx matrix
    e = ones((N-2)^2,1);
    esp = e; esp2 = e;
    esp(1:xSize-2:end) = 0;
    esp2(N-2:xSize-2:end) = 0;
    M = spdiags([-1*esp2 esp], [-1 1], (N-2)^2, (N-2)^2)/(2*h);
    
end     