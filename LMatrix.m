function M = LMatrix(N,h)

    
    I = eye(N-2);
    e = ones(N-2,1);
    T = spdiags([e -4*e e], [-1 0 1], N-2, N-2);
    S = spdiags([e e], [-1 1], N-2, N-2);
    M = (kron(I,T) + kron(S,I))/h^2;
end       