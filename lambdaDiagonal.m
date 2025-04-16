function L = lambdaDiagonal(lambda,m,n)

    L = spdiags(lambda*ones(m*n), 0, m*n, m*n);

end