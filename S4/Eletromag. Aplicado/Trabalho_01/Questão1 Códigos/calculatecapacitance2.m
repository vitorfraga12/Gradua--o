function C = calCap(d, l, x0, N, EO)
    RHO = calR(d, l, x0, N, EO);
    DL = l / N;
    Q = sum(RHO(1:N)) * DL;
    C = abs(Q) / 2;
end