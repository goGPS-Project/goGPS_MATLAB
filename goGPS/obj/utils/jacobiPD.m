function P = jacobiPD(n,a,b,x)
    P = 0;
    for  s = 0:n
        P = P + ...
            1/(factorial(s)*factorial(n+a-s)*factorial(b+s)*factorial(n-s)) * ...
            ((x-1)/2).^(n-s).*((x+1)/2).^s;
    end
    P = P*factorial(n+a) * factorial(n+b);
end
