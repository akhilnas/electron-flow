function f = gauss_distribution(x,a,b,c)
p1 = -.5 * ((x - b)/c) .^ 2;
f = exp(p1) .* a; 
end

