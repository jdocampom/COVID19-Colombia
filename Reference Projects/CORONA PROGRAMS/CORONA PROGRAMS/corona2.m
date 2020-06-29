function f=corona2(a,x)
f=(a(1)*tanh((x-a(2))/a(3)));
f=exp(f)-a(4);
