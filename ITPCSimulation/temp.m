
p = 0:1/1000:pi/4;
p1 = 0:1/100:pi/4;
p2 = 0:1/1000:pi/4;
abs(mean(exp(i*p1)))
abs(mean(exp(i*p2)))

p1 = rand(1,100)*pi/4;
p2 = rand(1,10000)*pi/4;
abs(mean(exp(i*p2)))

abs(mean(exp(i*p1)))
p1 = rand(1,100)*pi/4+pi/2;
p2 = rand(1,10000)*pi/4+pi/2;
abs(mean(exp(i*p1)))
abs(mean(exp(i*p2)))
p1 = rand(1,10)*pi/4+pi/2;
p2 = rand(1,1000)*pi/4+pi/2;
abs(mean(exp(i*p2)))
abs(mean(exp(i*p1)))
p2 = rand(1,1000)*pi/4;
p1 = rand(1,10)*pi/4;
abs(mean(exp(i*p1)))
abs(mean(exp(i*p2)))
p1 = rand(1,10)*pi;
abs(mean(exp(i*p1)))
p2 = rand(1,1000)*pi;
abs(mean(exp(i*p2)))
p1 = rand(1,300)*2*pi;
p2 = rand(1,600)*2*pi;
abs(mean(exp(i*p2)))
abs(mean(exp(i*p1)))