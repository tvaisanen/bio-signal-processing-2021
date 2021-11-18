

clf


x = 1:10;
x2 = -0.5:1:9.5;

y = 0.1:0.1:1;

size(x)
size(y)

plot(x,y) 
hold on;
plot(x2, interp1(x,y,x2, 'linear','extrap'),'o');
