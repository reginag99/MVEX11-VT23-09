<<<<<<< HEAD
clear,clc

dark_green = [30,120,130]/255;
dark_orange = [112,81,28]/255;
dark_purple = [101,120,163]/255;

light_green = [122,214,185]/255;
light_orange = [252,151,98]/255;
light_purple = [181,220,243]/255;

x = linspace(0,10);
y1 = 3+x;
y2 = 5+x; 


subplot(2,2,1)
plot(x,y1, 'color',light_green)
hold on
plot(x,y2, 'color',dark_green)
hold on

subplot(2,2,2)
plot(x,y1, 'color',light_orange)
hold on
plot(x,y2, 'color',dark_orange)
hold on

subplot(2,2,3)
plot(x,y1, 'color',light_purple)
hold on
plot(x,y2, 'color',dark_purple)
hold on


=======
clear,clc

dark_green = [30,120,130]/255;
dark_orange = [112,81,28]/255;
dark_purple = [101,120,163]/255;

light_green = [122,214,185]/255;
light_orange = [252,151,98]/255;
light_purple = [181,220,243]/255;

x = linspace(0,10);
y1 = 3+x;
y2 = 5+x; 


subplot(2,2,1)
plot(x,y1, 'color',light_green)
hold on
plot(x,y2, 'color',dark_green)
hold on

subplot(2,2,2)
plot(x,y1, 'color',light_orange)
hold on
plot(x,y2, 'color',dark_orange)
hold on

subplot(2,2,3)
plot(x,y1, 'color',light_purple)
hold on
plot(x,y2, 'color',dark_purple)
hold on


>>>>>>> fd47580eac778269c6a64a6efd89ecec96a69860
