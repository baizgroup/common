function [c, colormap] = create20colors()
% function returns different colors as array and as colormap. in the
% meantime it was extended with more than 20 colors however the old calls
% still work.
%
% example for usage:
% a = rand(10,5);
% [c, colormap] = create20colors();
% hold on;
% for ii = 1:5
%     plot(a(:,ii), 'color', c{ii});
% end
% 
% B. Knapp 2007-05-02

orange =    [0.9, 0.4, 0.1];
lightGray = [0.8, 0.8, 0.8];
darkRed =   [0.7, 0.1, 0.1];
pink =      [0.9, 0.0, 0.8];
lila =      [0.5, 0.0, 0.8];
lightLila = [0.8, 0.7, 0.8];
greenBlue = [0.0, 0.7, 0.7];
redBlue =   [0.7, 0.7, 0.0];
brown =     [0.7, 0.5, 0.5];
dirtyYellow=[0.9, 0.9, 0.5];
gold =      [0.9, 0.7, 0.2];
meadowGreen=[0.5, 0.7, 0.2];
bloodRed =  [0.9, 0.2, 0.4];
red =       [1.0, 0.0, 0.0];
green =     [0.0, 1.0, 0.0];
blue =      [0.0, 0.0, 1.0];
cyan =      [0.5, 1.0, 1.0];
mangenta  = [0.7    0.019608     0.73725];
yellow =    [1.0, 1.0, 0.0]; 
black  =    [0.0, 0.0, 0.0]; 
red1 =       [0.9, 0.0, 0.0];
red2 =       [0.8, 0.0, 0.0];
red3 =       [0.7, 0.0, 0.0];
red4 =       [0.6, 0.0, 0.0];
red5 =       [0.5, 0.0, 0.0];
red6 =       [0.4, 0.0, 0.0];
red7 =       [0.3, 0.0, 0.0];
red8 =       [0.2, 0.0, 0.0];
red9 =       [0.1, 0.0, 0.0];
green1 =     [0.0, 0.9, 0.0];
green2 =     [0.0, 0.8, 0.0];
green3 =     [0.0, 0.7, 0.0];
green4 =     [0.0, 0.6, 0.0];
green5 =     [0.0, 0.5, 0.0];
green6 =     [0.0, 0.4, 0.0];
green7 =     [0.0, 0.3, 0.0];
green8 =     [0.0, 0.2, 0.0];
green9 =     [0.0, 0.1, 0.0];
blue1 =      [0.0, 0.0, 0.9];
blue2 =      [0.0, 0.0, 0.8];
blue3 =      [0.0, 0.0, 0.7];
blue4 =      [0.0, 0.0, 0.6];
blue5 =      [0.0, 0.0, 0.5];
blue6 =      [0.0, 0.0, 0.4];
blue7 =      [0.0, 0.0, 0.3];
blue8 =      [0.0, 0.0, 0.2];
blue9 =      [0.0, 0.0, 0.1];
yellow1 =    [1.0, 1.0, 0.1]; 
yellow2 =    [1.0, 1.0, 0.2]; 
yellow3 =    [1.0, 1.0, 0.3]; 
yellow4 =    [1.0, 1.0, 0.4]; 
yellow5 =    [1.0, 1.0, 0.5]; 
yellow6 =    [1.0, 1.0, 0.6]; 
yellow7 =    [1.0, 1.0, 0.7]; 
yellow8 =    [1.0, 1.0, 0.8]; 
yellow9 =    [1.0, 1.0, 0.9]; 
cyan1 =      [0.1, 1.0, 1.0];
cyan2 =      [0.2, 1.0, 1.0];
cyan3 =      [0.3, 1.0, 1.0];
cyan4 =      [0.4, 1.0, 1.0];
cyan5 =      [0.5, 1.0, 1.0];
cyan6 =      [0.6, 1.0, 1.0];
cyan7 =      [0.7, 1.0, 1.0];
cyan8 =      [0.8, 1.0, 1.0];
cyan9 =      [0.9, 1.0, 1.0];





c = {blue, green, red, cyan, mangenta, yellow, black, orange, lightGray, darkRed, ...
          pink, lila, lightLila, greenBlue, redBlue, brown, dirtyYellow, gold, meadowGreen, bloodRed, ...
          red1, red2, red3, red4, red5, red6, red7, red8, red9, ...
          green1, green2, green3, green4, green5, green6, green7, green8, green9, ...
          blue1, blue2, blue3, blue4, blue5, blue6, blue7, blue8, blue9, ...
          yellow1, yellow2, yellow3, yellow4, yellow5, yellow6, yellow7, yellow8, yellow9, ...
          cyan1, cyan2, cyan3, cyan4, cyan5, cyan6, cyan7, cyan8, cyan9
          };
      
colormap(1,:) = blue;
colormap(2,:) = green;
colormap(3,:) = red;
colormap(4,:) = cyan;
colormap(5,:) = mangenta;
% colormap(6,:) = yellow;
% colormap(7,:) = black;
% colormap(8,:) = orange;
% colormap(9,:) = lightGray;
% colormap(10,:) = darkRed;
% colormap(11,:) = pink;
% colormap(12,:) = lila;
% colormap(13,:) = lightLila;
% colormap(14,:) = greenBlue;
% colormap(15,:) = redBlue;
% colormap(16,:) = brown;
% colormap(17,:) = dirtyYellow;
% colormap(18,:) = gold;
% colormap(19,:) = meadowGreen;
% colormap(20,:) = bloodRed;

