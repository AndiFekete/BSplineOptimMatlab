%Example for demonstrating the nullad algorithm
load('beat');
y=beat;
x=1:length(beat);
[knots,aprx,err]=gyorsnullad(x,y,20,true);
