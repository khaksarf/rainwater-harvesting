function [newCoor, X, Y] = rotation(x, y, theta)
angle = theta/180*pi;
newCoor = [ ...
    cos(angle) -sin(angle)
    sin(angle) cos(angle)
    ] * [x y]';
X = newCoor(1); Y = newCoor(2);

end