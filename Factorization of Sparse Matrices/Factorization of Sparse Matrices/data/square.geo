cl = 0.3;
Point(1) = {0, 0, 0, cl};
Point(2) = {1, 0, 0, cl};
Point(3) = {1, 1, 0, cl};
Point(4) = {0, 1, 0, cl};
Line(1) = {1, 2};
Line(2) = {3, 4};
Line(3) = {2, 3};
Line(4) = {1, 4};
Line Loop(6) = {1, 3, 2, -4};
Plane Surface(6) = {6};