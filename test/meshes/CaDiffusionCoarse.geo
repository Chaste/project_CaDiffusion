cl__1 = 0.5;
Point(7) = {0, 0, 0, 0.5};
Point(8) = {10, 0, 0, 0.5};
Point(9) = {-10, 0, 0, 0.5};
Point(10) = {0, 10, 0, 0.5};
Point(11) = {0, -10, 0, 0.5};
Point(12) = {10, 0, 1.5, 0.5};
Point(13) = {0, 0, 1.5, 0.5};
Point(14) = {0, 10, 1.5, 0.5};
Point(19) = {-10, 0, 1.5, 0.5};
Point(24) = {0, -10, 1.5, 0.5};
Circle(5) = {8, 7, 10};
Circle(6) = {10, 7, 9};
Circle(7) = {9, 7, 11};
Circle(8) = {11, 7, 8};
Circle(12) = {12, 13, 14};
Circle(13) = {14, 13, 19};
Circle(14) = {19, 13, 24};
Circle(15) = {24, 13, 12};
Line(17) = {8, 12};
Line(18) = {10, 14};
Line(22) = {9, 19};
Line(26) = {11, 24};
Line Loop(10) = {5, 6, 7, 8};
Plane Surface(10) = {10};
Line Loop(19) = {5, 18, -12, -17};
Ruled Surface(19) = {19};
Line Loop(23) = {6, 22, -13, -18};
Ruled Surface(23) = {23};
Line Loop(27) = {7, 26, -14, -22};
Ruled Surface(27) = {27};
Line Loop(31) = {8, 17, -15, -26};
Ruled Surface(31) = {31};
Line Loop(32) = {12, 13, 14, 15};
Plane Surface(32) = {32};
Surface Loop(1) = {31, 19, 27, 32, 10, 23};
Volume(1) = {1};

Physical Volume(1) = {1};
Physical Surface(2) = {31, 19, 27, 32, 10, 23};
