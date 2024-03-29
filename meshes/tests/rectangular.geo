cl__1 = 1.0;
xdim = 4;
ydim = 3;
zdim = 2.7;

Point(1) = {0, 0, 0, cl__1};
Point(2) = {xdim, 0, 0, cl__1};
Point(3) = {xdim, ydim, 0, cl__1};
Point(4) = {0, ydim, 0, cl__1};
Point(5) = {0, 0, zdim, cl__1};
Point(6) = {xdim, 0, zdim, cl__1};
Point(7) = {xdim, ydim, zdim, cl__1};
Point(8) = {0, ydim, zdim, cl__1};
Line(9) = {1, 2};
Line(10) = {2, 3};
Line(11) = {3, 4};
Line(12) = {4, 1};
Line(13) = {5, 6};
Line(14) = {6, 7};
Line(15) = {7, 8};
Line(16) = {8, 5};
Line(17) = {1, 5};
Line(18) = {2, 6};
Line(19) = {3, 7};
Line(20) = {4, 8};
Line Loop(22) = {16, 13, 14, 15};
Plane Surface(22) = {22};
Line Loop(24) = {14, -19, -10, 18};
Plane Surface(24) = {24};
Line Loop(26) = {10, 11, 12, 9};
Plane Surface(26) = {26};
Line Loop(28) = {12, 17, -16, -20};
Plane Surface(28) = {28};
Line Loop(30) = {20, -15, -19, 11};
Plane Surface(30) = {30};
Line Loop(32) = {18, -13, -17, 9};
Plane Surface(32) = {32};
Surface Loop(34) = {22, 28, 26, 24, 30, 32};
Volume(34) = {34};
Physical Surface("Frequency Independent",2) = {22, 24, 26, 28, 30, 32};
Physical Volume(10) = {34};
