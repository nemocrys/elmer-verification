

rhi = 0.4;
rho = 0.5;
    
rii = 0.9;
rio = 1;

cl = 0.049;

Point(0) = {0, 0, 0, cl};

Point(1) = {0, rhi, 0, cl};
Point(2) = {0, rho, 0, cl};
Point(3) = {0, rii, 0, cl};
Point(4) = {0, rio, 0, cl};

Point(5) = {rhi, 0, 0, cl};
Point(6) = {rho, 0, 0, cl};
Point(7) = {rii, 0, 0, cl};
Point(8) = {rio, 0, 0, cl};

Point(9) =  {0, -rhi, 0, cl};
Point(10) = {0, -rho, 0, cl};
Point(11) = {0, -rii, 0, cl};
Point(12) = {0, -rio, 0, cl};

Line(1) = {12,11};
Line(2) = {11,10};
Line(3) = {10,9};
Line(4) = {9,1};
Line(5) = {1,2};
Line(6) = {2,3};
Line(7) = {3,4};

Circle(8) = {1,0,5};
Circle(9) = {5,0,9};
Circle(10) = {2,0,6};
Circle(11) = {6,0,10};
Circle(12) = {3,0,7};
Circle(13) = {7,0,11};
Circle(14) = {4,0,8};
Circle(15) = {8,0,12};

Line Loop(20) = {5,10,11,3,-9,-8};  Plane Surface(21) = {20};
Line Loop(22) = {7,14,15,1,-13,-12};  Plane Surface(23) = {22};

Physical Surface ("volheat") = {21};
Physical Surface ("volins") = {23};

Physical Line ("surfheatin") = {8,9};
Physical Line ("surfheatout") = {10,11};
Physical Line ("surfinsin") = {12,13};
Physical Line ("surfinsout") = {14,15};











