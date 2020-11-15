// Gmsh project created on Mon Feb 03 10:13:39 2020
SetFactory("OpenCASCADE");
//+
Cylinder(1) = {0, 0, 0, 0, 0, 1, 1, 2*Pi};
//+
Physical Surface("PR", 1) = {2, 3, 1};
//+
Physical Volume("Volume", 10) = {1};

Mesh.CharacteristicLengthFactor = 1;
Mesh.MshFileVersion = 2.2;

Mesh 3;
OptimizeMesh "Netgen";
SetOrder 8;
// Save "3DCylinderN8ESF013.msh"

//For i In {1:3}
	//If(i==1)
		//Save Sprintf("ConvergenceTestLibP2\\3DCylinderN%g_Split%03g.msh",1,0);
		//For N In {2:8}
		//	SetOrder N;
		//	Save Sprintf("ConvergenceTestLibP2\\3DCylinderN%g_Split%03g.msh",N,0);
		//EndFor
		//SetOrder 1;
	//EndIf
	//RefineMesh;
	//Save Sprintf("ConvergenceTestLibP2\\3DCylinderN%g_Split%03g.msh",1,i);
	//For N In {2:8}
	//	SetOrder N;
	//	Save Sprintf("ConvergenceTestLibP2\\3DCylinderN%g_Split%03g.msh",N,i);
	//EndFor
	//SetOrder 1;
//EndFor


//RefineMesh;
//RefineMesh;
//RefineMesh;
//For N In {4:8}
//	SetOrder N;
//	Save Sprintf("ConvergenceTestLibP2\\3DCylinderN%g_Split%03g.msh",N,3);
//EndFor