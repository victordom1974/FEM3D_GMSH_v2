// Janus hexahedron particles: 11 two-material cubes floating in a box.
// Original: hexahedron.geo (Nov 2022, Alejandro & vD).
// 2026-10-12: meshing loop disabled (MesshingQ = 0), optimization enabled;
//             the refinement hierarchy is built with refineTtrh, not here.
// Generate the coarse mesh with:
//   gmsh -3 -order 2 -format msh41 meshes/janus.geo -o meshes/janus_p2.msh

// November 2022
//
// Adapted to my own 3D FEM code (Alejandro & vD)




//+
SetFactory("OpenCASCADE");


// My first Macro... Yahoo! 
Macro Particle
 
Point(j+1) = { l, l ,-l,h};
Point(j+2) = {-l, l ,-l,h};
Point(j+3) = {-l,-l ,-l,h};
Point(j+4) = { l,-l ,-l,h};

Point(j+5) = { l, l ,0,h};
Point(j+6) = {-l, l ,0,h};
Point(j+7) = {-l,-l ,0,h};
Point(j+8) = { l,-l ,0,h};

Point(j+9) = { l, l  , l,h};
Point(j+10) = {-l, l , l,h};
Point(j+11) = {-l,-l , l,h};
Point(j+12) = { l,-l , l,h};
 
 
 //XY-Edges 
 

//+
Line(j+1) = {j+1,j+2};
Line(j+2) = {j+2,j+3};
Line(j+3) = {j+3,j+4};
Line(j+4) = {j+4,j+1};

//+
Line(j+5) = {j+5,j+6};
Line(j+6) = {j+6,j+7};
Line(j+7) = {j+7,j+8};
Line(j+8) = {j+8,j+5};


//+
Line(j+9)  = {j+9 ,j+10};
Line(j+10) = {j+10,j+11};
Line(j+11) = {j+11,j+12};
Line(j+12) = {j+12,j+9};


 
 //Z-Edges 
 

//+
Line(j+13) = {j+1,j+5};
Line(j+14) = {j+2,j+6};
Line(j+15) = {j+3,j+7};
Line(j+16) = {j+4,j+8};

//+
Line(j+17) = {j+5,j+9};
Line(j+18) = {j+6,j+10};
Line(j+19) = {j+7,j+11};
Line(j+20) = {j+8,j+12};


// Faces
//
// XY 
Curve Loop(j+1) = {j+1, j+2, j+3, j+4};
Plane Surface(j+1) = {j+1};
Curve Loop(j+2) = {j+5, j+6, j+7, j+8};
Plane Surface(j+2) = {j+2};
Curve Loop(j+3) = {j+9, j+10, j+11, j+12};
Plane Surface(j+3) = {j+3};

// Z
//+
Curve Loop(j+4) = {j+1, j+14, -(j+5), -(j+13)};
Plane Surface(j+4) = {j+4};
//+
Curve Loop(j+5) = {j+2, j+15, -(j+6), -(j+14)};
Plane Surface(j+5) = {j+5};
//+
Curve Loop(j+6) = {j+3, j+16, -(j+7), -(j+15)};
Plane Surface(j+6) = {j+6};
//+
Curve Loop(j+7) = {j+4, j+13, -(j+8), -(j+16)};
Plane Surface(j+7) = {j+7};
//
//+
Curve Loop(j+8) = {j+5, j+18, -(j+9), -(j+17)};
Plane Surface(j+8) = {j+8};
//+
Curve Loop(j+9) = {j+6, j+19, -(j+10), -(j+18)};
Plane Surface(j+9) = {j+9};
//+
Curve Loop(j+10) = {j+7, j+20, -(j+11), -(j+19)};
Plane Surface(j+10) = {j+10};

//+
Curve Loop(j+11) = {j+8, j+17, -(j+12), -(j+20)};
Plane Surface(j+11) = {j+11};
//



//+
//Rotate {{d1,d2,d3},{0,0,0}, th} {
//  Volume{j+1}; Volume{j+2}; 
//  Surface{j+1,j+2, j+4, j+5,j+6,j+7};
//  Surface{j+2,j+3, j+8, j+9,j+10,j+11};
//}


// When rotating and translating the surfaces, the tags 
// for lines and points are changed!!!
//
// So, if we construct first the volumes and next rotate them
// we'd loss the control on the surfaces of these volumes
//
// The solution is clear: construct the surfaces, rotate and translate them
// and next construct the volumes. Doing so, the labels for the surfaces
// j+1:j+11 remain unchanged. 
//
// 

Rotate {{d1,d2,d3},{0,0,0}, th} {
  Surface{j+1:j+11}; 
}
Translate {tx, ty, tz} {
Surface{j+1:j+11}; 
}

//+
Surface Loop(j+21) = {j+1,j+2, j+4, j+5,j+6,j+7};
//+
Volume(j+1) = {j+21};
//+
Surface Loop(j+22) = {j+2,j+3, j+8, j+9,j+10,j+11};
//+
Volume(j+2) = {j+22};



Physical Volume(Str(name1),j+1) = {j+1};
Physical Volume(Str(name2),j+2) = {j+2};
Physical Surface(Str(nameSurface1),j+50) = {j+1,j+2, j+4, j+5,j+6,j+7};
Physical Surface(Str(nameSurface2),j+51) = {j+2,j+3, j+8, j+9,j+10,j+11};

Return
 
 
 // LET'S START
 
 
j = 100; // For points indices
h  = 0.3; 
l = 0.3;
th = Pi/3;  
name1 = "Janus 1a"; 
name2 = "Janus 1b";
nameSurface1 = "Janus 1a - Surface";
nameSurface2 = "Janus 1b - Surface";

d1 = 1; d2 = 1; d3 = 1;
tx = 0; ty = 0; tz = 1;
Call Particle;

j = 200; // For points indices
h  = 0.3; 
l = 0.3;
th = 0;  
d1 = 1; d2 = 1; d3 = 1;
tx = 0; ty = -1; tz = 0;
name1 = "Janus 2a"; 
name2 = "Janus 2b";
nameSurface1 = "Janus 2a - Surface";
nameSurface2 = "Janus 2b - Surface";
Call Particle;

j = 300; // For points indices
h  = 0.3; 
l = 0.3;
th = Pi/2;  
d1 = 1; d2 = 1; d3 = 1;
tx = -1; ty = 0; tz = 0;
name1 = "Janus 3a"; 
name2 = "Janus 3b";
nameSurface1 = "Janus 3a - Surface";
nameSurface2 = "Janus 3b - Surface";
Call Particle;

j = 400; // For points indices
h  = 0.3; 
l = 0.3;
th = Pi/6;  
d1 = 1; d2 = 1; d3 = 1;
tx = 1; ty = 1; tz = 1;
name1 = "Janus 4a" ;
name2 = "Janus 4b";
nameSurface1 = "Janus 4a - Surface";
nameSurface2 = "Janus 4b - Surface";
Call Particle;

j = 500; // For points indices
h  = 0.3; 
l = 0.3;
th = Pi/6;  
d1 = 1; d2 = 1; d3 = 1;
tx = -1; ty = -1; tz = -1;
name1 = "Janus 5a" ;
name2 = "Janus 5b";
nameSurface1 = "Janus 5a - Surface";
nameSurface2 = "Janus 5b - Surface";
Call Particle;


j = 600; // For points indices
h  = 0.3; 
l = 0.3;
th = Pi/4;  
d1 = 1; d2 = 1; d3 = 1;
tx = 1; ty = -1; tz = 1;
name1 = "Janus 6a"; 
name2 = "Janus 6b";
nameSurface1 = "Janus 6a - Surface";
nameSurface2 = "Janus 6b - Surface";
Call Particle;

j = 700; // For points indices
h  = 0.3; 
l = 0.3;
th = Pi/4;  
d1 = 1; d2 = 1; d3 = 1;
tx = 1; ty = -1; tz = -1;
name1 = "Janus 7a" ;
name2 = "Janus 7b";
nameSurface1 = "Janus 7a - Surface";
nameSurface2 = "Janus 7b - Surface";
Call Particle;

j = 800; // For points indices
h  = 0.3; 
l = 0.3;
th = Pi/3;  
d1 = 0; d2 = 1; d3 = 1;
tx = -0.6; ty = 1; tz = -1;
name1 = "Janus 8a" ;
name2 = "Janus 8b";
nameSurface1 = "Janus 8a - Surface";
nameSurface2 = "Janus 8b - Surface";
Call Particle;

j = 900; // For points indices
h  = 0.3; 
l = 0.3;
th = Pi/3;  
d1 = 1; d2 = 1; d3 = 1;
tx = -1; ty = 0.7; tz = 1;
name1 = "Janus 9a"; 
name2 = "Janus 9b";
nameSurface1 = "Janus 9a - Surface";
nameSurface2 = "Janus 9b - Surface";
Call Particle;

j = 1000; // For points indices
h  = 0.3; 
l = 0.3;
th = Pi/3;  
d1 = 1; d2 = 1; d3 = 1;
tx =  1; ty = 1; tz = -0.9;
name1 = "Janus 10a" ;
name2 = "Janus 10b";
nameSurface1 = "Janus 10a - Surface";
nameSurface2 = "Janus 10b - Surface";
Call Particle;

j = 1100; // For points indices
h  = 0.3; 
l = 0.3;
th = -Pi/3;  
d1 = 1; d2 = 1; d3 = 1;
tx =  -1; ty = -1; tz = 0.7;
name1 = "Janus 11a" ;
name2 = "Janus 11b";
nameSurface1 = "Janus 11a - Surface";
nameSurface2 = "Janus 11b - Surface";
Call Particle;



h  = 0.5;
lx = 4; 
ly = 4;
lz = 3; 
Point(1) = { lx, ly ,-lz,h};
Point(2) = {-lx, ly ,-lz,h};
Point(3) = {-lx,-ly ,-lz,h};
Point(4) = { lx,-ly ,-lz,h};


Point(5) = { lx, ly ,lz,h};
Point(6) = {-lx, ly ,lz,h};
Point(7) = {-lx,-ly ,lz,h};
Point(8) = { lx,-ly ,lz,h};


//+
Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};
//+
Line(5) = {5,6};
Line(6) = {6,7};
Line(7) = {7,8};
Line(8) = {8,5};
//
Line(9) =  {1,5};
Line(10) = {2,6};
Line(11) = {3,7};
Line(12) = {4,8};



Curve Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};
Curve Loop(2) = {5, 6, 7, 8}; 
Plane Surface(2) = {2};
Curve Loop(3) = {1,10,-5,-9};
Plane Surface(3) = {3};
Curve Loop(4) = {2,11,-6,-10};
Plane Surface(4) = {4};
Curve Loop(5) = {3,12,-7,-11};
Plane Surface(5) = {5};
Curve Loop(6) = {4,9,-8,-12};
Plane Surface(6) = {6};


//+
Surface Loop(1) = {1:6};
Volume(3) = {1};
//+

//BooleanDifference{ Volume{1}; Delete; }{ Volume{701}; Volume{702}; Volume{1002}; Volume{1001}; Volume{202}; Volume{201}; Volume{502}; Volume{501}; Volume{801}; Volume{802}; Volume{301}; Volume{302}; Volume{1101}; Volume{1102}; Volume{602}; Volume{601}; Volume{402}; Volume{401}; Volume{102}; Volume{101}; Volume{902}; Volume{901};}

p = newv;
BooleanDifference(1) = { Volume{3};  }{ Volume{101:1101:100};Volume{102:1102:100};};


Physical Surface("OuterBoundary",100) = {1:6};
  


//+
Recursive Delete {
  Volume{3}; 
}

Physical Volume("OuterVolume",1) = {1};

//Physical Surface(Str(nameSurface2),j+51) = {j+2,j+3, j+8, j+9,j+10,j+11};


//Mesh.MeshSizeMin =.3;

// Used up to Nov 15 2022
// Mesh.MeshSizeMax =.6;
// Used up to Nov 22 2022
 Mesh.MeshSizeMax =1;
Mesh.Optimize = 1;
Mesh.OptimizeNetgen = 1;
Mesh.OptimizeThreshold = 0.3;
 
Mesh.MshFileVersion = 4.1; 
Mesh.ElementOrder = 1;

// Let's talk about 3D meshes.... 
MesshingQ = 0;  //Meshing "a saco paco"?  If so, go ahead. Be cautious with values greater than 3 

// Meshing can be done from the gui interface... 


If (MesshingQ  > 0)
	Mesh 2; // Just for seeing the (2D) mesh on the boundaries. 

 	Mesh 3; 
 	Printf("First Mesh");
	Save "JanusHexahedronP1Ref01.msh"; 
	
 	Printf("First Mesh - Order 2");
	Mesh.ElementOrder = 2;
	Mesh 3; 
	Save "JanusHexahedronP2Ref01.msh";
	 
 	Printf("First Mesh - Order 3");
	Mesh.ElementOrder = 3;
	Mesh 3; 
	Save "JanusHexahedronP3Ref01.msh"; 
	
	 
 	Printf("First Mesh - Order 4");
	Mesh.ElementOrder = 4;
	Mesh 3; 
	Save "JanusHexahedronP4Ref01.msh"; 
EndIf

If (MesshingQ  > 1)
	Printf("Second Mesh");
	Mesh.ElementOrder = 1;
	RefineMesh;  
		 
 	Printf("First Mesh");
	Save "JanusHexahedronP1Ref02.msh"; 
	
 	Printf("First Mesh - Order 2");
	Mesh.ElementOrder = 2;
	Mesh 3; 
	Save "JanusHexahedronP2Ref02.msh";
	 
 	Printf("First Mesh - Order 3");
	Mesh.ElementOrder = 3;
	Mesh 3; 
	Save "JanusHexahedronP3Ref02.msh"; 
	
	 
 	Printf("First Mesh - Order 4");
	Mesh.ElementOrder = 4;
	Mesh 3; 
	Save "JanusHexahedronP4Ref02.msh"; 
EndIf 

If (MesshingQ  > 2)
	Printf("Third Mesh");
	Mesh.ElementOrder = 1;
	RefineMesh;  
	
	 Printf("First Mesh");
	Save "JanusHexahedronP1Ref03.msh"; 
	
 	Printf("First Mesh - Order 2");
	Mesh.ElementOrder = 2;
	Mesh 3; 
	Save "JanusHexahedronP2Ref03.msh";
	 
 	Printf("First Mesh - Order 3");
	Mesh.ElementOrder = 3;
	Mesh 3; 
	Save "JanusHexahedronP3Ref03.msh"; 
	
	 
 	Printf("First Mesh - Order 4");
	Mesh.ElementOrder = 4;
	Mesh 3; 
	Save "JanusHexahedronP4Ref03.msh"; 
EndIf

If (MesshingQ  > 3)
	Printf("Fourth Mesh");
	Mesh.ElementOrder = 1;
	RefineMesh;  
	Printf("First Mesh");
	Save "JanusHexahedronP1Ref04.msh"; 
	
	Printf("First Mesh - Order 2");
	Mesh.ElementOrder = 2;
	Mesh 3; 
	Save "JanusHexahedronP2Ref04.msh";
 
EndIf


If (MesshingQ  > 4)
	Printf("Fifth Mesh");
	Mesh.ElementOrder = 1;
	RefineMesh;  
	Save "JanusHexahedron05.msh"; 
EndIf

//+
//+
//+
