Merge "stator.step";
Merge "rotor_animation.msh";

// Enable OpenGL transparency
General.AlphaBlending = 1;

// Stator STEP: show as transparent surface
Geometry.Surfaces = 1;
Geometry.SurfaceType = 2;
Geometry.Light = 1;

// Per-surface color with alpha (RGBA, alpha 0-255)
// Apply to all surfaces of the STEP model
Color {200, 200, 220, 50} { Surface{:}; }

// Rotor mesh
Mesh.SurfaceFaces = 1;
Mesh.VolumeEdges = 0;
Mesh.SurfaceEdges = 0;
Mesh.NumSubEdges = 4;
Mesh.ColorCarousel = 0;
Mesh.Color.Triangles = {0, 200, 50, 255};
Mesh.Color.Tetrahedra = {0, 200, 50, 255};

// Displacement animation
View[0].VectorType = 5;
View[0].DisplacementFactor = 1.0;
View[0].ShowElement = 1;
View[0].ShowScale = 0;
