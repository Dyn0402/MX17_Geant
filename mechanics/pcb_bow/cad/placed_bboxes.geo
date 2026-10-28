SetFactory("OpenCASCADE");
Geometry.OCCImportLabels = 1;
v() = ShapeFromFile("/home/dylan/x17/mx17_cad/MX17 2026 v1.3 05 06 2026.stp");
For i In {0:#v()-1}
  bb() = BoundingBox Volume{v(i)};
  Printf("VOL %g %g %g %g %g %g %g", v(i), bb(0), bb(1), bb(2), bb(3), bb(4), bb(5));
EndFor
