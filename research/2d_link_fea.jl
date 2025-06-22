using Pkg; Pkg.activate(@__DIR__)
using Makie
using Gmsh

using Gmsh
using GridapGmsh
# TODO:
# - add center point to the 

link_lengh = 1.0
link_width = 0.2

mm = 1e-3

gmsh.initialize()

gmsh.option.setNumber("General.Terminal", 0);
ms_max = 1;
ms_min = 0.1;
MSFC = 6;
EL_ORDER = 1
# gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", MSFC);
# gmsh.option.setNumber("Mesh.MeshSizeMax", diam/2);
# gmsh.option.setNumber("Mesh.MeshSizeMin", diam/5);
gmsh.option.setNumber("Mesh.MaxNumThreads3D", 8);
gmsh.option.setNumber("Mesh.ElementOrder", EL_ORDER);
# gmsh.option.setNumber("Mesh.SubdivisionAlgorithm", 2);
gmsh.option.setNumber("General.Verbosity", 5)
gmsh.model.add("L_2_1")

Lc1 = 20mm;



factory = gmsh.model.geo;



lhc_p = factory.add_point(0, 0, 0, Lc1); #left hinge center
lb_p = factory.add_point(0, -link_width / 2, 0, Lc1);
lm_p = factory.add_point(-link_width / 2, 0, 0, Lc1);
lt_p = factory.add_point(0, link_width / 2, 0, Lc1);
rhc_p = factory.add_point(link_lengh, 0, 0, Lc1);
rb_p = factory.add_point(link_lengh, -link_width / 2, 0, Lc1);
rm_p = factory.add_point(link_lengh + link_width / 2, 0, 0, Lc1);
rt_p = factory.add_point(link_lengh, link_width / 2, 0, Lc1);


lb_arc = factory.addCircleArc(lb_p, lhc_p, lm_p);
lt_arc = factory.addCircleArc(lm_p, lhc_p, lt_p);
t_line = factory.addLine(lt_p, rt_p);
rt_arc = factory.addCircleArc(rt_p, rhc_p, rm_p);
rb_arc = factory.addCircleArc(rm_p, rhc_p, rb_p);
b_line = factory.addLine(rb_p, lb_p);

# factory.addPhysicalGroup(1, [12], -1, "left");
# factory.addPhysicalGroup(1, [6], -1, "right");
# factory.addPhysicalGroup(1, [1, 2, 3, 4, 5, 7, 8, 9, 10, 11], -1, "free");

# factory.synchronize();

outter_bound = factory.addCurveLoop([lb_arc, lt_arc, t_line, rt_arc, rb_arc, b_line]);

factory.addPlaneSurface([outter_bound], 6)


factory.synchronize()

# factory.addPhysicalGroup(2, [6], -1, "test_sample");

# factory.synchronize()

gmsh.model.mesh.generate(2)


name = (@__DIR__) * "\\" * "test_sample.msh";
print(name)
gmsh.write(name)

gmsh.fltk.run()

gmsh.finalize()
