# Nearest plans

1. 2D problems
   1. Implement hinge joint
   2. Change autodiff jacobian into analytical
   3. Implement Slider Joint
   4. Make Module and split into scripts
   5. Add tests
   6. Visualization (GLMakie)
2. 3D problems
   1. 3D body
   2. Hinge, Sдider, Spherical joints
   3. Implement Springs and Damphers
   4. Visualization (GLMakie)
3. HIL - look for the RT hardware to run
   1. RT - Linux (LinuxCNC) 
      1. Xenomai
      2. RTLinux
   2. SoC? 
4. Flexible Bodies
   1. Full mesh based on FFRF
   2. Model Reduction
   3. Hyperreduction
5. ???
6. Perfomance:
   1. Add type speciication for `rhs` and `state` parameters of functions `add_joint_to_rhs!` and `add_joint_to_rhs!`