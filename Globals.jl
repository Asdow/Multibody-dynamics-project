# Global variables
# gravity
global const g = sa.SVector{3,Float64}(0.0, 0.0, -9.81);
# Timestep
global const Δt = 0.01;
# End of simulation time
global const simTime = 15.0;



# Global temp. arrays etc.
global const vec3 = [sa.MVector{3,Float64}(0.0,0.0,0.0) for i in 1:3];
global const gmat3 = [sa.MMatrix{3,3, Float64, 9}(zeros(3,3)) for i in 1:3];
global const gmat4 = [sa.MMatrix{4,4, Float64, 16}(zeros(4,4)) for i in 1:3];
global const gmat4f0 = [sa.MMatrix{4,4, Float32, 16}(zeros(Float32,4,4)) for i in 1:3];
