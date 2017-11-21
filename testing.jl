# File for testing and interactive development
#######################################
include("Multiprocess.jl")
include("Packages.jl")
include("Globals.jl")
include("Quaternions.jl")
include("Datatypes.jl")
include("Shapes.jl")
include("apufunktiot.jl")
include("Visualization.jl")
include("Body_Init_functions.jl")
include("Moment_of_inertia_functions.jl")
include("Solids_functions.jl")
include("Collision_detection.jl")
include("Contact_modeling.jl")
include("Coordinate_transformations.jl")
include("Dynamics_functions.jl")
include("EOMs.jl")

function test()
    # Rsys = init_RSys([create_cube(ones(3)..., 15.0, 0.6, 0.45) for i in 1:5]);
    Rsys = init_RSys([create_sphere(1.0, 100.0, 0.6, 0.4) for i in 1:2]);
    nb = length(Rsys);
    bodies = Rsys.bodies;
    randomize!(Rsys)
    for i in 1:Rsys.nb
        Rsys.bodies[i].sv.x.z += 4.0;
    end
    # Create plane so body doesn't fall into emptiness
    plane = Plane([0.0,0.0, 1.0], point3D(0.0, 0.0, 0.0))

    window = init_screen();
    origo = origin();
    bodyvis = [body_vis(bodies[i]) for i in 1:nb];

    gl._view(origo, window, camera=:perspective);
    gl._view(bodyvis, window, camera=:perspective);

    sleep(1)
    while Rsys.t < simTime
        # NOTE Plane check is a special case at the moment. Needs to not be.
        list = checkcoll(Rsys, plane);
        Flist = resolvecoll(Rsys, list, plane);

        collList = checkcoll(Rsys);
        FlistColl = resolvecoll(Rsys, collList);

        eval_forces!(Rsys, Flist);
        dynamics!(Rsys)

        body_vis!(bodyvis, Rsys)
        Rsys.t += Δt;
        sleep(0.001);
    end
    return nothing
end



mutable struct rbody{T<:Real}
      x::sa.SVector{3,T} # Reference point location in global coordinates
      ẋ::sa.SVector{3,T}
      a::sa.SVector{3,T}
      F::sa.SVector{3,T}
      T::sa.SVector{3,T}
      ω::sa.SVector{3,T}
      α::sa.SVector{3,T}
      q::sa.SVector{4,T}
      q̇::sa.SVector{4,T}
      R::sa.SArray{Tuple{3,3},T,2,9}
      R_inv::sa.SArray{Tuple{3,3},T,2,9}
      Jb::sa.SArray{Tuple{3,3},T,2,9}
      Jb_inv::sa.SArray{Tuple{3,3},T,2,9}
      J::sa.SArray{Tuple{3,3},T,2,9}
      J_inv::sa.SArray{Tuple{3,3},T,2,9}
      ωs::sa.SArray{Tuple{4,4},T,2,16}
      m_inv::T
end
function translation!(body::rbody{T}) where {T}
    body.a = body.F.*body.m_inv
    body.ẋ += body.a .* Δt
    body.x += body.ẋ .* Δt
    return nothing
end
function skew4(o::sa.SVector{3,T}) where {T}
    @inbounds begin
        o1 = o[1];
        o11 = -o1;
        o2 = o[2];
        o21 = -o2;
        o3 = o[3];
        o31 = -o3;
        T0 = zero(T);
        # os = sa.SArray{Tuple{4,4},T,2,16}(T0, o1, o2, o3, o11, T0, o3, o21, o21, o31, T0, o1, o31, o2, o11, T0)
        os = sa.SMatrix{4,4,T,16}(T0, o1, o2, o3, o11, T0, o3, o21, o21, o31, T0, o1, o31, o2, o11, T0);
    end
    return os
end
function J_local2global(rotmat_body, Jbody)
      @inbounds begin
            if iszero(Jbody) == true
                  J = sa.SMatrix{3,3,Float64,9}(zeros(3,3))
            else
                  RotJ11 = 0.0
                  RotJ21 = 0.0
                  RotJ31 = 0.0
                  RotJ12 = 0.0
                  RotJ22 = 0.0
                  RotJ32 = 0.0
                  RotJ13 = 0.0
                  RotJ23 = 0.0
                  RotJ33 = 0.0
                  for i in 1:3
                        RotJ11 += rotmat_body[1,i]*Jbody[i,1]
                        RotJ21 += rotmat_body[2,i]*Jbody[i,1]
                        RotJ31 += rotmat_body[3,i]*Jbody[i,1]
                        RotJ12 += rotmat_body[1,i]*Jbody[i,2]
                        RotJ22 += rotmat_body[2,i]*Jbody[i,2]
                        RotJ32 += rotmat_body[3,i]*Jbody[i,2]
                        RotJ13 += rotmat_body[1,i]*Jbody[i,3]
                        RotJ23 += rotmat_body[2,i]*Jbody[i,3]
                        RotJ33 += rotmat_body[3,i]*Jbody[i,3]
                  end
                  J11 = RotJ11*rotmat_body[1,1] + RotJ12*rotmat_body[1,2] + RotJ13*rotmat_body[1,3]
                  J21 = RotJ21*rotmat_body[1,1] + RotJ22*rotmat_body[1,2] + RotJ23*rotmat_body[1,3]
                  J31 = RotJ31*rotmat_body[1,1] + RotJ32*rotmat_body[1,2] + RotJ33*rotmat_body[1,3]

                  J12 = RotJ11*rotmat_body[2,1] + RotJ12*rotmat_body[2,2] + RotJ13*rotmat_body[2,3]
                  J22 = RotJ21*rotmat_body[2,1] + RotJ22*rotmat_body[2,2] + RotJ23*rotmat_body[2,3]
                  J32 = RotJ31*rotmat_body[2,1] + RotJ32*rotmat_body[2,2] + RotJ33*rotmat_body[2,3]

                  J13 = RotJ11*rotmat_body[3,1] + RotJ12*rotmat_body[3,2] + RotJ13*rotmat_body[3,3]
                  J23 = RotJ21*rotmat_body[3,1] + RotJ22*rotmat_body[3,2] + RotJ23*rotmat_body[3,3]
                  J33 = RotJ31*rotmat_body[3,1] + RotJ32*rotmat_body[3,2] + RotJ33*rotmat_body[3,3]
                  J = sa.SMatrix{3,3,Float64,9}(J11, J21, J31, J12, J22, J32, J13, J23, J33)
            end
      end
      return J
end
function rotmat(q::sa.SVector{4,T}) where {T<:Real}
      @inbounds begin
            e1_toiseen = q[1]*q[1]
            e2_toiseen = q[2]*q[2]
            e3_toiseen = q[3]*q[3]
            e4_toiseen = q[4]*q[4]
            e2e3 = q[2]*q[3]
            e1e4 = q[1]*q[4]
            e2e4 = q[2]*q[4]
            e1e3 = q[1]*q[3]
            e3e4 = q[3]*q[4]
            e1e2 = q[1]*q[2]
            # Wikipedia
            R11=e1_toiseen + e2_toiseen - e3_toiseen - e4_toiseen
            R21=2.0*(e2e3 + e1e4)
            R31=2.0*(e2e4 - e1e3)

            R12=2.0*(e2e3 - e1e4)
            R22=e1_toiseen - e2_toiseen + e3_toiseen - e4_toiseen
            R32=2.0*(e3e4 + e1e2)

            R13=2.0*(e2e4 + e1e3)
            R23=2.0*(e3e4 - e1e2)
            R33=e1_toiseen - e2_toiseen - e3_toiseen + e4_toiseen

            R = sa.SArray{Tuple{3,3},T,2,9}(R11, R21, R31, R12, R22, R32, R13, R23, R33)
      end
      return R
end
function rotation!(body::rbody{T}) where {T}
    body.α = body.J_inv*(body.T - cross(body.ω, (body.J*body.ω)))
    body.ω += body.α .* Δt
    # Quaternions time derivative
    body.q̇ = 0.5.*cross(body.ω,body.q)
    # Update & normalize Euler parameters
    body.q += body.q̇ .* Δt
    body.q = normalize(body.q);
    # Update Rotation matrices
    body.R = rotmat(body.q)
    body.R_inv = inv(body.R)
    # Update mass moment of inertias in global frame
    body.J = J_local2global(body.R, body.Jb)
    body.J_inv = J_local2global(body.R, body.Jb_inv)
    return nothing
end

function testia()
    stp = 1000000;

    voima = rand(3);
    vääntö = rand(3);
    a = create_sphere(1.0, 100.0, 0.6, 0.4);
    a.f.F[:] = voima;
    a.f.T[:] = vääntö;

    b = rbody([sa.@SVector(zeros(3)) for i in 1:7]..., [sa.@SVector(zeros(4)) for i in 1:2]..., [sa.@SMatrix(eye(3)) for i in 1:6]..., sa.@SMatrix(eye(4)), inv(100.0))
    b.q = [1.0, zeros(3)...];
    b.J = eye(3)*40.0;
    b.Jb = eye(3)*40.0;
    b.J_inv = inv(b.J)
    b.Jb_inv = inv(b.Jb)
    b.F = voima;
    b.T = vääntö;

    @time for i in 1:stp
        translation!(a)
        rotation!(a)
    end
    @time for i in 1:stp
        translation!(b)
        rotation!(b)
    end
    return nothing
end


import Base.cross
function cross(a::sa.SVector{3,T}, b::sa.SVector{4,T}) where {T}
    @inbounds begin
        a1 = a[1]
        a2 = a[2]
        a3 = a[3]
        b1 = b[1]
        b2 = b[2]
        b3 = b[3]
        b4 = b[4]
    end
    c1 = -a1*b2 - a2*b3 - a3*b4
    c2 = a1*b1 - a3*b3 + a2*b4
    c3 = a2*b1 + a3*b2 - a1*b4
    c4 = a3*b1 - a2*b2 + a1*b3
    c = sa.SVector{4,T}(c1, c2, c3, c4)
end

function inverse!(minv::sa.SArray{Tuple{3,3},T,2,9}, m::sa.SArray{Tuple{3,3},T,2,9}) where {T}
    minv = inv(m);
    return minv
end
