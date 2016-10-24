# This file contains functions to assemble the system matrices and EOMs for rigid bodies.
############################################################
# Build inverted system Mass Matrix
function func_build_sysM_inv(Nbodies)
      n = length(Nbodies)
      SysM_inv = zeros(Float64, 6*n, 6*n)

      # Masses
      for ii = 1:n
            for i = 1:3
                  SysM_inv[i+6*(ii-1), i+6*(ii-1)] = Nbodies[ii].md.massa_inv
            end
      end
      # Inertias
      for ii = 1:n
            for j = 1:3
                  for i = 1:3
                        SysM_inv[i+3+6*(ii-1), j+3+6*(ii-1)] = Nbodies[ii].J.world_inv[i,j]
                  end
            end
      end
return SysM_inv
end
function func_build_sysM_inv!(Nbodies, SysM_inv)
      n = length(Nbodies)
      # Masses
      for ii = 1:n
            for i = 1:3
                  SysM_inv[i+6*(ii-1), i+6*(ii-1)] = Nbodies[ii].md.massa_inv
            end
      end
      # Inertias
      for ii = 1:n
            for j = 1:3
                  for i = 1:3
                        SysM_inv[i+3+6*(ii-1), j+3+6*(ii-1)] = Nbodies[ii].J.world_inv[i,j]
                  end
            end
      end
return nothing
end
###################
# Build system vector of forces
function func_build_sysF(Nbodies)
      n = length(Nbodies)
      #sysF = Array(6*n)
      sysF = Float64[]
      for i = 1:n
            push!(sysF, Nbodies[i].f.force_world...)
            push!(sysF, Nbodies[i].f.torque_world...)
      end
      return sysF
end
function func_build_sysF!(Nbodies, sysF)
      empty!(sysF)
      n = length(Nbodies)
      for i = 1:n
            push!(sysF, Nbodies[i].f.force_world...)
            push!(sysF, Nbodies[i].f.torque_world...)
      end
      return nothing
end
function func_build_sysF2(Nbodies)
      n = length(Nbodies)
      sysF = Array(Float64, 6*n)
      for ii = 1:n
            for i = 1:3
                  sysF[i+6*(ii-1)] = Nbodies[ii].f.force_world[i]
                  sysF[i+3+6*(ii-1)] = Nbodies[ii].f.torque_world[i]
            end
      end
      return sysF
end
function func_build_sysF2!(Nbodies, sysF)
      n = length(Nbodies)
      for ii = 1:n
            for i = 1:3
                  sysF[i+6*(ii-1)] = Nbodies[ii].f.force_world[i]
                  sysF[i+3+6*(ii-1)] = Nbodies[ii].f.torque_world[i]
            end
      end
      return nothing
end

function timeit(Nbodies, sysF)
      @time for i = 1:10^3
            func_build_sysF!(Nbodies, sysF)
      end
      @time for i = 1:10^3
            func_build_sysF2!(Nbodies, sysF)
      end
      return nothing
end
###################
# Build system state vector
function func_build_sysSv(Nbodies)
      n = length(Nbodies)
      sysSv = Array(Float64, 6*n)
      for ii = 1:n
            for i = 1:3
                  sysSv[i+6*(ii-1)] = Nbodies[ii].sv.P[i]
                  sysSv[i+3+6*(ii-1)] = Nbodies[ii].sv.L[i]
            end
      end
      return sysSv
end
function func_build_sysSV!(Nbodies, sysSv)
      n = length(Nbodies)
      for ii = 1:n
            for i = 1:3
                  sysSv[i+6*(ii-1)] = Nbodies[ii].sv.P[i]
                  sysSv[i+3+6*(ii-1)] = Nbodies[ii].sv.L[i]
            end
      end
      return nothing
end
