# Functions for the rigid body dynamics
# Funktio joka laskee omegalle 4x4 kokoisen skew symmetric matriisin globaalissa referenssikoordinaatistossa
function func_omega_tilde4_world!(omega_tilde4, omega)
    @inbounds begin
    omega_tilde4[1,1] = 0.0
    omega_tilde4[2,1] = omega[1]
    omega_tilde4[3,1] = omega[2]
    omega_tilde4[4,1] = omega[3]

    omega_tilde4[1,2] = -omega[1]
    omega_tilde4[2,2] = 0.0
    omega_tilde4[3,2] = omega[3]
    omega_tilde4[4,2] = -omega[2]

    omega_tilde4[1,3] = -omega[2]
    omega_tilde4[2,3] = -omega[3]
    omega_tilde4[3,3] = 0.0
    omega_tilde4[4,3] = omega[1]

    omega_tilde4[1,4] = -omega[3]
    omega_tilde4[2,4] = omega[2]
    omega_tilde4[3,4] = -omega[1]
    omega_tilde4[4,4] = 0.0
    end
    return nothing
end
# velocity for Quaternions
function func_euler_vel!(euler_vel, euler, omega_tilde4)
    @simd for i in range(1,4)
        @inbounds euler_vel[i] = 0.5*(omega_tilde4[i,1]*euler[1] + omega_tilde4[i,2]*euler[2] + omega_tilde4[i,3]*euler[3] + omega_tilde4[i,4]*euler[4])
    end
    return nothing
end
# Quaternions
function func_euler!(euler_acc, euler_vel, euler, delta_t)
    for i in range(1,4)
        euler[i] = euler[i] + euler_vel[i]*delta_t + 0.5*euler_acc[i]*(delta_t^2)
    end
    return nothing
end
# Function that clears a body's forces
function func_clear_forces!(kpl::Kappale)
      fill!(kpl.f.force_body, 0.0)
      fill!(kpl.f.force_world, 0.0)
      fill!(kpl.f.torque_body, 0.0)
      fill!(kpl.f.torque_world, 0.0)
      return nothing
end
