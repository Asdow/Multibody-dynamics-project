# Functions related to quaternions / Euler parameters 
#######################################
# Calculates the rotation matrix based on quaternions
function func_rotation_matrix!(Rotation_matrix, euler)
    e1_toiseen = euler[1]^2
    e2_toiseen = euler[2]^2
    e3_toiseen = euler[3]^2
    e4_toiseen = euler[4]^2
    e2e3 = euler[2]*euler[3]
    e1e4 = euler[1]*euler[4]
    e2e4 = euler[2]*euler[4]
    e1e3 = euler[1]*euler[3]
    e3e4 = euler[3]*euler[4]
    e1e2 = euler[1]*euler[2]
    # Wikipedia
    Rotation_matrix[1,1]=e1_toiseen + e2_toiseen - e3_toiseen - e4_toiseen
    Rotation_matrix[2,1]=2.0*(e2e3 + e1e4)
    Rotation_matrix[3,1]=2.0*(e2e4 - e1e3)

    Rotation_matrix[1,2]=2.0*(e2e3 - e1e4)
    Rotation_matrix[2,2]=e1_toiseen - e2_toiseen + e3_toiseen - e4_toiseen
    Rotation_matrix[3,2]=2.0*(e3e4 + e1e2)

    Rotation_matrix[1,3]=2.0*(e2e4 + e1e3)
    Rotation_matrix[2,3]=2.0*(e3e4 - e1e2)
    Rotation_matrix[3,3]=e1_toiseen - e2_toiseen - e3_toiseen + e4_toiseen

    return nothing
end
