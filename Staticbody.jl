mutable struct Static{T<:Real}
    # Body shape described as a HomogenousMesh in local coordinates
    mesh::gt.HomogenousMesh{gt.Point{3,Float32},gt.Face{3,gt.OffsetInteger{-1,UInt32}},gt.Normal{3,Float32},Void,Void,Void,Void}
    # Coll shape
    coll::Collshape{T}
end    
