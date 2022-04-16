module AcousticRayTracers

using UnderwaterAcoustics
using UnderwaterAcoustics: check, envrealtype, fast_absorption, RayArrival

export RaySolver

include("RaySolver.jl")

function __init__()
  UnderwaterAcoustics.addmodel!(RaySolver)
end

end # module
