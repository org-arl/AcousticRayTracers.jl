module AcousticRayTracers

using UnderwaterAcoustics
using UnderwaterAcoustics: AbstractRayPropagationModel, AbstractModePropagationModel
using UnderwaterAcoustics: AbstractAcousticSource, AbstractAcousticReceiver, AcousticReceiverGrid2D
using UnderwaterAcoustics: RayArrival, SampledFieldX, SampledFieldZ, xyz, tmap, env_type
using UnderwaterAcoustics: is_range_dependent, is_constant, value, in_units, db2amp

export RaySolver, BellhopJL

include("RaySolver.jl")

# BellhopJL lives in an internal submodule; note that unlike the rest of this
# package (MIT), src/bellhopjl/** is GPL-3.0-or-later (see LICENSE-GPL-3.0)
include("bellhopjl/BellhopJLCore.jl")
using .BellhopJLCore: BellhopJL

end # module
