using TestItems

@testitem "aqua" begin
  using AcousticRayTracers
  using Aqua
  Aqua.test_all(AcousticRayTracers; persistent_tasks=false)
end
