using AMDGPU

# Disable scalar indexing of device arrays: any scalar index into a ROCArray inside a package method
# now throws, turning the no-scalar-indexing invariant from a convention into a checked guarantee.
AMDGPU.allowscalar(false)

const GPU_BACKEND = "ROCm / AMDGPU"
device_functional() = AMDGPU.functional()
to_device(host_array) = AMDGPU.ROCArray(host_array)

include(joinpath(@__DIR__, "..", "suite.jl"))
