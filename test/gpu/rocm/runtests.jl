using AMDGPU

# Disable scalar indexing of device arrays: any scalar index into a ROCArray inside a package method
# now throws, turning the no-scalar-indexing invariant from a convention into a checked guarantee.
AMDGPU.allowscalar(false)

const GPU_BACKEND = "ROCm / AMDGPU"
device_functional() = AMDGPU.functional()
to_device(host_array) = AMDGPU.ROCArray(host_array)
const DEVICE_FAMILY = AMDGPU.ROCArray
DEVICE_CODE_LLVM(io, f) = AMDGPU.@device_code_llvm io=io f()   # the device LLVM of every kernel f launches

include(joinpath(@__DIR__, "..", "suite.jl"))
