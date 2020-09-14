PROJECT(SuiteSparse_GPURuntime)

if(SuiteSparse_GPU)

find_package(CUDA)
include_directories("${CUDA_INCLUDE_DIRS}")
include_directories("${suitesparseconfig_SOURCE_DIR}")

# Headers:
include_directories("${SuiteSparse_GPURuntime_SOURCE_DIR}/Include")
FILE(GLOB LIBHDRS "${SuiteSparse_GPURuntime_SOURCE_DIR}/Include/*.h")
FILE(GLOB LIBSRCS "${SuiteSparse_GPURuntime_SOURCE_DIR}/Source/*.cpp")

# Create the library
declare_suitesparse_library(suitesparse_gpu_runtime "${LIBSRCS}" "${LIBHDRS}")

# We need to set GPU_BLAS for #ifdef statements
target_compile_definitions(suitesparse_gpu_runtime PUBLIC GPU_BLAS)

endif(SuiteSparse_GPU)

