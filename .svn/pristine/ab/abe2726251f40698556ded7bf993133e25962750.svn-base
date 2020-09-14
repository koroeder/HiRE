PROJECT(GPUQREngine)

if(SuiteSparse_GPU)

include_directories("${suitesparseconfig_SOURCE_DIR}")
include_directories("${AMD_SOURCE_DIR}/Include")
include_directories("${CAMD_SOURCE_DIR}/Include")
include_directories("${COLAMD_SOURCE_DIR}/Include")
include_directories("${CCOLAMD_SOURCE_DIR}/Include")

# Headers:
include_directories("${GPUQREngine_SOURCE_DIR}/Include")
include_directories("${SuiteSparse_GPURuntime_SOURCE_DIR}/Include")
FILE(GLOB LIBHDRS "${GPUQREngine_SOURCE_DIR}/Include/*.h")
FILE(GLOB LIBSRCS "${GPUQREngine_SOURCE_DIR}/Source/*.cpp"
                  "${GPUQREngine_SOURCE_DIR}/Source/BucketList/*.cpp"
                  "${GPUQREngine_SOURCE_DIR}/Source/LLBundle/*.cpp"
                  "${GPUQREngine_SOURCE_DIR}/Source/Scheduler/*.cpp"
                  "${GPUQREngine_SOURCE_DIR}/Source/TaskDescriptor/*.cpp")

# Build the kernels 
find_package(CUDA)
set(CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS};-O3 -gencode=arch=compute_20,code=sm_20 -gencode=arch=compute_30,code=sm_30 -gencode=arch=compute_35,code=sm_35)
cuda_include_directories("${suitesparseconfig_SOURCE_DIR}")
cuda_include_directories("${GPUQREngine_SOURCE_DIR}/Include")
cuda_include_directories("${SuiteSparse_GPURuntime_SOURCE_DIR}/Include")
# Compile kernels
cuda_compile(GPUQREngine_CU_FILES "${GPUQREngine_SOURCE_DIR}/Source/GPUQREngine_UberKernel.cu")
cuda_add_library(gpuqr_cuda ${GPUQREngine_CU_FILES})
cuda_add_cublas_to_target(gpuqr_cuda)
install_suitesparse_project(gpuqr_cuda "")
# Create the library
include_directories("${CUDA_INCLUDE_DIRS}")
declare_suitesparse_library(gpu_qr_engine "${LIBSRCS}" "${LIBHDRS}"
   TARGET_PRIVATE_LINK	${SuiteSparse_LINKER_LAPACK_BLAS_LIBS}
   TARGET_PUBLIC_LINK 	${SuiteSparse_LINKER_METIS_LIBS} amd camd colamd ccolamd gpuqr_cuda ${CUDA_CUBLAS_LIBRARIES}
)
# We need to set GPU_BLAS for #ifdef statements in the cholmod code
target_compile_definitions(gpu_qr_engine PUBLIC GPU_BLAS)

endif(SuiteSparse_GPU)

