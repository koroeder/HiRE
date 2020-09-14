PROJECT(UMFPACK)

include_directories("${suitesparseconfig_SOURCE_DIR}")
include_directories("${AMD_SOURCE_DIR}/Include")
include_directories("${CHOLMOD_SOURCE_DIR}/Include")

# Headers:
include_directories("${UMFPACK_SOURCE_DIR}/Include")
FILE(GLOB LIBHDRS "${UMFPACK_SOURCE_DIR}/Include/*.h")
FILE(GLOB LIBSRCS "${UMFPACK_SOURCE_DIR}/SourceWrappers/*.c")

if(SuiteSparse_GPU)
# GPU build with CUDA
   find_package(CUDA)
# Make sure to include CUDA headers
   include_directories("${CUDA_INCLUDE_DIRS}")
   declare_suitesparse_library(umfpack "${LIBSRCS}" "${LIBHDRS}"
      TARGET_PRIVATE_LINK	${SuiteSparse_LINKER_LAPACK_BLAS_LIBS}
      TARGET_PUBLIC_LINK 	amd camd colamd ccolamd cholmod_cuda ${CUDA_CUBLAS_LIBRARIES} ${SuiteSparse_LINKER_PUBLIC_LIBS}
   )
# We need to set GPU_BLAS for #ifdef statements in the cholmod code
   target_compile_definitions(umfpack PUBLIC GPU_BLAS)
else()
# CPU-only build
   declare_suitesparse_library(umfpack "${LIBSRCS}" "${LIBHDRS}"
   	  TARGET_PRIVATE_LINK	${SuiteSparse_LINKER_LAPACK_BLAS_LIBS}
	  TARGET_PUBLIC_LINK 	amd camd colamd ccolamd cholmod ${SuiteSparse_LINKER_PUBLIC_LIBS}
   )
endif(SuiteSparse_GPU)

