# routines for myblas
# for now we add an extra library which is compiled within the project.
# The only reason for doing that is to ensure that compile flags are identical as
# in the rest of the project. It should in principle work without this trick
# for BLAS/MYBLAS (see BLAS/LAPACK  system libraries). My feeling is that 
# the compiler flag issue is a CHARMM artifact (vr274)

option(WITH_MYBLAS "Compile own blas  (needed for charmm, can cause problems with gfortran 4.7)" ON)
if(WITH_MYBLAS)
    file(GLOB MYBLAS_SOURCES ${HiRE_ROOT}/BLAS/*.f)
    file(GLOB NOT_MYBLAS_SOURCES 
      ${HiRE_ROOT}/BLAS/zhpmv.f
      ${HiRE_ROOT}/BLAS/ztpsv.f
      ${HiRE_ROOT}/BLAS/ztrsv.f
      ${HiRE_ROOT}/BLAS/zgbmv.f
      ${HiRE_ROOT}/BLAS/zgerc.f
      ${HiRE_ROOT}/BLAS/zgemv.f
      ${HiRE_ROOT}/BLAS/ztbsv.f
      ${HiRE_ROOT}/BLAS/ztrsm.f
      ${HiRE_ROOT}/BLAS/zher2k.f
      ${HiRE_ROOT}/BLAS/zher2.f
      ${HiRE_ROOT}/BLAS/zgemm.f
      ${HiRE_ROOT}/BLAS/zhbmv.f
      ${HiRE_ROOT}/BLAS/zhpr2.f
      ${HiRE_ROOT}/BLAS/zhpr.f
      ${HiRE_ROOT}/BLAS/ztbmv.f
      ${HiRE_ROOT}/BLAS/zherk.f
      ${HiRE_ROOT}/BLAS/ztrmm.f
      ${HiRE_ROOT}/BLAS/ztpmv.f
      ${HiRE_ROOT}/BLAS/zher.f
      ${HiRE_ROOT}/BLAS/ztrmv.f
      ${HiRE_ROOT}/BLAS/zhemm.f
      ${HiRE_ROOT}/BLAS/zhemv.f
    )
    
    list(REMOVE_ITEM MYBLAS_SOURCES ${NOT_MYBLAS_SOURCES})
    
    add_library(myblas ${MYBLAS_SOURCES})
    
    message("${HiRE_ROOT}/CMakeModules/FindMYBLAS.cmake: creating BLAS library.")         
    
    SET(MYBLAS_LIBS myblas)
    MARK_AS_ADVANCED(MYBLAS_LIBS)
else(WITH_MYBLAS)
    find_package(BLAS REQUIRED)
    message("using system blas: ${BLAS_LIBRARIES}")
    SET(MYBLAS_LIBS ${BLAS_LIBRARIES})
    if(${MYBLAS_LIBS})
    else(${MYBLAS_LIBS})
       set(MYBLAS_LIBS /state/partition1/software/chemistry/blas/3.7.1/lib/librefblas.a)
    endif(${MYBLAS_LIBS})
endif(WITH_MYBLAS)
