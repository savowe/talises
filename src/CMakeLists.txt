ADD_EXECUTABLE( raman raman.cpp  )
TARGET_LINK_LIBRARIES( raman myutils ${MUPARSER_LIBRARY} ${GSL_LIBRARY_1} ${GSL_LIBRARY_2})

ADD_LIBRARY( myutils cft_1d.cpp cft_2d.cpp cft_3d.cpp misc.cpp ParameterHandler.cpp pugixml.cpp )
TARGET_LINK_LIBRARIES( myutils m gomp ${FFTW_LIBRARY_1} ${FFTW_LIBRARY_2} )

ADD_EXECUTABLE( gen_psi_0 gen_psi_0.cpp )
TARGET_LINK_LIBRARIES( gen_psi_0 myutils ${MUPARSER_LIBRARY} )