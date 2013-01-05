#!/bin/zsh


mv src/lib/ofdm/Makefile.am src/lib/ofdm/Makefile.am.bk
sed 's/\$(GNURADIO_CORE_LA)/\$(GNURADIO_CORE_LA) \\/g; /\$(GNURADIO_CORE_LA)/ a \ 
\ \ \ \ \ \ \ \ \$(GR_DIGITAL_LA)' src/lib/ofdm/Makefile.am.bk > src/lib/ofdm/Makefile.am


mv config/Makefile config/Makefile.bk
sed 's/LIBS =  -L\/usr\/local\/lib/LIBS =  -L\/usr\/local\/lib -lgnuradio-digital/g' config/Makefile.bk > config/Makefile
mv src/python/Makefile src/python/Makefile.bk
sed 's/LIBS =  -L\/usr\/local\/lib/LIBS =  -L\/usr\/local\/lib -lgnuradio-digital/g' src/python/Makefile.bk > src/python/Makefile
mv src/lib/ofdm/Makefile src/lib/ofdm/Makefile.bk
sed 's/LIBS =  -L\/usr\/local\/lib/LIBS =  -L\/usr\/local\/lib -lgnuradio-digital/g' src/lib/ofdm/Makefile.bk > src/lib/ofdm/Makefile
mv src/lib/qam/Makefile src/lib/qam/Makefile.bk
sed 's/LIBS =  -L\/usr\/local\/lib/LIBS =  -L\/usr\/local\/lib -lgnuradio-digital/g' src/lib/qam/Makefile.bk > src/lib/qam/Makefile
mv src/lib/rscode/Makefile src/lib/rscode/Makefile.bk
sed 's/LIBS =  -L\/usr\/local\/lib/LIBS =  -L\/usr\/local\/lib -lgnuradio-digital/g' src/lib/rscode/Makefile.bk > src/lib/rscode/Makefile
mv src/lib/spiral/Makefile src/lib/spiral/Makefile.bk
sed 's/LIBS =  -L\/usr\/local\/lib/LIBS =  -L\/usr\/local\/lib -lgnuradio-digital/g' src/lib/spiral/Makefile.bk > src/lib/spiral/Makefile
mv src/lib/util/Makefile src/lib/util/Makefile.bk
sed 's/LIBS =  -L\/usr\/local\/lib/LIBS =  -L\/usr\/local\/lib -lgnuradio-digital/g' src/lib/util/Makefile.bk > src/lib/util/Makefile
mv src/lib/Makefile src/lib/Makefile.bk
sed 's/LIBS =  -L\/usr\/local\/lib/LIBS =  -L\/usr\/local\/lib -lgnuradio-digital/g' src/lib/Makefile.bk > src/lib/Makefile
mv src/Makefile src/Makefile.bk
sed 's/LIBS =  -L\/usr\/local\/lib/LIBS =  -L\/usr\/local\/lib -lgnuradio-digital/g' src/Makefile.bk > src/Makefile
mv Makefile Makefile.bk
sed 's/LIBS =  -L\/usr\/local\/lib/LIBS =  -L\/usr\/local\/lib -lgnuradio-digital/g' Makefile.bk > Makefile


mv Makefile.common Makefile.common.bk
sed 's/-I\$(GNURADIO_CORE_INCLUDEDIR)\/swig/-I\$(GNURADIO_CORE_INCLUDEDIR)\/swig \\/g; -I\$(GNURADIO_CORE_INCLUDEDIR)\/swig/ a \
\ \ \ \ \ \ -I\/usr\/local\/include\/gruel\/swig' Makefile.common.bk > Makefile.common






