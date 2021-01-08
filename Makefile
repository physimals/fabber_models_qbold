include ${FSLCONFDIR}/default.mk

PROJNAME = fabber_qbold
LIBS = -lfsl-fabberexec -lfsl-fabbercore -lfsl-newimage \
       -lfsl-miscmaths -lfsl-utils -lfsl-NewNifti \
       -lfsl-cprob -lfsl-znz -ldl
XFILES = fabber_qbold
SOFILES = libfsl-fabber_models_qbold.so

# Forward models
OBJS =  fwdmodel_qbold.o

# For debugging:
#OPTFLAGS = -ggdb

# Pass Git revision details
GIT_SHA1:=$(shell git describe --match=NeVeRmAtCh --always --abbrev=40 --dirty)
GIT_DATE:=$(shell git log -1 --format=%ad --date=local)
CXXFLAGS += -DGIT_SHA1=\"${GIT_SHA1}\" -DGIT_DATE="\"${GIT_DATE}\""

#
# Build
#
all: ${XFILES} ${SOFILES}

# models in a library
libfsl-fabber_models_qbold.so : ${OBJS}
	${CXX} ${CXXFLAGS} -shared -o $@ $^ ${LDFLAGS}

# fabber built from the FSL fabbercore library including the models specifieid in this project
fabber_qbold : fabber_main.o libfsl-fabber_models_qbold.so
	${CXX} ${CXXFLAGS} -o $@ $< -lfsl-fabber_models_qbold ${LDFLAGS}
