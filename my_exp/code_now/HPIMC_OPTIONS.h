#ifndef HPIMC_OPTIONS_H
#define HPIMC_OPTIONS_H
#include "PACKAGES_CONFIG.h"
#include "CPP_OPTIONS.h"

CBOP
C !ROUTINE: HPIMC_OPTIONS.h
C !INTERFACE:
C #include "HPIMC_OPTIONS.h"

C !DESCRIPTION:
C *==================================================================*
C | CPP options file for pkg "HPIMC":
C | Control which optional features to compile in this package code.
C *==================================================================*
CEOP

#ifdef ALLOW_HPIMC
C Place CPP define/undef flag here

C to reduce memory storage, disable unused array with those CPP flags :
#undef HPIMC_3D_STATE
#define HPIMC_2D_STATE
#undef HPIMC_TENDENCY

#undef HPIMC_SPECIAL_COMPILE_OPTION1

#define HPIMC_SPECIAL_COMPILE_OPTION2

#endif /* ALLOW_HPIMC */
#endif /* HPIMC_OPTIONS_H */
