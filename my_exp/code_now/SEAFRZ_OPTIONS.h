#ifndef SEAFRZ_OPTIONS_H
#define SEAFRZ_OPTIONS_H
#include "PACKAGES_CONFIG.h"
#include "CPP_OPTIONS.h"

CBOP
C !ROUTINE: SEAFRZ_OPTIONS.h
C !INTERFACE:
C #include "SEAFRZ_OPTIONS.h"

C !DESCRIPTION:
C *==================================================================*
C | CPP options file for pkg "seafrz":
C | Control which optional features to compile in this package code.
C *==================================================================*
CEOP

#ifdef ALLOW_SEAFRZ
C Place CPP define/undef flag here

C to reduce memory storage, disable unused array with those CPP flags :
#undef SEAFRZ_3D_STATE
#define SEAFRZ_2D_STATE
#undef SEAFRZ_TENDENCY

#undef SFZ_SPECIAL_COMPILE_OPTION1

#define SFZ_SPECIAL_COMPILE_OPTION2

#endif /* ALLOW_SEAFRZ */
#endif /* SEAFRZ_OPTIONS_H */
