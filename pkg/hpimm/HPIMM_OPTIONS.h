#ifndef HPIMM_OPTIONS_H
#define HPIMM_OPTIONS_H
#include "PACKAGES_CONFIG.h"
#include "CPP_OPTIONS.h"

CBOP
C !ROUTINE: HPIMM_OPTIONS.h
C !INTERFACE:
C #include "HPIMM_OPTIONS.h"

C !DESCRIPTION:
C *==================================================================*
C | CPP options file for pkg "HPIMM":
C | Control which optional features to compile in this package code.
C *==================================================================*
CEOP

#ifdef ALLOW_HPIMM
C Place CPP define/undef flag here

C to reduce memory storage, disable unused array with those CPP flags :
#undef HPIMM_3D_STATE
#define HPIMM_2D_STATE
#undef HPIMM_TENDENCY

#undef HPIMM_SPECIAL_COMPILE_OPTION1

#define HPIMM_SPECIAL_COMPILE_OPTION2

#endif /* ALLOW_HPIMM */
#endif /* HPIMM_OPTIONS_H */
