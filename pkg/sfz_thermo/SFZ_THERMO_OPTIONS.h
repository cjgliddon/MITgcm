#ifndef SFZ_THERMO_OPTIONS_H
#define SFZ_THERMO_OPTIONS_H
#include "PACKAGES_CONFIG.h"
#include "CPP_OPTIONS.h"

CBOP
C !ROUTINE: SFZ_THERMO_OPTIONS_H.h
C !INTERFACE:
C #include "SFZ_THERMO_OPTIONS_H.h"

C !DESCRIPTION:
C *==================================================================*
C | CPP options file for pkg "SFz_Thermo":
C | Control which optional features to compile in this package code.
C *==================================================================*
CEOP

#ifdef ALLOW_SFZ_THERMO
C Place CPP define/undef flag here

C to reduce memory storage, disable unused array with those CPP flags :
#define SFZ_THERMO_3D_STATE
#define SFZ_THERMO_2D_STATE
#define SFZ_THERMO_TENDENCY

#undef MYPA_SPECIAL_COMPILE_OPTION1

#define MYPA_SPECIAL_COMPILE_OPTION2

#endif /* ALLOW_SFZ_THERMO */
#endif /* SFZ_THERMO_OPTIONS_H */
