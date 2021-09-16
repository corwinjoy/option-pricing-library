# Microsoft Developer Studio Generated NMAKE File, Format Version 4.00
# ** DO NOT EDIT **

# TARGTYPE "Win32 (x86) Static Library" 0x0104

!IF "$(CFG)" == ""
CFG=options - Win32 Debug
!MESSAGE No configuration specified.  Defaulting to options - Win32 Debug.
!ENDIF 

!IF "$(CFG)" != "options - Win32 Release" && "$(CFG)" !=\
 "options - Win32 Debug"
!MESSAGE Invalid configuration "$(CFG)" specified.
!MESSAGE You can specify a configuration when running NMAKE on this makefile
!MESSAGE by defining the macro CFG on the command line.  For example:
!MESSAGE 
!MESSAGE NMAKE /f "options.mak" CFG="options - Win32 Debug"
!MESSAGE 
!MESSAGE Possible choices for configuration are:
!MESSAGE 
!MESSAGE "options - Win32 Release" (based on "Win32 (x86) Static Library")
!MESSAGE "options - Win32 Debug" (based on "Win32 (x86) Static Library")
!MESSAGE 
!ERROR An invalid configuration is specified.
!ENDIF 

!IF "$(OS)" == "Windows_NT"
NULL=
!ELSE 
NULL=nul
!ENDIF 
################################################################################
# Begin Project
CPP=cl.exe

!IF  "$(CFG)" == "options - Win32 Release"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 0
# PROP BASE Output_Dir "Release"
# PROP BASE Intermediate_Dir "Release"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 0
# PROP Output_Dir "Release"
# PROP Intermediate_Dir "Release"
# PROP Target_Dir ""
OUTDIR=.\Release
INTDIR=.\Release

ALL : "$(OUTDIR)\options.lib"

CLEAN : 
	-@erase ".\Release\options.lib"
	-@erase ".\Release\spread.obj"
	-@erase ".\Release\asian.obj"
	-@erase ".\Release\basket.obj"
	-@erase ".\Release\binary.obj"
	-@erase ".\Release\boot.obj"
	-@erase ".\Release\cmpd.obj"
	-@erase ".\Release\date.obj"
	-@erase ".\Release\euro.obj"
	-@erase ".\Release\excg.obj"
	-@erase ".\Release\faure.obj"
	-@erase ".\Release\optutil.obj"
	-@erase ".\Release\amer.obj"

"$(OUTDIR)" :
    if not exist "$(OUTDIR)/$(NULL)" mkdir "$(OUTDIR)"

# ADD BASE CPP /nologo /W3 /GX /D "WIN32" /D "NDEBUG" /D "_WINDOWS" /YX /c
# ADD CPP /nologo /W3 /GX /D "WIN32" /D "NDEBUG" /D "_WINDOWS" /YX /c
CPP_PROJ=/nologo /ML /W3 /GX /D "WIN32" /D "NDEBUG" /D "_WINDOWS"\
 /Fp"$(INTDIR)/options.pch" /YX /Fo"$(INTDIR)/" /c 
CPP_OBJS=.\Release/
CPP_SBRS=
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
BSC32_FLAGS=/nologo /o"$(OUTDIR)/options.bsc" 
BSC32_SBRS=
LIB32=link.exe -lib
# ADD BASE LIB32 /nologo
# ADD LIB32 /nologo
LIB32_FLAGS=/nologo /out:"$(OUTDIR)/options.lib" 
LIB32_OBJS= \
	"$(INTDIR)/spread.obj" \
	"$(INTDIR)/asian.obj" \
	"$(INTDIR)/basket.obj" \
	"$(INTDIR)/binary.obj" \
	"$(INTDIR)/boot.obj" \
	"$(INTDIR)/cmpd.obj" \
	"$(INTDIR)/date.obj" \
	"$(INTDIR)/euro.obj" \
	"$(INTDIR)/excg.obj" \
	"$(INTDIR)/faure.obj" \
	"$(INTDIR)/optutil.obj" \
	"$(INTDIR)/amer.obj"

"$(OUTDIR)\options.lib" : "$(OUTDIR)" $(DEF_FILE) $(LIB32_OBJS)
    $(LIB32) @<<
  $(LIB32_FLAGS) $(DEF_FLAGS) $(LIB32_OBJS)
<<

!ELSEIF  "$(CFG)" == "options - Win32 Debug"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 1
# PROP BASE Output_Dir "Debug"
# PROP BASE Intermediate_Dir "Debug"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 1
# PROP Output_Dir "Debug"
# PROP Intermediate_Dir "Debug"
# PROP Target_Dir ""
OUTDIR=.\Debug
INTDIR=.\Debug

ALL : "$(OUTDIR)\options.lib"

CLEAN : 
	-@erase ".\Debug\options.lib"
	-@erase ".\Debug\spread.obj"
	-@erase ".\Debug\asian.obj"
	-@erase ".\Debug\basket.obj"
	-@erase ".\Debug\binary.obj"
	-@erase ".\Debug\boot.obj"
	-@erase ".\Debug\cmpd.obj"
	-@erase ".\Debug\date.obj"
	-@erase ".\Debug\euro.obj"
	-@erase ".\Debug\excg.obj"
	-@erase ".\Debug\faure.obj"
	-@erase ".\Debug\optutil.obj"
	-@erase ".\Debug\amer.obj"

"$(OUTDIR)" :
    if not exist "$(OUTDIR)/$(NULL)" mkdir "$(OUTDIR)"

# ADD BASE CPP /nologo /W3 /GX /Z7 /Od /D "WIN32" /D "_DEBUG" /D "_WINDOWS" /YX /c
# ADD CPP /nologo /W3 /GX /Z7 /Od /D "WIN32" /D "_DEBUG" /D "_WINDOWS" /YX /c
CPP_PROJ=/nologo /MLd /W3 /GX /Z7 /Od /D "WIN32" /D "_DEBUG" /D "_WINDOWS"\
 /Fp"$(INTDIR)/options.pch" /YX /Fo"$(INTDIR)/" /c 
CPP_OBJS=.\Debug/
CPP_SBRS=
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
BSC32_FLAGS=/nologo /o"$(OUTDIR)/options.bsc" 
BSC32_SBRS=
LIB32=link.exe -lib
# ADD BASE LIB32 /nologo
# ADD LIB32 /nologo
LIB32_FLAGS=/nologo /out:"$(OUTDIR)/options.lib" 
LIB32_OBJS= \
	"$(INTDIR)/spread.obj" \
	"$(INTDIR)/asian.obj" \
	"$(INTDIR)/basket.obj" \
	"$(INTDIR)/binary.obj" \
	"$(INTDIR)/boot.obj" \
	"$(INTDIR)/cmpd.obj" \
	"$(INTDIR)/date.obj" \
	"$(INTDIR)/euro.obj" \
	"$(INTDIR)/excg.obj" \
	"$(INTDIR)/faure.obj" \
	"$(INTDIR)/optutil.obj" \
	"$(INTDIR)/amer.obj"

"$(OUTDIR)\options.lib" : "$(OUTDIR)" $(DEF_FILE) $(LIB32_OBJS)
    $(LIB32) @<<
  $(LIB32_FLAGS) $(DEF_FLAGS) $(LIB32_OBJS)
<<

!ENDIF 

.c{$(CPP_OBJS)}.obj:
   $(CPP) $(CPP_PROJ) $<  

.cpp{$(CPP_OBJS)}.obj:
   $(CPP) $(CPP_PROJ) $<  

.cxx{$(CPP_OBJS)}.obj:
   $(CPP) $(CPP_PROJ) $<  

.c{$(CPP_SBRS)}.sbr:
   $(CPP) $(CPP_PROJ) $<  

.cpp{$(CPP_SBRS)}.sbr:
   $(CPP) $(CPP_PROJ) $<  

.cxx{$(CPP_SBRS)}.sbr:
   $(CPP) $(CPP_PROJ) $<  

################################################################################
# Begin Target

# Name "options - Win32 Release"
# Name "options - Win32 Debug"

!IF  "$(CFG)" == "options - Win32 Release"

!ELSEIF  "$(CFG)" == "options - Win32 Debug"

!ENDIF 

################################################################################
# Begin Source File

SOURCE=.\spread.cpp
DEP_CPP_SPREA=\
	".\option.h"\
	".\date.h"\
	

"$(INTDIR)\spread.obj" : $(SOURCE) $(DEP_CPP_SPREA) "$(INTDIR)"


# End Source File
################################################################################
# Begin Source File

SOURCE=.\asian.cpp
DEP_CPP_ASIAN=\
	".\option.h"\
	".\date.h"\
	

"$(INTDIR)\asian.obj" : $(SOURCE) $(DEP_CPP_ASIAN) "$(INTDIR)"


# End Source File
################################################################################
# Begin Source File

SOURCE=.\basket.cpp
DEP_CPP_BASKE=\
	".\option.h"\
	".\date.h"\
	

"$(INTDIR)\basket.obj" : $(SOURCE) $(DEP_CPP_BASKE) "$(INTDIR)"


# End Source File
################################################################################
# Begin Source File

SOURCE=.\binary.cpp
DEP_CPP_BINAR=\
	".\option.h"\
	".\date.h"\
	

"$(INTDIR)\binary.obj" : $(SOURCE) $(DEP_CPP_BINAR) "$(INTDIR)"


# End Source File
################################################################################
# Begin Source File

SOURCE=.\boot.cpp
DEP_CPP_BOOT_=\
	".\option.h"\
	".\date.h"\
	

"$(INTDIR)\boot.obj" : $(SOURCE) $(DEP_CPP_BOOT_) "$(INTDIR)"


# End Source File
################################################################################
# Begin Source File

SOURCE=.\cmpd.cpp
DEP_CPP_CMPD_=\
	".\option.h"\
	".\date.h"\
	

"$(INTDIR)\cmpd.obj" : $(SOURCE) $(DEP_CPP_CMPD_) "$(INTDIR)"


# End Source File
################################################################################
# Begin Source File

SOURCE=.\date.cpp
DEP_CPP_DATE_=\
	".\date.h"\
	

"$(INTDIR)\date.obj" : $(SOURCE) $(DEP_CPP_DATE_) "$(INTDIR)"


# End Source File
################################################################################
# Begin Source File

SOURCE=.\date.h

!IF  "$(CFG)" == "options - Win32 Release"

!ELSEIF  "$(CFG)" == "options - Win32 Debug"

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=.\euro.cpp
DEP_CPP_EURO_=\
	".\option.h"\
	".\date.h"\
	

"$(INTDIR)\euro.obj" : $(SOURCE) $(DEP_CPP_EURO_) "$(INTDIR)"


# End Source File
################################################################################
# Begin Source File

SOURCE=.\excg.cpp
DEP_CPP_EXCG_=\
	".\option.h"\
	".\date.h"\
	

"$(INTDIR)\excg.obj" : $(SOURCE) $(DEP_CPP_EXCG_) "$(INTDIR)"


# End Source File
################################################################################
# Begin Source File

SOURCE=.\faure.cpp
DEP_CPP_FAURE=\
	".\option.h"\
	".\date.h"\
	

"$(INTDIR)\faure.obj" : $(SOURCE) $(DEP_CPP_FAURE) "$(INTDIR)"


# End Source File
################################################################################
# Begin Source File

SOURCE=.\option.h

!IF  "$(CFG)" == "options - Win32 Release"

!ELSEIF  "$(CFG)" == "options - Win32 Debug"

!ENDIF 

# End Source File
################################################################################
# Begin Source File

SOURCE=.\optutil.cpp
DEP_CPP_OPTUT=\
	".\option.h"\
	".\date.h"\
	

"$(INTDIR)\optutil.obj" : $(SOURCE) $(DEP_CPP_OPTUT) "$(INTDIR)"


# End Source File
################################################################################
# Begin Source File

SOURCE=.\amer.cpp
DEP_CPP_AMER_=\
	".\option.h"\
	".\date.h"\
	

"$(INTDIR)\amer.obj" : $(SOURCE) $(DEP_CPP_AMER_) "$(INTDIR)"


# End Source File
# End Target
# End Project
################################################################################
