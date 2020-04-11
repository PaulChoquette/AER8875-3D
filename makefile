#----------------------------------------------------+
# C++ makefile with recursive directory file search
# by Matthieu Parenteau
# PI4 2018
#----------------------------------------------------+

# Operating System
SHELL = /bin/sh
OS := $(shell uname)

# Text Variables
comma := ,
empty :=
space := $(empty) $(empty)
sharp := "\#"
quote := \"
newline := ""

#-------------------------------------------------------------------------------------------------------------------+
# Directory Variables

DeactivatedFiles := $(shell grep -v '^//' .makeignore | grep -v '^$$')

# Names of all root folders (these contains source code in any tree structure of visible folders)
SourceDirs := $(shell tree -dfi --noreport)

# Determine the path for the source code folders in the .debug folder (will mimic the visible root tree structure)
DebugDirs := $(subst . ,,$(subst ./,.debug/,$(SourceDirs)))

# Determine the path for the source code folders in the .release folder (will mimic the visible root tree structure)
ReleaseDirs := $(subst . ,,$(subst ./,.release/,$(SourceDirs)))

# Concatenate all of the previous 3 folder name variables
AllDirs := $(SourceDirs) $(ReleaseDirs) $(DebugDirs)

#-------------------------------------------------------------------------------------------------------------------+
#-----------------------------------------------------------------------------------------------------+
# File Variables

# Shared Library
current_dir := $(shell basename $(CURDIR))

Executable = $(current_dir)

ExecutableSourceFile = main.cpp

# Executable debug object file
ExecutableDebugObjectFile = $(subst .cpp,.o,$(addprefix .debug/,$(ExecutableSourceFile)))

# Executable release object file
ExecutableReleaseObjectFile = $(subst .cpp,.o,$(addprefix .release/,$(ExecutableSourceFile)))

# Executable shared object file
ExecutableSharedObjectFile = $(subst .cpp,.o,$(addprefix .shared/,$(ExecutableSourceFile)))

# Source files (.cpp)
SourceFiles := $(filter-out $(ExecutableSourceFile),$(subst ./,,$(subst .//,,$(shell find ./ -regex .*.cpp))))

# Filtered Source files (.cpp with deactivated files)
ActiveSourceFiles = $(filter-out $(DeactivatedFiles),$(SourceFiles))

# Debug object files (names of the object files that are produced by make debug)
DebugObjectFiles := $(addprefix .debug/,$(subst .cpp,.o,$(ActiveSourceFiles)))

# Release object files (names of the object files that are produced by make release)
ReleaseObjectFiles := $(addprefix .release/,$(subst .cpp,.o,$(ActiveSourceFiles)))

# Shared object files (names of the object files that are produced by make shared)
SharedObjectFiles := $(addprefix .shared/,$(subst .cpp,.o,$(ActiveSourceFiles)))

# Concatenate all object files
AnyObjectFiles := $(notdir $(shell find ./ -regex .*.o))

#-----------------------------------------------------------------------------------------------------+
#--------------------------------------------------------------------------------------------------------------------------------------+
# Compilation Variables

# Default file search locations (All directories)
VPATH := $(AllDirs)

# Default compilation configuration
CXX = mpic++
CXXFLAGS += -std=c++11 -Wall -Wno-unused-function -Wno-strict-overflow -Wno-unused-variable

DEBUGFLAGS += -Og -g -pg
#RELEASEFLAGS += -O3 -fopenmp
SHAREDFLAGS += -shared -Wl,--export-dynamic,--no-as-needed -fpic

# Include search paths
METIS_INCLUDE_PATH = $(METISROOT)/include
TECIO_INCLUDE_PATH = $(TECIOROOT)/include

# Included Libraries
METISLIBS += -L$(METISROOT)/lib -lmetis
TECIOLIBS += -L$(TECIOROOT)/lib -ltecio -Lm -pthread
#--------------------------------------------------------------------------------------------------------------------------------------+
#---------------------------------------------------------------------------------------------------+
# Targets

all : release $(ReleaseObjectFiles)

debug : .debug  begun $(DebugObjectFiles) $(ExecutableDebugObjectFile)
	@printf '   Linking Debug...'
	@$(CXX) $(CXXFLAGS) $(DEBUGFLAGS) $(DebugObjectFiles) $(ExecutableDebugObjectFile) -o $(addprefix bin/,$(Executable)) $(METISLIBS) $(TECIOLIBS)
	@printf 'Done'
	@printf '\n'

release : .release begun $(ReleaseObjectFiles) $(ExecutableReleaseObjectFile)
	@printf '   Linking Release...'
	@$(CXX) $(CXXFLAGS) $(RELEASEFLAGS) $(ReleaseObjectFiles) $(ExecutableReleaseObjectFile) -o $(addprefix bin/,$(Executable)) $(METISLIBS) $(TECIOLIBS)
	@printf 'Done'
	@printf '\n'

shared : .shared begun $(SharedObjectFiles) $(ExecutableSharedObjectFile)
	@printf '   Linking Shared...'
	@$(CXX) $(CXXFLAGS) $(RELEASEFLAGS) $(SHAREDFLAGS) $(SharedObjectFiles) $(ExecutableSharedObjectFile) -o $(addprefix lib/,$(Executable).so) $(METISLIBS) $(TECIOLIBS)
	@printf 'Done'
	@printf '\n'

reset : clean
	@$(shell reset)

verify : release $(ReleaseObjectFiles)

#---------------------------------------------------------------------------------------------------+
#---------------------------------------------------------------------------------------+
# Pattern Rules

.debug/%.o : %.cpp
	@$(CXX) -c $(CXXFLAGS) $(DEBUGFLAGS)  -I$(subst $(space), -I,$(AllDirs)) -I$(METIS_INCLUDE_PATH) -I$(TECIO_INCLUDE_PATH) $< -o $@
	@echo '   Pattern Rule | Compiling | '$(CXXFLAGS) $(DEBUGFLAGS) ' | ' $<' ... Done'

.release/%.o : %.cpp
	@$(CXX) -c $(CXXFLAGS) $(RELEASEFLAGS)  -I$(subst $(space), -I,$(AllDirs)) -I$(METIS_INCLUDE_PATH) -I$(TECIO_INCLUDE_PATH) $< -o $@
	@echo '   Pattern Rule | Compiling | '$(CXXFLAGS) $(RELEASEFLAGS) ' | ' $<' ... Done '

.shared/%.so : %.o
	@$(CXX) $(SHAREDFLAGS)  $< -o $@ $(METISLIBS) -I$(METIS_INCLUDE_PATH) $< -o $@

.shared/%.o : %.cpp
	@$(CXX) -c $(CXXFLAGS) $(RELEASEFLAGS) $(SHAREDFLAGS)  -I$(subst $(space), -I,$(AllDirs)) -I$(METIS_INCLUDE_PATH) -I$(TECIO_INCLUDE_PATH) $< -o $@
	@echo '   Pattern Rule | Compiling | '$(CXXFLAGS) $(RELEASEFLAGS) $(SHAREDFLAGS)' | ' $<' ... Done '


#---------------------------------------------------------------------------------------+
#--------------------------------------------------------------------------------+
# Phony Targets

.PHONY : clean cleandebug cleanrelease cleanshared begin end

begun :

clean : cleandebug cleanrelease cleanshared
	@-rm -rf $(AnyObjectFiles)
	@-rm -f bin/$(current_dir)
	@-rm -f lib/$(current_dir).so

cleandebug :
	@-rm -rf .debug .debugtimestamp
	@-rm -f $(addprefix .debug/,$(Executable))

cleanrelease :
	@-rm -rf .release .releasetimestamp
	@-rm -f $(addprefix .release/,$(Executable))

cleanshared :
	@-rm -rf .shared .sharedtimestamp
	@-rm -f $(addprefix .shared/,$(Executable))

.debug : .debugtimestamp

.debugtimestamp :
	@mkdir -p .debug $(DebugDirs)
	@mkdir -p bin
#	@touch .debugtimestamp

.release : .releasetimestamp

.releasetimestamp :
	@mkdir -p .release $(ReleaseDirs)
	@mkdir -p bin

.shared : .sharedtimestamp

.sharedtimestamp :
	@mkdir -p .shared $(SharedDirs)
	@mkdir -p lib

#--------------------------------------------------------------------------------+
#--------------------------------------+
# Information

sourcedirs :
	@printf '%s\n' $(SourceDirs)

sourcefiles :
	@printf '%s\n' $(SourceFiles)

activesourcefiles :
	@printf '%s\n' $(ActiveSourceFiles)

deactivatedfiles :
	@printf '%s\n' $(DeactivatedFiles)

debugobjectfiles :
	@printf '%s\n' $(DebugObjectFiles)

releaseobjectfiles :
	@printf '%s\n' $(ReleaseObjectFiles)

os :
	@echo $(OS)

executable:
	@echo $(current_dir)
#--------------------------------------+