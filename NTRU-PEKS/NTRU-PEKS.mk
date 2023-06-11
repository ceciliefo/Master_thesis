##
## Auto Generated makefile by CodeLite IDE
## any manual changes will be erased      
##
## Debug
ProjectName            :=NTRU-PEKS
ConfigurationName      :=Debug
WorkspaceConfiguration :=Debug
WorkspacePath          :=/Users/cecilie/Documents/Finished_Code_Masterproject/NTRU-PEKS
ProjectPath            :=/Users/cecilie/Documents/Finished_Code_Masterproject/NTRU-PEKS
IntermediateDirectory  :=build-$(WorkspaceConfiguration)
OutDir                 :=$(IntermediateDirectory)
CurrentFileName        :=
CurrentFilePath        :=
CurrentFileFullPath    :=
User                   :=Cecilie Fougner
Date                   :=11/06/2023
CodeLitePath           :="/Users/cecilie/Library/Application Support/CodeLite"
MakeDirCommand         :=mkdir -p
LinkerName             :=g++
SharedObjectLinkerName :=g++ -dynamiclib -fPIC
ObjectSuffix           :=.o
DependSuffix           :=.o.d
PreprocessSuffix       :=.o.i
DebugSwitch            :=-gstab
IncludeSwitch          :=-I
LibrarySwitch          :=-l
OutputSwitch           :=-o 
LibraryPathSwitch      :=-L
PreprocessorSwitch     :=-D
SourceSwitch           :=-c 
OutputDirectory        :=/Users/cecilie/Documents/Finished_Code_Masterproject/NTRU-PEKS/build-$(WorkspaceConfiguration)/bin
OutputFile             :=build-$(WorkspaceConfiguration)/bin/$(ProjectName)
Preprocessors          :=
ObjectSwitch           :=-o 
ArchiveOutputSwitch    := 
PreprocessOnlySwitch   :=-E 
ObjectsFileList        :=$(IntermediateDirectory)/ObjectsList.txt
PCHCompileFlags        :=
LinkOptions            :=  -lntl -lgmp
IncludePath            :=  $(IncludeSwitch). $(IncludeSwitch). 
IncludePCH             := 
RcIncludePath          := 
Libs                   := 
ArLibs                 :=  
LibPath                := $(LibraryPathSwitch). 

##
## Common variables
## AR, CXX, CC, AS, CXXFLAGS and CFLAGS can be overridden using an environment variable
##
AR       := ar rcus
CXX      := g++
CC       := gcc
CXXFLAGS := -Wall -O2 -std=c++11 -Wall   $(Preprocessors)
CFLAGS   :=  -O2 -Wall $(Preprocessors)
ASFLAGS  := 
AS       := as


##
## User defined environment variables
##
CodeLiteDir:=/Applications/codelite.app/Contents/SharedSupport/
Objects0=$(IntermediateDirectory)/Sampling.cc$(ObjectSuffix) $(IntermediateDirectory)/io.cc$(ObjectSuffix) $(IntermediateDirectory)/FFT.cc$(ObjectSuffix) $(IntermediateDirectory)/Biometrics.cc$(ObjectSuffix) $(IntermediateDirectory)/Random.cc$(ObjectSuffix) $(IntermediateDirectory)/Scheme.cc$(ObjectSuffix) $(IntermediateDirectory)/Algebra.cc$(ObjectSuffix) $(IntermediateDirectory)/PEKS.cc$(ObjectSuffix) 



Objects=$(Objects0) 

##
## Main Build Targets 
##
.PHONY: all clean PreBuild PrePreBuild PostBuild MakeIntermediateDirs
all: MakeIntermediateDirs $(OutputFile)

$(OutputFile): $(IntermediateDirectory)/.d $(Objects) 
	@$(MakeDirCommand) "$(IntermediateDirectory)"
	@echo "" > $(IntermediateDirectory)/.d
	@echo $(Objects0)  > $(ObjectsFileList)
	$(LinkerName) $(OutputSwitch)$(OutputFile) @$(ObjectsFileList) $(LibPath) $(Libs) $(LinkOptions)

MakeIntermediateDirs:
	@$(MakeDirCommand) "$(IntermediateDirectory)"
	@$(MakeDirCommand) "$(OutputDirectory)"

$(IntermediateDirectory)/.d:
	@$(MakeDirCommand) "$(IntermediateDirectory)"

PreBuild:


##
## Objects
##
$(IntermediateDirectory)/Sampling.cc$(ObjectSuffix): Sampling.cc 
	$(CXX) $(IncludePCH) $(SourceSwitch) "/Users/cecilie/Documents/Finished_Code_Masterproject/NTRU-PEKS/Sampling.cc" $(CXXFLAGS) $(ObjectSwitch)$(IntermediateDirectory)/Sampling.cc$(ObjectSuffix) $(IncludePath)
$(IntermediateDirectory)/Sampling.cc$(PreprocessSuffix): Sampling.cc
	$(CXX) $(CXXFLAGS) $(IncludePCH) $(IncludePath) $(PreprocessOnlySwitch) $(OutputSwitch) $(IntermediateDirectory)/Sampling.cc$(PreprocessSuffix) Sampling.cc

$(IntermediateDirectory)/io.cc$(ObjectSuffix): io.cc 
	$(CXX) $(IncludePCH) $(SourceSwitch) "/Users/cecilie/Documents/Finished_Code_Masterproject/NTRU-PEKS/io.cc" $(CXXFLAGS) $(ObjectSwitch)$(IntermediateDirectory)/io.cc$(ObjectSuffix) $(IncludePath)
$(IntermediateDirectory)/io.cc$(PreprocessSuffix): io.cc
	$(CXX) $(CXXFLAGS) $(IncludePCH) $(IncludePath) $(PreprocessOnlySwitch) $(OutputSwitch) $(IntermediateDirectory)/io.cc$(PreprocessSuffix) io.cc

$(IntermediateDirectory)/FFT.cc$(ObjectSuffix): FFT.cc 
	$(CXX) $(IncludePCH) $(SourceSwitch) "/Users/cecilie/Documents/Finished_Code_Masterproject/NTRU-PEKS/FFT.cc" $(CXXFLAGS) $(ObjectSwitch)$(IntermediateDirectory)/FFT.cc$(ObjectSuffix) $(IncludePath)
$(IntermediateDirectory)/FFT.cc$(PreprocessSuffix): FFT.cc
	$(CXX) $(CXXFLAGS) $(IncludePCH) $(IncludePath) $(PreprocessOnlySwitch) $(OutputSwitch) $(IntermediateDirectory)/FFT.cc$(PreprocessSuffix) FFT.cc

$(IntermediateDirectory)/Biometrics.cc$(ObjectSuffix): Biometrics.cc 
	$(CXX) $(IncludePCH) $(SourceSwitch) "/Users/cecilie/Documents/Finished_Code_Masterproject/NTRU-PEKS/Biometrics.cc" $(CXXFLAGS) $(ObjectSwitch)$(IntermediateDirectory)/Biometrics.cc$(ObjectSuffix) $(IncludePath)
$(IntermediateDirectory)/Biometrics.cc$(PreprocessSuffix): Biometrics.cc
	$(CXX) $(CXXFLAGS) $(IncludePCH) $(IncludePath) $(PreprocessOnlySwitch) $(OutputSwitch) $(IntermediateDirectory)/Biometrics.cc$(PreprocessSuffix) Biometrics.cc

$(IntermediateDirectory)/Random.cc$(ObjectSuffix): Random.cc 
	$(CXX) $(IncludePCH) $(SourceSwitch) "/Users/cecilie/Documents/Finished_Code_Masterproject/NTRU-PEKS/Random.cc" $(CXXFLAGS) $(ObjectSwitch)$(IntermediateDirectory)/Random.cc$(ObjectSuffix) $(IncludePath)
$(IntermediateDirectory)/Random.cc$(PreprocessSuffix): Random.cc
	$(CXX) $(CXXFLAGS) $(IncludePCH) $(IncludePath) $(PreprocessOnlySwitch) $(OutputSwitch) $(IntermediateDirectory)/Random.cc$(PreprocessSuffix) Random.cc

$(IntermediateDirectory)/Scheme.cc$(ObjectSuffix): Scheme.cc 
	$(CXX) $(IncludePCH) $(SourceSwitch) "/Users/cecilie/Documents/Finished_Code_Masterproject/NTRU-PEKS/Scheme.cc" $(CXXFLAGS) $(ObjectSwitch)$(IntermediateDirectory)/Scheme.cc$(ObjectSuffix) $(IncludePath)
$(IntermediateDirectory)/Scheme.cc$(PreprocessSuffix): Scheme.cc
	$(CXX) $(CXXFLAGS) $(IncludePCH) $(IncludePath) $(PreprocessOnlySwitch) $(OutputSwitch) $(IntermediateDirectory)/Scheme.cc$(PreprocessSuffix) Scheme.cc

$(IntermediateDirectory)/Algebra.cc$(ObjectSuffix): Algebra.cc 
	$(CXX) $(IncludePCH) $(SourceSwitch) "/Users/cecilie/Documents/Finished_Code_Masterproject/NTRU-PEKS/Algebra.cc" $(CXXFLAGS) $(ObjectSwitch)$(IntermediateDirectory)/Algebra.cc$(ObjectSuffix) $(IncludePath)
$(IntermediateDirectory)/Algebra.cc$(PreprocessSuffix): Algebra.cc
	$(CXX) $(CXXFLAGS) $(IncludePCH) $(IncludePath) $(PreprocessOnlySwitch) $(OutputSwitch) $(IntermediateDirectory)/Algebra.cc$(PreprocessSuffix) Algebra.cc

$(IntermediateDirectory)/PEKS.cc$(ObjectSuffix): PEKS.cc 
	$(CXX) $(IncludePCH) $(SourceSwitch) "/Users/cecilie/Documents/Finished_Code_Masterproject/NTRU-PEKS/PEKS.cc" $(CXXFLAGS) $(ObjectSwitch)$(IntermediateDirectory)/PEKS.cc$(ObjectSuffix) $(IncludePath)
$(IntermediateDirectory)/PEKS.cc$(PreprocessSuffix): PEKS.cc
	$(CXX) $(CXXFLAGS) $(IncludePCH) $(IncludePath) $(PreprocessOnlySwitch) $(OutputSwitch) $(IntermediateDirectory)/PEKS.cc$(PreprocessSuffix) PEKS.cc

##
## Clean
##
clean:
	$(RM) -r $(IntermediateDirectory)


