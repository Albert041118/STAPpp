#pragma once

#include <string>

// Convert Bridge-1 Abaqus .inp into an internal STAP++ .dat file.  The
// converted file is then solved by the ordinary CDomain/ElementGroup/Solver
// pipeline.
bool ConvertBridgeInpToStapDat(const std::string& inpFile, const std::string& datFile);
