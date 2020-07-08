/******************************************************************************
 * A program for geological map topology analysis.
 *
 * Author: Vitaliy Ogarko, vogarko@gmail.com
 *******************************************************************************/

#include <iostream>
#include <stdexcept>

#include "file_utils.h"
#include "parameters_reader.h"
#include "tester.h"
#include "topology_analyzer.h"

const bool RUN_TESTS = false;

// Main entry point.
int main(int argc, char *argv[]) {
  if (RUN_TESTS) {
    ConverterLib::Test_FaultAndPolygonIntersecting();
    ConverterLib::Test_FaultsAreIntersecting();
    return 0;
  }

  try {
    // Retrieve parameters filename from a command line.
    // std::string parfile_path;
    // if (argc < 2) {
    //   std::cout << "Error: parameter file is not specified!" << std::endl;
    //   std::cout << "Usage: map2model <Parfile>" << std::endl;
    //   exit(0);
    // } else {
    //   parfile_path = argv[1];
    // }

    // Reading the parameters file.
    ConverterLib::Parameters par;
    // const std::string parameter_lines = par.Read(parfile_path);
    const std::string parameter_lines = par.directRead(
        "../model-test/output", "../model-test/tmp/geology.csv",
        "../model-test/tmp/faults.csv", "../model-test/tmp/mindep.csv");

    // Create the output data folder.
    FileUtils::CreateDirectory(par.path_output.c_str());

    // Copy parameters file to the output folder (for reference).
    // FileUtils::CopyFile(parfile_path.c_str(),
    //                     (par.path_output + "/Parfile").c_str());

    //----------------------------------------------------------------------------------
    // Main calculations start here.
    ConverterLib::TopologyAnalyzer topo_analyzer;
    topo_analyzer.Initialize(par);

    if (par.subregion_size_x > 0 && par.subregion_size_y > 0) {
      topo_analyzer.AnalyzeLocalTopology(par);
    } else {
      topo_analyzer.AnalyzeGlobalTopology(par, parameter_lines);
    }
  } catch (const std::exception &e) {
    std::cerr << "Unexpected exception in main caught: " << e.what()
              << std::endl;
  }

  std::cout << "End." << std::endl;
  return 0;
}
