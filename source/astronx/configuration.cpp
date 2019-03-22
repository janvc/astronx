/*
 * Copyright 2015 Jan von Cosel
 *
 * This file is part of astronx.
 *
 * astronx is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * astronx is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have recieved a copy of the GNU General Public License
 * along with astronx. If not, see <http://www.gnu.org/licenses/>.
 *
 */


#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <algorithm>
#include <boost/program_options.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/tokenizer.hpp>
#include "configuration.h"

namespace Astronx
{

int Configuration::init(int argnum, char *arguments[])
{
    namespace po = boost::program_options;

    /*
     * parse the command line options
     */
    try
    {
        po::options_description desc("Allowed options");
        desc.add_options()
                ("help,h", "produce this help message")
                ("verbose,v", "be verbose during calculation")
                ("overwrite,w", "overwrite existing files")
                ("input,i", po::value<std::string>(&m_inputFileName), "the input file")
                ;

        po::positional_options_description pos_desc;
        pos_desc.add("input", -1);

        po::command_line_parser parser(argnum, arguments);
        parser.options(desc).positional(pos_desc).allow_unregistered();
        po::parsed_options options = parser.run();


        po::variables_map vm;
        po::store(options, vm);
        po::notify(vm);

        if (vm.count("help"))
        {
            std::cout << desc << std::endl;
            return -1;
        }
        if (vm.count("input"))
        {
            m_inputFileName = vm["input"].as<std::string>();
        }

        if (vm.count("verbose"))
        {
            m_verbose = true;
        }
        else
        {
            m_verbose = false;
        }

        if (vm.count("overwrite"))
        {
            m_overwrite = true;
        }
        else
        {
            m_overwrite = false;
        }
    }
    catch (const po::error &ex)
    {
        std::cerr << ex.what() << std::endl;
        return -2;
    }

    /*
     * read and parse the input file
     */
    std::ifstream inputFile(m_inputFileName);
    if (! inputFile.good())
    {
        std::cerr << "Input file '" << m_inputFileName << "' does not exist or is not accessible.\n";
        return -3;
    }

    int lastDot = m_inputFileName.find_last_of(".");
    if (lastDot == std::string::npos)
    {
        m_baseName = m_inputFileName;
    }
    else
    {
        m_baseName = m_inputFileName.substr(0, lastDot);
    }

    setDefaults();

    std::string currentLine;
    int lineNumber = 0;
    int coordStartLine = 0;
    while (std::getline(inputFile, currentLine))
    {
        lineNumber++;
        currentLine = currentLine.substr(0, currentLine.find_first_of("#"));

        if (currentLine.size() > 0)
        {
            try
            {
                int status = parseInputLine(currentLine);

                if (status == 1)
                {
                    coordStartLine = lineNumber;
                    break;
                }
            }
            catch (const std::exception &e)
            {
                std::cerr << "Error reading input line number " << lineNumber << "\n" << e.what() << std::endl;
                return -3;
            }
        }
    }

    if (m_tFinal == 0.0)
    {
        std::cerr << "Error: tfinal not specified\n";
        return -4;
    }
    if (m_tOut == 0.0)
    {
        std::cerr << "Error: tout not specified\n";
        return -4;
    }
    if ((m_IntType == Rk4FixM || m_IntType == Rk4FixT) && m_steps)
    {
        std::cerr << "Note: steps file makes no sense with a fixed-step integrator.\n";
        m_steps = false;
    }
    if (m_InitStep == 0.0)
        m_InitStep = m_tOut;

    /*
     * Read the initial coordinates and velocities of the objects
     */
    inputFile.seekg(0);
    for (int i = 0; i < coordStartLine; i++)
        inputFile.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

    while (std::getline(inputFile, currentLine))
    {
        lineNumber++;
        currentLine = currentLine.substr(0, currentLine.find_first_of("#"));

        if (currentLine.find("end_coords") == std::string::npos)    // proceed if this is NOT the end of the coordinates
        {
            if (currentLine.size() > 0)
            {
                try
                {
                    // tokenize string
                    boost::char_separator<char> sep(" ");
                    boost::tokenizer<boost::char_separator<char>> tokens(currentLine, sep);
                    boost::tokenizer<boost::char_separator<char>>::iterator iter = tokens.begin();

                    // read the object data
                    m_names.push_back(*iter);
                    iter++;

                    m_masses.push_back(std::stod(*iter));
                    iter++;

                    m_X0.push_back(std::stod(*iter));
                    iter++;
                    m_Y0.push_back(std::stod(*iter));
                    iter++;
                    m_Z0.push_back(std::stod(*iter));
                    iter++;

                    m_VX0.push_back(std::stod(*iter));
                    iter++;
                    m_VY0.push_back(std::stod(*iter));
                    iter++;
                    m_VZ0.push_back(std::stod(*iter));
                }
                catch (...)
                {
                    std::cerr << "Error reading initial configuration in line number " << lineNumber << "\n" << std::endl;
                    return -5;
                }
            }
        }
        else
            break;
    }

    /*
     * check vector lengths for consistency
     */
    int nnm = m_names.size();
    int nma = m_masses.size();
    int nx0 = m_X0.size();
    int ny0 = m_Y0.size();
    int nz0 = m_Z0.size();
    int nvx = m_VX0.size();
    int nvy = m_VY0.size();
    int nvz = m_VZ0.size();

    if (! (nnm == nma && nnm == nx0 && nnm == ny0 && nnm == nz0 && nnm == nvx && nnm == nvy && nnm == nvz))
    {
        std::cerr << "error in the coordinate section\n";
        return -4;
    }

    m_Nobj = m_masses.size();

    // calculate the total mass
//    m_Mtot = 0.0;
//    for (int i = 0; i < m_Nobj; i++)
//        m_Mtot += m_masses[i];

    // open the output files
    m_outputFileName = m_baseName + ".out";
    m_outputFile = std::ofstream(m_outputFileName);
    m_binTrjFile = std::ofstream(m_baseName + ".bin.trj", std::ios::binary);

    if (m_textTrj)
    {
        m_txtTrjFile = std::ofstream(m_baseName + ".txt.trj");
    }

    if (m_steps)
    {
        m_stepsFile = std::ofstream(m_baseName + ".txt.trj");
    }

    if (m_restart)
    {
        m_restartFile = std::ofstream(m_baseName + ".txt.trj");
    }

    m_N_BS_LargeStep = 0;

    return 0;
}

void Configuration::setDefaults()
{
    m_MaxSubStep = 12;
    m_IncThres   =  8;
    m_Ndigit     = 10;
    m_Nstep      =  5;
//    m_tFinal     =  0.0;
//    m_tOut       =  0.0;
//    m_InitStep   =  0.0;
    m_eps        =  1.0e-6;
    m_epsThres   =  0.9;
    m_MinStep    =  1.0e2;
    m_MaxInc     = 10.0;
    m_RedMin     =  0.9;
    m_RedMax     =  0.01;
    m_ShiftCOM   = false;
    m_ShiftMom   = false;
    m_restart    = false;
    m_steps      = false;
    m_textTrj    = false;
    m_UnResProp  = false;
    m_IntType    = BS;
}

int Configuration::parseInputLine(std::string &inputLine)
{
    namespace ba = boost::algorithm;

    std::string::size_type separator = inputLine.find("=");
    std::string keystring = ba::to_lower_copy(ba::trim_copy(inputLine.substr(0, separator)));
    std::string valuestring = ba::to_lower_copy(inputLine.substr(separator + 1));

    // remove all whitespace before and after the keyword and the value
    keystring.erase(std::remove_if(keystring.begin(), keystring.end(), isspace), keystring.end());
    valuestring.erase(std::remove_if(valuestring.begin(), valuestring.end(), isspace), valuestring.end());

    if (keystring == "tfinal")
    {
        m_tFinal = std::stod(valuestring);
    }
    else if (keystring == "tout")
    {
        m_tOut = std::stod(valuestring);
    }
    else if (keystring == "eps")
    {
        m_eps = std::stod(valuestring);
    }
    else if (keystring == "eps_thres")
    {
        m_epsThres = std::stod(valuestring);
    }
    else if (keystring == "initstep")
    {
        m_InitStep = std::stod(valuestring);
    }
    else if (keystring == "maxsubstep")
    {
        m_MaxSubStep = std::stoi(valuestring);
    }
    else if (keystring == "minstep")
    {
        m_MinStep = std::stod(valuestring);
    }
    else if (keystring == "maxinc")
    {
        m_MaxInc = std::stod(valuestring);
    }
    else if (keystring == "redmin")
    {
        m_RedMin = std::stod(valuestring);
    }
    else if (keystring == "redmax")
    {
        m_RedMax = std::stod(valuestring);
    }
    else if (keystring == "nsteps")
    {
        m_Nstep = std::stoi(valuestring);
    }
    else if (keystring == "ndigit")
    {
        m_Ndigit = std::stoi(valuestring);
    }
    else if (keystring == "shift_cog")
    {
        m_ShiftCOM = true;
    }
    else if (keystring == "shift_mom")
    {
        m_ShiftMom = true;
    }
    else if (keystring == "restart")
    {
        m_restart = true;
    }
    else if (keystring == "steps")
    {
        m_steps = true;
    }
    else if (keystring == "text_trj")
    {
        m_textTrj = true;
    }
    else if (keystring == "proptype")
    {
        if (valuestring == "unrestricted")
            m_UnResProp = true;
        else if (valuestring == "normal")
            m_UnResProp = false;
        else
            throw std::invalid_argument("unknown propagation type:" + valuestring);
    }
    else if (keystring == "inttype")
    {
        if (valuestring == "bs")
            m_IntType = BS;
        else if (valuestring == "rkqs")
            m_IntType = RkQS;
        else if (valuestring == "rk4fixm")
            m_IntType = Rk4FixM;
        else if (valuestring == "rk4fixt")
            m_IntType = Rk4FixT;
        else
            throw std::invalid_argument("unknown integrator type:" + valuestring);
    }
    else if (keystring == "begin_coords")
    {
        return 1;
    }
    else if (keystring == "end_coords")
    {
        return -1;
    }
    else
        throw std::invalid_argument("unknown keyword: " + keystring);

    return 0;
}

std::string &Configuration::inputFile()
{
    return m_inputFileName;
}

std::string Configuration::baseName()
{
    return m_baseName;
}

double Configuration::tfinal()
{
    return m_tFinal;
}

double Configuration::tout()
{
    return m_tOut;
}

double Configuration::minStep()
{
    return m_MinStep;
}

double Configuration::Eps()
{
    return m_eps;
}

double Configuration::RedMin()
{
    return m_RedMin;
}

double Configuration::RedMax()
{
    return m_RedMax;
}

double Configuration::MaxInc()
{
    return m_MaxInc;
}

std::ofstream &Configuration::outputFile()
{
    return m_outputFile;
}

std::ofstream &Configuration::binTrjFile()
{
    return m_binTrjFile;
}

std::ofstream &Configuration::txtTrjFile()
{
    return m_txtTrjFile;
}

std::ofstream &Configuration::stepsFile()
{
    return m_stepsFile;
}

std::ofstream &Configuration::restartFile()
{
    return m_restartFile;
}

int Configuration::Nobj()
{
    return m_Nobj;
}

int Configuration::MaxSubStep()
{
    return m_MaxSubStep;
}

int Configuration::IncThres()
{
    return m_IncThres;
}

int Configuration::Ndigit()
{
    return m_Ndigit;
}

std::vector<double> Configuration::XX0()
{
    return m_X0;
}

std::vector<double> Configuration::XY0()
{
    return m_Y0;
}

std::vector<double> Configuration::XZ0()
{
    return m_Z0;
}

std::vector<double> Configuration::VX0()
{
    return m_VX0;
}

std::vector<double> Configuration::VY0()
{
    return m_VY0;
}

std::vector<double> Configuration::VZ0()
{
    return m_VZ0;
}

std::vector<double> Configuration::masses()
{
    return m_masses;
}

std::vector<std::string> Configuration::names()
{
    return m_names;
}

bool Configuration::Verbose()
{
    return m_verbose;
}

bool Configuration::ShiftCOM()
{
    return m_ShiftCOM;
}

bool Configuration::ShiftMom()
{
    return m_ShiftMom;
}

bool Configuration::TextTrj()
{
    return m_textTrj;
}

bool Configuration::Restart()
{
    return m_restart;
}

bool Configuration::Steps()
{
    return m_steps;
}

bool Configuration::UnRes()
{
    return m_UnResProp;
}

double Configuration::TotMass()
{
    return m_Mtot;
}

void Configuration::listParas()
{
    m_outputFile << "---------------------\n";
    m_outputFile << "SIMULATION PARAMETERS\n";
    m_outputFile << "---------------------\n\n";

    m_outputFile << std::setprecision(3);
    m_outputFile << "eps         " << std::setw(11) << m_eps        << "\n";
    m_outputFile << "eps_thres   " << std::setw(11) << m_epsThres   << "\n";
    m_outputFile << "tfinal      " << std::setw(11) << m_tFinal     << "\n";
    m_outputFile << "tout        " << std::setw(11) << m_tOut       << "\n";
    m_outputFile << "init_step   " << std::setw(11) << m_InitStep   << "\n";
    m_outputFile << "maxsubstep  " << std::setw(3)  << m_MaxSubStep << "\n";
    m_outputFile << "inc_thres   " << std::setw(3)  << m_IncThres   << "\n";
    m_outputFile << "int_type    " << std::setw(2)  << m_IntType    << "\n";
    m_outputFile << "nstep       " << std::setw(5)  << m_Nstep      << "\n";
    m_outputFile << "min_step    " << std::setw(11) << m_MinStep    << "\n";
    m_outputFile << "maxinc      " << std::setw(11) << m_MaxInc     << "\n";
    m_outputFile << "redmin      " << std::setw(11) << m_RedMin     << "\n";
    m_outputFile << "redmax      " << std::setw(11) << m_RedMax     << "\n";

    m_outputFile << "shift_cog     " << (m_ShiftCOM  ? "yes" : "no") << "\n";
    m_outputFile << "shift_mom     " << (m_ShiftMom  ? "yes" : "no") << "\n";
    m_outputFile << "restart       " << (m_restart   ? "yes" : "no") << "\n";
    m_outputFile << "steps         " << (m_steps     ? "yes" : "no") << "\n";
    m_outputFile << "text_trj      " << (m_textTrj   ? "yes" : "no") << "\n";
    m_outputFile << "prop_type     " << (m_UnResProp ? "unrestricted" : "normal") << "\n";

    m_outputFile << "\n\n\n";
}

IntType Configuration::intType()
{
    return m_IntType;
}

double Configuration::initStep()
{
    return m_InitStep;
}

void Configuration::inc_BSL()
{
    m_N_BS_LargeStep++;
}

void Configuration::inc_BSO()
{
    m_N_BS_OneStep++;
}

void Configuration::inc_BSS()
{
    m_N_BS_SubSteps++;
}

void Configuration::inc_BSE()
{
    m_N_BS_Extrapolate++;
}

} // namespace astronx
