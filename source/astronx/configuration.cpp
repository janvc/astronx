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

        if (currentLine.find("end_coords") == std::string::npos)    // proceed if this is NOT the end of the coordinates
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

    // open the output file
    int lastDot = m_inputFileName.find_last_of(".");
    if (lastDot == std::string::npos)
        m_outputFileName = m_inputFileName + ".out";

    m_outputFileName = m_inputFileName.substr(0, lastDot) + ".out";

    std::cout << m_outputFileName;

    m_outputFile = std::ofstream(m_outputFileName);

    return 0;
}

void Configuration::setDefaults()
{
    m_MaxSubStep = 12;
    m_IncThres   =  8;
    m_Ndigit     = 10;
    m_Nstep      =  5;
    m_tFinal     =  0.0;
    m_tOut       =  0.0;
    m_InitStep   =  0.0;
    m_eps        =  1.0e-6;
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
    else if (keystring == "shiftcom")
    {
        m_ShiftCOM = true;
    }
    else if (keystring == "shiftmom")
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
    else if (keystring == "texttrj")
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

double Configuration::tfinal()
{
    return m_tFinal;
}

std::ofstream &Configuration::outputFile()
{
    return m_outputFile;
}

} // namespace astronx
