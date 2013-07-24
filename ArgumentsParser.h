//=============================================================================
//
//  Copyright (c) Kitware, Inc.
//  All rights reserved.
//  See LICENSE.txt for details.
//
//  This software is distributed WITHOUT ANY WARRANTY; without even
//  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
//  PURPOSE.  See the above copyright notice for more information.
//
//  Copyright 2012 Sandia Corporation.
//  Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
//  the U.S. Government retains certain rights in this software.
//
//=============================================================================
#ifndef __argumentsParser_h
#define __argumentsParser_h

#include <string>

namespace dax { namespace testing {

class ArgumentsParser
{
public:
  ArgumentsParser();
  virtual ~ArgumentsParser();

  bool parseArguments(int argc, char* argv[]);

  enum PipelineMode
    {
    THRESHOLD = 1,
    MARCHING_CUBES = 2

    };
  PipelineMode pipeline() const
    { return this->Pipeline; }

  std::string file() const
    { return this->File; }

  int cores() const
    { return this->Cores; }

  double resampleRate() const
  { return this->ResampleRate; }


private:
  std::string File;
  PipelineMode Pipeline;
  int Cores;
  double ResampleRate;


};

}}
#endif
