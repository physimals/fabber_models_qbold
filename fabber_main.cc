/*  FABBER - Fast ASL and BOLD Bayesian Estimation Routine

 Adrian Groves and Michael Chappell, FMRIB Image Analysis & IBME QuBIc groups

 Copyright (C) 2007-2015 University of Oxford  */

/*  CCOPYRIGHT */

#include "fabber_core/fabber_core.h"

#include "fwdmodel_qbold.h"

int main(int argc, char **argv)
{
    // Add the Qbold models - these will autoregister at this point
    R2primeFwdModel::NewInstance();
    return execute(argc, argv);
}
