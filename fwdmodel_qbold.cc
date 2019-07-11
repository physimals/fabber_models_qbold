/**
 * fwdmodel_qbold.cc - ASE qBOLD curve fitting model
 *
 * Matthew Cherukara, IBME
 *
 * Copyright (C) 2019 University of Oxford  
 */

#include "fwdmodel_qbold.h"

#include <fabber_core/fwdmodel.h>
#include <fabber_core/priors.h>

#include <newmat.h>

#include <vector>
#include <string>
#include <complex>

using namespace std;
using namespace NEWMAT;

FactoryRegistration<FwdModelFactory, R2primeFwdModel>
    R2primeFwdModel::registration("qboldR2p");

FwdModel *R2primeFwdModel::NewInstance() // unchanged
{
    return new R2primeFwdModel();
}

string R2primeFwdModel::GetDescription() const 
{
    return "ASE qBOLD model R2-prime version";
}

static OptionSpec OPTIONS[] = {
    { "inferOEF", OPT_BOOL, "", OPT_NONREQ, "" },
    { "inferR2p", OPT_BOOL, "", OPT_NONREQ, "" },
    { "inferDBV", OPT_BOOL, "", OPT_NONREQ, "" },
    { "inferR2t", OPT_BOOL, "", OPT_NONREQ, "" },
    { "inferS0", OPT_BOOL, "", OPT_NONREQ, "" },
    { "inferHct", OPT_BOOL, "", OPT_NONREQ, "" },
    { "inferR2e", OPT_BOOL, "", OPT_NONREQ, "" },
    { "inferdF", OPT_BOOL, "", OPT_NONREQ, "" },
    { "inferlam", OPT_BOOL, "", OPT_NONREQ, "" },
    { "include_intra", OPT_BOOL, "", OPT_NONREQ, "" },
    { "include_csf", OPT_BOOL, "", OPT_NONREQ, "" },
    { "m_ignore_t1", OPT_BOOL, "", OPT_NONREQ, "" },
    { "motional_narrowing", OPT_BOOL, "", OPT_NONREQ, "" },
    { "tau<n>", OPT_FLOAT, "", OPT_REQ, "" },
    { "TR", OPT_FLOAT, "", OPT_REQ, "3.0" },
    { "TI", OPT_FLOAT, "", OPT_REQ, "0.0" },
    { "" },
};

void R2primeFwdModel::GetOptions(vector<OptionSpec> &opts) const
{
    for (int i = 0; OPTIONS[i].name != ""; i++)
    {
        opts.push_back(OPTIONS[i]);
    }
}

string R2primeFwdModel::ModelVersion() const
{
    string version = "fwdmodel_qbold.cc";
#ifdef GIT_SHA1
    version += string(" Revision ") + GIT_SHA1;
#endif
#ifdef GIT_DATE
    version += string(" Last commit ") + GIT_DATE;
#endif
    return version;
}

void R2primeFwdModel::Initialize(FabberRunData &rundata)
{
    m_infer_oef = rundata.ReadBool("inferOEF");
    m_infer_r2p = rundata.ReadBool("inferR2p");
    m_infer_dbv = rundata.ReadBool("inferDBV");
    m_infer_r2t = rundata.ReadBool("inferR2t");
    m_infer_s0  = rundata.ReadBool("inferS0");
    m_infer_hct = rundata.ReadBool("inferHct");
    m_infer_r2e = rundata.ReadBool("inferR2e");
    m_infer_df  = rundata.ReadBool("inferdF");
    m_infer_lam = rundata.ReadBool("inferlam");

    m_inc_intra = rundata.ReadBool("include_intra");
    m_inc_csf = rundata.ReadBool("include_csf");
    m_ignore_t1 = rundata.ReadBool("m_ignore_t1");

    m_motion_narr = rundata.ReadBool("motional_narrowing");

    // since we can't do both, OEF will take precidence over R2p
    if (m_infer_oef)
    {
        m_infer_r2p = false;
    }

    if (m_motion_narr)
    {
        m_inc_intra = true;
    }

    if (m_infer_lam || m_infer_df )
    {
        m_inc_csf = true;
    }

    // temporary holders for input values
    string tau_temp;
    string TE_temp; 

    // allow for manual entry of prior precisions
    prec_R2p = rundata.GetDoubleDefault("precR2p", 1e-3);
    prec_DBV = rundata.GetDoubleDefault("precDBV", 1e-1);
    prec_CSF = rundata.GetDoubleDefault("precCSF", 1e-1);
    prec_OEF = rundata.GetDoubleDefault("precOEF", 1e-1);
    prec_DF  = rundata.GetDoubleDefault("precD", 1e-4);

    // First read tau values, since these will always be specified
    vector<double> tausv = rundata.GetDoubleList("tau");
    m_taus.ReSize(tausv.size());
    for (unsigned int i=0; i<tausv.size(); i++) m_taus(i+1) = tausv[i];

    // Then read TE values. There might be a sequence or just one
    vector<double> tes = rundata.GetDoubleList("TE");
    m_tes.ReSize(tausv.size());
    if (tes.size() == 1) m_tes = tes[0];
    else
    {
        for (unsigned int i=0; i<tes.size(); i++)
        {
            m_tes(i+1) = tes[i];
        }
    }

    // read TR and TI
    m_tr = rundata.GetDoubleDefault("TR", 3.0);
    m_ti = rundata.GetDoubleDefault("TI", 0.0);

    // add information to the log
    LOG << "Inference using development model" << endl;
    LOG << "Using TR = " << m_tr << "s, and TI = " << m_ti << "s" << endl;
    for (int ii = 1; ii <= m_taus.Nrows(); ii++)
    {
        LOG << "    TE(" << ii << ") = " << m_tes(ii) << "    tau(" << ii << ") = " << m_taus(ii) << endl;
    }
    if (m_inc_intra)
    {
        if (m_motion_narr)
        {
            LOG << "Using two-compartment model with motion narrowing intravascular signal " << endl;
        }
        else
        {
            LOG << "Using two-compartment model with static (powder) intravascular signal " << endl;
        }
    }
    else
    {
        LOG << "Using single-compartment model (ignoring intravascular signal) " << endl;
    }
    if (m_infer_oef)
    {
        LOG << "Inferring on OEF " << endl;
    }
    if (m_infer_r2p)
    {
        LOG << "Inferring on R2p " << endl;
    }
    if (m_infer_dbv)
    {
        LOG << "Inferring on DBV " << endl;
    }
    if (m_infer_r2t)
    {
        LOG << "Inferring on R2/T2 of tissue " << endl;
    }
    if (m_infer_s0)
    {
        LOG << "Inferring on scaling parameter S0 " << endl;
    }
    if (m_infer_hct)
    {
        LOG << "Inferring on fractional hematocrit " << endl;
    }
    if (m_infer_r2e)
    {
        LOG << "Inferring on R2 of CSF " << endl;
    }
    if (m_infer_df)
    {
        LOG << "Inferring on CSF frequency shift dF " << endl;
    }
    if (m_infer_lam)
    {
        LOG << "Inferring on CSF volume fraction lambda " << endl;
    }
}

void R2primeFwdModel::GetParameterDefaults(std::vector<Parameter> &params) const
{
    params.clear();

    int p=0;
    if (m_infer_oef) params.push_back(Parameter(p++, "OEF", DistParams(0.4, 10), DistParams(0.4, 10), PRIOR_NORMAL, TRANSFORM_IDENTITY()));
    if (m_infer_r2p) params.push_back(Parameter(p++, "R2p", DistParams(4.0, 1e3), DistParams(4.0, 1e2), PRIOR_NORMAL, TRANSFORM_IDENTITY()));
    if (m_infer_dbv) params.push_back(Parameter(p++, "DBV", DistParams(0.03, 10), DistParams(0.03, 10), PRIOR_NORMAL, TRANSFORM_IDENTITY()));
    if (m_infer_r2t) params.push_back(Parameter(p++, "R2t", DistParams(1/0.087, 1e2), DistParams(1/0.087, 1e2), PRIOR_NORMAL, TRANSFORM_IDENTITY()));
    if (m_infer_s0) params.push_back(Parameter(p++, "S0", DistParams(500.0, 1e6), DistParams(500.0, 100), PRIOR_NORMAL, TRANSFORM_IDENTITY()));
    if (m_infer_hct) params.push_back(Parameter(p++, "Hct", DistParams(0.40, 1e-3), DistParams(0.4, 1e-3), PRIOR_NORMAL, TRANSFORM_IDENTITY()));
    if (m_infer_r2e) params.push_back(Parameter(p++, "R2e", DistParams(0.5, 1e2), DistParams(0.5, 1e2), PRIOR_NORMAL, TRANSFORM_IDENTITY()));
    if (m_infer_df) params.push_back(Parameter(p++, "dF", DistParams(5.0, 1e4), DistParams(5.0, 100), PRIOR_NORMAL, TRANSFORM_IDENTITY()));
    if (m_infer_lam) params.push_back(Parameter(p++, "VC", DistParams(0.05, 10), DistParams(0.05, 10), PRIOR_NORMAL, TRANSFORM_IDENTITY()));
}

void R2primeFwdModel::Evaluate(const ColumnVector &params, ColumnVector &result) const
{
    // Check we have been given the right number of parameters
    assert(params.Nrows() == NumParams());
    result.ReSize(data.Nrows());

    // Calculated parameters
    double Ss;  // static dephasing signal (will be used to make St and Se)
    double St;  // tissue signal
    double Sb;  // blood signal
    double Se;  // extracellular signal
    complex<double> Sec; // complex version of the extracellular signal (may be unnecessary)
    complex<double> Sbc; // complex version of the blood signal (for powder model)
    complex<double> i(0,1);

    // Derived parameters
    double dw;          // characteristic time (protons in water)
    double tc;
    double lam0;        // apparent lambda (opposite way round from literature!)
    double mt;
    double me;
    double mb;

    // Default values for parameters - some of these may be inferred
    double OEF = 0.3;
    double R2p = 2.5;
    double DBV = 0.03;
    double R2t = 11.5;
    double S0 = 265.0;
    double Hct = 0.40;
    double R2e = 4.0;
    double dF = 5.00;
    double lam = 0.0;
    double CBV;

    // Assign values to parameters
    if (m_infer_dbv)
    {
        DBV = (params(DBV_index()));
        // this bit makes sure the value for DBV isn't ridiculous
        if (DBV < 0.0001)
        {
            DBV = 0.0001;
        }
        else if (DBV > 0.5)
        {
            DBV = 0.5;
        }
    }
    
    if (m_infer_r2t)
    {
        R2t = (params(R2t_index()));
    }
    
    if (m_infer_s0)
    {
        S0 = (params(S0_index()));
    }
    
    if (m_infer_hct)
    {
        Hct = (params(Hct_index()));
    }
    
    if (m_infer_r2e)
    {
        R2e = (params(R2e_index()));
    }
    
    if (m_infer_df)
    {
        dF = (params(dF_index()));
    }
    
    if (m_infer_lam)
    {
        lam = (params(lam_index()));
    }

    // This one is a little bit different
    if (m_infer_oef)
    {
        OEF = (params(OEF_index()));
        dw = 887.4082*Hct*OEF;
        R2p = dw*DBV;
    }
    else if (m_infer_r2p)
    {
        R2p = (params(R2p_index()));
        if (R2p < 0.01)
        {
            R2p = 0.01;
        }

        OEF = R2p/(887.4082*DBV*Hct);
        dw = 887.4082*Hct*OEF;
    }

    // Calculate tc and threshold it if necessary
    tc = 1.76/dw;

    // Evaluate blood relaxation rates (for linear model)
    double R2b  = ( 4.5 + (16.4*Hct)) + ( ((165.2*Hct) + 55.7)*pow(OEF,2.0) );
    double R2bp = (10.2 - ( 1.5*Hct)) + ( ((136.9*Hct) - 13.9)*pow(OEF,2.0) );
    
    // Simulated data doesn't have any T1 contrast
    if (m_ignore_t1)
    {
        lam0 = lam;
        CBV = DBV;
    }
    else
    {
        // Here are some more constants we will need
        double T1t = 1.20;
        double T1e = 3.87;
        double T1b = 1.58;

        double nt = 0.723;
        double ne = 1.000;
        double nb = 0.775;

        // Calculate magnetizations
        mt = 1.0 - ( ( 2.0 - exp(-(m_tr-m_ti)/T1t) ) * exp(-m_ti/T1t) );
        mb = 1.0 - ( ( 2.0 - exp(-(m_tr-m_ti)/T1b) ) * exp(-m_ti/T1b) );
        me = 1.0 - ( ( 2.0 - exp(-(m_tr-m_ti)/T1e) ) * exp(-m_ti/T1e) );

        // Calculate tissue compartment weightings
        lam0 = (ne*me*lam) / ( (nt*mt*(1-lam)) + (ne*me*lam) );
        CBV = nb*mb*(1-lam0)*DBV;
    }
    
    // loop through taus
    result.ReSize(m_taus.Nrows());
    for (int ii = 1; ii <= m_taus.Nrows(); ii++)
    {
        double tau = m_taus(ii);
        double TE = m_tes(ii);

        // Calculate tissue signal
        if (tau < -tc)
        {
            Ss = exp(DBV + (R2p*tau));          // SDR model
        }
        else if (tau > tc)
        {
            Ss = exp(DBV - (R2p*tau));          // SDR model
        }
        else
        {
            Ss = exp(-0.3*pow(R2p*tau,2.0)/DBV);          // SDR model
        }

        // Add T2 effect to tissue compartment
        St = Ss*exp(-R2t*TE);

        // Calculate intravascular signal
        if (m_motion_narr)
        {
            double td   = 0.0045067;       // (based on rc=2.6 um and D=1.5 um^2 / ms)
            double gm   = 2.67513e8;
            double dChi = ((0.27*OEF) + 0.14)*1e-6;
            double G0   = (4/45)*Hct*(1-Hct)*pow((dChi*3.0),2.0);
            double kk   = 0.5*pow(gm,2.0)*G0*pow(td,2.0);
            R2b = 5.291;    // fixed value (Berman, 2017) 

            // motion narrowing model
            Sb = exp(-kk* ( (TE/td) + pow((0.25 + (TE/td)),0.5) + 1.5 - 
                            (2*pow((0.25 + (pow((TE+tau),2.0)/td) ),0.5)) - 
                            (2*pow((0.25 + (pow((TE-tau),2.0)/td) ),0.5)) ) );

            // T2 effect 
            Sb *= exp(-R2b*TE); 

        }
        else if (m_inc_intra)
        {
            // linear model
            //Sb = exp(-R2b*TE)*exp(-R2bp*abs(tau));

            // powder model

            // threshold constant (similar to tc but different?)
            double pp = 1.5*dw*tau;

            // signum of tau
            double st = (tau > 0) - (tau < 0);

            if (abs(pp) > 1)
            {
                // large pp
                Sbc = 0.5*sqrt(M_PI/abs(pp))*exp(i*((pp/3.0)-(st*M_PI/4.0)));
            }
            else
            {
                // small pp
                Sbc = 1.0 - ((2.0/45.0)*pow(pp,2.0)) + ((8.0*i*pow(pp,3.0))/2835.0);
            }

            Sb = real(Sbc);

            // T2 effect
            Sb *= exp(-R2b*TE);
        }
        else
        {
            Sb = 0.0;
        }

        // Calculate CSF signal
        if (m_inc_csf)
        {
            Sec = Ss*exp(-R2e*TE)*exp(-2.0*i*M_PI*dF*abs(tau));
            Se = abs(Sec);
        }
        else
        {
            Se = 0.0;
        }

        // Add up the compartments
        result(ii) = S0*(((1-CBV-lam0)*St) + (CBV*Sb) + (lam0*Se)); 
    }
}
