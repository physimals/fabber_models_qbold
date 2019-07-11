/**
 * fwdmodel_qbold_R2p.h - ASE qBOLD curve fitting model
 *
 * Matthew Cherukara, IBME
 *
 * Copyright (C) 2019 University of Oxford  
 */
#pragma once

#include <fabber_core/fwdmodel.h>

#include <newmat.h>

#include <string>
#include <vector>

class R2primeFwdModel : public FwdModel {
public:
    static FwdModel* NewInstance();

    std::string ModelVersion() const;
    void GetOptions(std::vector<OptionSpec> &opts) const;
    std::string GetDescription() const;
    void Initialize(ArgsType &args);
    void GetParameterDefaults(std::vector<Parameter> &params) const;
    
    void Evaluate(const NEWMAT::ColumnVector &params, NEWMAT::ColumnVector &result) const;

protected:

    // Scan Parameters
    double m_tr;
    double m_ti;
    NEWMAT::ColumnVector m_taus;
    NEWMAT::ColumnVector m_tes;

    // Bayesian inference parameters
    double prec_R2p;
    double prec_DBV;
    double prec_CSF;
    double prec_OEF;
    double prec_DF;

    // Lookup starting indices of parameters
    int OEF_index() const
    {
        return (m_infer_oef ? 1 : 0);
    }
    
    int R2p_index() const
    {
        return OEF_index() + (m_infer_r2p ? 1 : 0);
    }

    int DBV_index() const
    {
        return R2p_index() + (m_infer_dbv ? 1 : 0);
    }

    int R2t_index() const
    {
        return DBV_index() + (m_infer_r2t ? 1 : 0);
    }

    int S0_index() const
    {
        return R2t_index() + (m_infer_s0 ? 1 : 0);
    }

    int Hct_index() const
    {
        return S0_index() + (m_infer_hct ? 1 : 0);
    }

    int R2e_index() const
    {
        return Hct_index() + (m_infer_r2e ? 1 : 0);
    }

    int dF_index() const
    {
        return R2e_index() + (m_infer_df ? 1 : 0);
    }

    int lam_index() const
    {
        return dF_index() + (m_infer_lam ? 1 : 0);
    }

    // Which parameters will we infer on
    bool m_infer_oef;
    bool m_infer_r2p;
    bool m_infer_dbv;
    bool m_infer_r2t;
    bool m_infer_s0;
    bool m_infer_hct;
    bool m_infer_r2e;
    bool m_infer_df;
    bool m_infer_lam;

    // Bunch of random booleans for choosing exactly which model we want to run
    bool m_motion_narr;
    bool m_inc_intra;
    bool m_inc_csf;
    bool m_ignore_t1;

private:
    static FactoryRegistration<FwdModelFactory, R2primeFwdModel> registration;
};
