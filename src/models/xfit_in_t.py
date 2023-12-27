#! /usr/bin/python

#
# Description:
# ================================================================
# Time-stamp: "2023-12-27 05:31:47 trottar"
# ================================================================
#
# Author:  Richard L. Trotta III <trotta@cua.edu>
#
# Copyright (c) trottar
#
import ROOT
from ROOT import TFile, TNtuple, TText
import math

################################################################################################################################################
'''
ltsep package import and pathing definitions
'''

# Import package for cuts
from ltsep import Root
# Import package for progress bar
from ltsep import Misc

lt=Root(os.path.realpath(__file__),"Plot_Prod")

# Add this to all files for more dynamic pathing
USER=lt.USER # Grab user info for file finding
HOST=lt.HOST
REPLAYPATH=lt.REPLAYPATH
UTILPATH=lt.UTILPATH
LTANAPATH=lt.LTANAPATH
ANATYPE=lt.ANATYPE
OUTPATH=lt.OUTPATH
CACHEPATH=lt.CACHEPATH

################################################################################################################################################

'''
1) ?????????grab file src/kaon/xsects/x_sep.pl_{q_set}_{eps} -> sigL,dsigL,sigT,dsigT,sigTT, dsigTT,sigLT,dsigLT,q2,tm
   Bill grabs the following from this file==> t:t_e:l:l_e:lt:lt_e:tt:tt_e:chi:u:u_min:t:w:q2
   Is he actually grabbing src/kaon/xsects/x_unsep.pl_{q_set}_{eps} -> x_real,dx_real,x_mod,eps_mod, th_mod*180./3.14159,phi*180./3.14159,tm,w,q2

   * Need the values t:t_e:l:l_e:lt:lt_e:tt:tt_e:chi:u:u_min:t:w:q2 
      May need to adjust what is written to x_sep.pl_* but what is the error?
      Perhaps, this is actually the separated values NOT from model. In this case, Rosenbluth sep script is needed.
   ** Such a script would involve Fitting results of L,T,LT,TT for each t-bin??

2) Read in param files from src/kaon/parameters/par.pl_{q_set}.dat

3) Read in src/kaon/averages/avek.{q_set}.dat
   Bill grabs the following from this file==> iss >> w >> w_e >> q2 >> q2_e >> tt >>  tt_e >> thetacm >> it 

***COMPLETE REST OF SCRIPT THEN WORRY ABOUT ABOVE***

4)  t0  = prv_par_vec[0]; 
    t1  = prv_par_vec[1]; 
    t2  = prv_par_vec[2]; 
    t3  = prv_par_vec[3]; 
    l0  = prv_par_vec[4]; 
    l1  = prv_par_vec[5]; 
    l2  = prv_par_vec[6]; 
    l3  = prv_par_vec[7]; 
    lt0 = prv_par_vec[8]; 
    lt1 = prv_par_vec[9]; 
    lt2 = prv_par_vec[10]; 
    lt3 = prv_par_vec[11]; 
    tt0 = prv_par_vec[12]; 
    tt1 = prv_par_vec[13]; 
    tt2 = prv_par_vec[14]; 
    tt3 = prv_par_vec[15]; 


5) Define g = 1 /pow(pow(w,2)-pow(m_p,2),2);

6) Define plots...

    TGraph* g_sigt_prv = new TGraph();
    TGraph* g_sigl_prv = new TGraph();
    TGraph* g_siglt_prv = new TGraph();
    TGraph* g_sigtt_prv = new TGraph();

    TGraphErrors* g_sigt_fit  = new TGraphErrors();
    TGraphErrors* g_sigl_fit  = new TGraphErrors();
    TGraphErrors* g_siglt_fit = new TGraphErrors();
    TGraphErrors* g_sigtt_fit = new TGraphErrors();

    TGraph* g_sigt_fit_tot  = new TGraph();
    TGraph* g_sigl_fit_tot  = new TGraph();
    TGraph* g_siglt_fit_tot = new TGraph();
    TGraph* g_sigtt_fit_tot = new TGraph();

    TCanvas* c1 =  new TCanvas("c1", "c1", 800, 800); 
    c1->Divide(2,2);

    TCanvas* c2 =  new TCanvas("c2", "c2", 800, 800); 
    c2->Divide(2,2);

    n1->Draw("u:u_min", "", "goff");

   // Create a TGraph to map u_min values to corresponding u values
   TGraph* u_umin_map = new TGraph(n1->GetSelectedRows(), n1->GetV1(), n1->GetV2());

   // Extract u and u_min lists from the mapping
   Double_t* u_list = u_umin_map->GetX();
   Double_t* u_min_list = u_umin_map->GetY();

'''
'''

1) Iterate for transverse params

    c1->cd(1)->SetLeftMargin(0.12);
    n1->Draw("t:u:t_e", "", "goff");
    TGraphErrors* g_sigt = new TGraphErrors(n1->GetSelectedRows(), n1->GetV2(), n1->GetV1(), 0, n1->GetV3());

    TF1* f_sigT_pre = new TF1("sig_T_pre", fun_Sig_T, 0, 0.5, 2);
    f_sigT_pre->SetParameters(t0, t1);

    for(int i =0; i < w_vec.size(); i++) {

            Float_t sigt_X_pre;

            Float_t q2_term;
            q2_term = t2 / logq2_vec[i] + t3 * g_sigt->GetX()[i] / logq2_vec[i];
            q2_dep = sqrt(q2_vec[i]);

            sigt_X_pre = (f_sigT_pre->Eval(g_sigt->GetX()[i]) / q2_dep + q2_term) * g_vec[i] ; 
            g_sigt_prv->SetPoint(i, g_sigt->GetX()[i], sigt_X_pre);

            Float_t sigt_X_fit, sigt_X_fit_err;
            sigt_X_fit     = ((g_sigt->GetY()[i]) / g_vec[i]   - q2_term) * q2_dep;
            sigt_X_fit_err = g_sigt->GetEY()[i] / g_vec[i] * q2_dep;

            g_sigt_fit->SetPoint(i, g_sigt->GetX()[i], sigt_X_fit);
            g_sigt_fit->SetPointError(i, 0, sigt_X_fit_err);
}

2) 	g_max = g_sigt->GetYaxis()->GetXmax();
    gp_max = TMath::MaxElement(g_sigt_prv->GetN(), g_sigt_prv->GetY());
    g_min = g_sigt->GetYaxis()->GetXmin();
    gp_min = TMath::MinElement(g_sigt_prv->GetN(), g_sigt_prv->GetY());

    Double_t difff = (g_max - g_min)/5;

    if (g_max < gp_max) {
            g_sigt->SetMaximum(gp_max + difff);
    } 

    if (g_min > gp_min) {
            g_sigt->SetMinimum(gp_min - difff);
    } 


3) 	g_sigt->SetTitle("Sig T");

    g_sigt->SetMarkerStyle(5);
    g_sigt->Draw("AP");

    g_sigt->GetXaxis()->SetTitle("#it{-u} [GeV^{2}]");
    g_sigt->GetXaxis()->CenterTitle();
    g_sigt->GetYaxis()->SetTitle("#left(#frac{#it{d#sigma}}{#it{du}}#right)_{T} [#mub/GeV^{2}]");
    g_sigt->GetYaxis()->SetTitleOffset(1.5);
    g_sigt->GetYaxis()->SetTitleSize(0.035);
    g_sigt->GetYaxis()->CenterTitle();

    g_sigt_prv->SetMarkerColor(4);
    g_sigt_prv->SetMarkerStyle(21);
    g_sigt_prv->Draw("P");

    c2->cd(1);
    g_sigt_fit->SetTitle("Sigma T Model Fit");
    g_sigt_fit->Draw("A*");

    TF1* f_sigT = new TF1("sig_T", fun_Sig_T, 0, 0.5, 2); 
    f_sigT->SetParameters(t0, t1);

4) 	TFitResultPtr fit_t_result	= g_sigt_fit->Fit(f_sigT, "S");

    for(int i =0; i < w_vec.size(); i++) {

            Float_t sigt_X;

            Float_t q2_term;
            q2_term = t2 / logq2_vec[i] + t3 * g_sigt->GetX()[i] / logq2_vec[i];

            q2_dep = sqrt(q2_vec[i]);

            sigt_X = (f_sigT->Eval(g_sigt->GetX()[i]) / q2_dep + q2_term) * g_vec[i] ; 

            cout << f_sigT->Eval(g_sigt->GetX()[i]) << "      " << sigt_X << endl;

            g_sigt_fit_tot->SetPoint(i, g_sigt->GetX()[i], sigt_X);

    }

5) fit_status->DrawTextNDC(0.35, 0.8, " Fit Status: " + gMinuit->fCstatu);

    c1->cd(1);

    g_sigt_fit_tot->SetMarkerStyle(26);
    g_sigt_fit_tot->SetMarkerColor(2);
    g_sigt_fit_tot->SetLineColor(2);
    g_sigt_fit_tot->Draw("LP");

6) t0 = f_sigT->GetParameter(0);
    t1 = f_sigT->GetParameter(1);

    par_vec.push_back(t0);
    par_vec.push_back(t1);
    par_vec.push_back(t2);
    par_vec.push_back(t3);

    par_err_vec.push_back(f_sigT->GetParError(0));
    par_err_vec.push_back(f_sigT->GetParError(1));
    par_err_vec.push_back(0.0);
    par_err_vec.push_back(0.0);

    par_chi2_vec.push_back(f_sigT->GetChisquare());
    par_chi2_vec.push_back(f_sigT->GetChisquare());
    par_chi2_vec.push_back(f_sigT->GetChisquare());
    par_chi2_vec.push_back(f_sigT->GetChisquare());

'''

'''

1) Iterate for longitudinal params

    c1->cd(2)->SetLeftMargin(0.12);
    n1->Draw("l:u:l_e", "", "goff");

    TF1* f_sigL_pre = new TF1("sig_L", fun_Sig_L, 0, 0.5, 2);
    f_sigL_pre->SetParameters(l0, l1);

    TGraphErrors* g_sigl = new TGraphErrors(n1->GetSelectedRows(), n1->GetV2(), n1->GetV1(), 0, n1->GetV3());

    for(int i =0; i < w_vec.size(); i++) {

            Float_t sigl_X_pre;

            Float_t q2_term;
            q2_term = l2 / q2_vec[i] + l3 * g_sigl->GetX()[i] / q2_vec[i];
            q2_dep = q2_vec[i]*q2_vec[i];

            sigl_X_pre = (f_sigL_pre->Eval(g_sigl->GetX()[i])/q2_dep + q2_term) * g_vec[i]; 
            g_sigl_prv->SetPoint(i, g_sigl->GetX()[i], sigl_X_pre);

            Float_t sigl_X_fit, sigl_X_fit_err;
            sigl_X_fit     = (g_sigl->GetY()[i] / g_vec[i]  - q2_term) * q2_dep;
            sigl_X_fit_err = g_sigl->GetEY()[i] / g_vec[i] * q2_dep;

            cout << "aa  " << sigl_X_fit  << "  " << sigl_X_fit_err << endl;

            g_sigl_fit->SetPoint(i, g_sigl->GetX()[i], sigl_X_fit);
            g_sigl_fit->SetPointError(i, 0, sigl_X_fit_err);

            cout << l0 << "  " << l1 << "  " << q2_dep << "  " << g_sigt->GetX()[i] << "  "<< w_vec[i] << "  " << g_vec[i] << "  " <<  sigl_X_pre << endl;

    }


2)  g_max = g_sigl->GetYaxis()->GetXmax();
    gp_max = TMath::MaxElement(g_sigl_prv->GetN(), g_sigl_prv->GetY());

    g_min = g_sigl->GetYaxis()->GetXmin();
    gp_min = TMath::MinElement(g_sigl_prv->GetN(), g_sigl_prv->GetY());

    difff = (g_max - g_min)/5;

    if (g_max < gp_max) {
            g_sigl->SetMaximum(gp_max + difff);
    } 

    if (g_min > gp_min) {
            g_sigl->SetMinimum(gp_min - difff);
    } 

3) g_sigl->SetTitle("Sig L");

    g_sigl->SetMarkerStyle(5);
    g_sigl->Draw("AP");

    g_sigl->GetXaxis()->SetTitle("#it{-u} [GeV^{2}]");
    g_sigl->GetXaxis()->CenterTitle();
    g_sigl->GetYaxis()->SetTitle("#left(#frac{#it{d#sigma}}{#it{du}}#right)_{L} [#mub/GeV^{2}]");
    g_sigl->GetYaxis()->SetTitleOffset(1.5);
    g_sigl->GetYaxis()->SetTitleSize(0.035);
    g_sigl->GetYaxis()->CenterTitle();

    g_sigl_prv->SetMarkerColor(4);
    g_sigl_prv->SetMarkerStyle(21);
    g_sigl_prv->Draw("P");

    c2->cd(2);
    g_sigl_fit->SetTitle("Sigma L Model Fit");
    g_sigl_fit->Draw("A*");

    TF1* f_sigL = new TF1("sig_L", fun_Sig_L, 0, 0.5, 2);
    f_sigL->SetParameters(l0, l1);

4) g_sigl_fit->Fit(f_sigL);

    for(int i =0; i < w_vec.size(); i++) {

            Float_t sigl_X;
            Float_t q2_term;
            q2_term = l2 * q2_vec[i] + l3 * g_sigl->GetX()[i] * q2_vec[i];

            q2_dep = q2_vec[i]*q2_vec[i];

            sigl_X = (f_sigL->Eval(g_sigl->GetX()[i]) / q2_dep + q2_term) * g_vec[i]; 
            g_sigl_fit_tot->SetPoint(i, g_sigl->GetX()[i], sigl_X);

    }

5) fit_status->DrawTextNDC(0.35, 0.8, " Fit Status: " + gMinuit->fCstatu);

    c1->cd(2);

    g_sigl_fit_tot->SetMarkerStyle(26);
    g_sigl_fit_tot->SetMarkerColor(2);
    g_sigl_fit_tot->SetLineColor(2);
    g_sigl_fit_tot->Draw("LP");

6) l0 = f_sigL->GetParameter(0);
    l1 = f_sigL->GetParameter(1);

    par_vec.push_back(l0);
    par_vec.push_back(l1);
    par_vec.push_back(l2);
    par_vec.push_back(l3);

    par_err_vec.push_back(f_sigL->GetParError(0)); 
    par_err_vec.push_back(f_sigL->GetParError(1)); 
    par_err_vec.push_back(0);
    par_err_vec.push_back(0);

    par_chi2_vec.push_back(f_sigL->GetChisquare());
    par_chi2_vec.push_back(f_sigL->GetChisquare());
    par_chi2_vec.push_back(f_sigL->GetChisquare());
    par_chi2_vec.push_back(f_sigL->GetChisquare());

'''

'''

1) Iterate for LT params

    c1->cd(3)->SetLeftMargin(0.12);
    n1->Draw("lt:u:lt_e", "", "goff");

    TF1* f_sigLT_pre = new TF1("sig_LT", fun_Sig_LT, 0, 0.5, 2);
    f_sigLT_pre->SetParameters(lt0, lt1);

    TGraphErrors* g_siglt = new TGraphErrors(n1->GetSelectedRows(), n1->GetV2(), n1->GetV1(), 0, n1->GetV3());

    for(int i =0; i < w_vec.size(); i++) {

            Float_t siglt_X_pre;
            Float_t q2_term;
            q2_term = lt2 / logq2_vec[i] + lt3 * g_siglt->GetX()[i] / logq2_vec[i];

            q2_dep = q2_vec[i];

            siglt_X_pre =  (f_sigLT_pre->Eval(g_siglt->GetX()[i]) / q2_dep + q2_term) * g_vec[i] * sin(th_vec[i] * pi / 180); 
            g_siglt_prv->SetPoint(i, g_sigl->GetX()[i], siglt_X_pre);

            Float_t siglt_X_fit, siglt_X_fit_err;

            cout << th_vec[i] << endl;

            cout << "!!!!!!!!!!!!!!!!!!     " << q2_term << endl;

            if( th_vec[i]  == 180) { 

                    siglt_X_fit     = 0.0;
                    siglt_X_fit_err = 1.0;

            } else {

                    siglt_X_fit     = (g_siglt->GetY()[i]  / g_vec[i] / sin(th_vec[i] * pi / 180)  - q2_term) * q2_dep;
                    siglt_X_fit_err = g_siglt->GetEY()[i] / g_vec[i] / sin(th_vec[i] * pi / 180) * q2_dep;

            }

            g_siglt_fit->SetPoint(i, g_siglt->GetX()[i], siglt_X_fit);
            g_siglt_fit->SetPointError(i, 0, siglt_X_fit_err);

    }

2) 	g_max = g_siglt->GetYaxis()->GetXmax();
    gp_max = TMath::MaxElement(g_siglt_prv->GetN(), g_siglt_prv->GetY());

    g_min = g_siglt->GetYaxis()->GetXmin();
    gp_min = TMath::MinElement(g_siglt_prv->GetN(), g_siglt_prv->GetY());

    difff = (g_max - g_min)/5;

    if (g_max < gp_max) {
            g_siglt->SetMaximum(gp_max + difff);
    } 

    if (g_min > gp_min) {
            g_siglt->SetMinimum(gp_min - difff);
    } 


3) 	g_siglt->SetTitle("Sig LT");

    g_siglt->SetMarkerStyle(5);
    g_siglt->Draw("AP");

    g_siglt->GetXaxis()->SetTitle("#it{-u} [GeV^{2}]");
    g_siglt->GetXaxis()->CenterTitle();
    g_siglt->GetYaxis()->SetTitle("#left(#frac{#it{d#sigma}}{#it{du}}#right)_{LT} [#mub/GeV^{2}]");
    g_siglt->GetYaxis()->SetTitleOffset(1.5);
    g_siglt->GetYaxis()->SetTitleSize(0.035);
    g_siglt->GetYaxis()->CenterTitle();

    g_siglt_prv->SetMarkerColor(4);
    g_siglt_prv->SetMarkerStyle(21);
    g_siglt_prv->Draw("P");

    c2->cd(3);
    g_siglt_fit->SetTitle("Sigma LT Model Fit");
    g_siglt_fit->Draw("A*");

    TF1* f_sigLT = new TF1("sig_LT", fun_Sig_LT, 0, 0.5, 2); 
    f_sigLT->SetParameters(lt0, lt1);

4) 	g_siglt_fit->Fit(f_sigLT);

    for(int i =0; i < w_vec.size(); i++) {

            Float_t siglt_X;
            Float_t q2_term;
            q2_term = lt2 / logq2_vec[i] + lt3 * g_siglt->GetX()[i] / logq2_vec[i];

            q2_dep = q2_vec[i];

            siglt_X = (f_sigLT->Eval(g_siglt->GetX()[i]) / q2_dep + q2_term) * g_vec[i]   * sin(th_vec[i] * pi / 180); 
            g_siglt_fit_tot->SetPoint(i, g_siglt->GetX()[i], siglt_X);

    }

5) 	fit_status->DrawTextNDC(0.35, 0.8, " Fit Status: " + gMinuit->fCstatu);

    c1->cd(3);

    g_siglt_fit_tot->SetMarkerStyle(26);
    g_siglt_fit_tot->SetMarkerColor(2);
    g_siglt_fit_tot->SetLineColor(2);
    g_siglt_fit_tot->Draw("LP");

6) 	lt0 = f_sigLT->GetParameter(0);
    lt1 = f_sigLT->GetParameter(1);

    par_vec.push_back(lt0);
    par_vec.push_back(lt1);
    par_vec.push_back(lt2);
    par_vec.push_back(lt3);


    par_err_vec.push_back(f_sigLT->GetParError(0));
    par_err_vec.push_back(f_sigLT->GetParError(1));
    par_err_vec.push_back(0.0);
    par_err_vec.push_back(0.0);

    par_chi2_vec.push_back(f_sigLT->GetChisquare());
    par_chi2_vec.push_back(f_sigLT->GetChisquare());
    par_chi2_vec.push_back(f_sigLT->GetChisquare());
    par_chi2_vec.push_back(f_sigLT->GetChisquare());

'''

'''

1) Iterate for TT params

    c1->cd(4)->SetLeftMargin(0.12);
    n1->Draw("tt:u:tt_e", "", "goff");

    TF1* f_sigTT_pre = new TF1("sig_TT", fun_Sig_TT, 0, 0.5, 2);
    f_sigTT_pre->SetParameters(tt0, tt1);

    TGraphErrors* g_sigtt = new TGraphErrors(n1->GetSelectedRows(), n1->GetV2(), n1->GetV1(), 0, n1->GetV3());

    for(int i =0; i < w_vec.size(); i++) {

            Float_t sigtt_X_pre;
            Float_t q2_term;
            q2_term = tt2 / logq2_vec[i] + tt3 * g_sigtt->GetX()[i] / logq2_vec[i];
            q2_dep = q2_vec[i];

            sigtt_X_pre = (f_sigTT_pre->Eval(g_sigtt->GetX()[i]) / q2_dep + q2_term) * g_vec[i]  * sin(th_vec[i] * pi / 180) * sin(th_vec[i] * pi / 180); 

            g_sigtt_prv->SetPoint(i, n1->GetV2()[i], sigtt_X_pre);

            Float_t sigtt_X_fit, sigtt_X_fit_err;

            if( th_vec[i]  == 180) { 

                    sigtt_X_fit     = 0.0;
                    sigtt_X_fit_err = 1.0;

            } else {

                    sigtt_X_fit     = (g_sigtt->GetY()[i] / g_vec[i] / sin(th_vec[i] * pi / 180) / sin(th_vec[i] * pi / 180) - q2_term) * q2_dep;
                    sigtt_X_fit_err = g_sigtt->GetEY()[i] / g_vec[i] / sin(th_vec[i] * pi / 180) / sin(th_vec[i] * pi / 180) * q2_dep;

            }

            g_sigtt_fit->SetPoint(i, g_sigtt->GetX()[i], sigtt_X_fit);
            g_sigtt_fit->SetPointError(i, 0, sigtt_X_fit_err);

    }

2)	g_max = g_sigtt->GetYaxis()->GetXmax();
    gp_max = TMath::MaxElement(g_sigtt_prv->GetN(), g_sigtt_prv->GetY());

    g_min = g_sigtt->GetYaxis()->GetXmin();
    gp_min = TMath::MinElement(g_sigtt_prv->GetN(), g_sigtt_prv->GetY());

    difff = (g_max - g_min)/5;

    if (g_max < gp_max) {
            g_sigtt->SetMaximum(gp_max + difff);
    } 

    if (g_min > gp_min) {
            g_sigtt->SetMinimum(gp_min - difff);
    } 

3)	g_sigtt->SetTitle("Sig TT");

    g_sigtt->SetMarkerStyle(5);
    g_sigtt->Draw("AP");

    g_sigtt->GetXaxis()->SetTitle("#it{-u} [GeV^{2}]");
    g_sigtt->GetXaxis()->CenterTitle();
    g_sigtt->GetYaxis()->SetTitle("#left(#frac{#it{d#sigma}}{#it{du}}#right)_{TT} [#mub/GeV^{2}]");
    g_sigtt->GetYaxis()->SetTitleOffset(1.5);
    g_sigtt->GetYaxis()->SetTitleSize(0.035);
    g_sigtt->GetYaxis()->CenterTitle();

    g_sigtt_prv->SetMarkerColor(4);
    g_sigtt_prv->SetMarkerStyle(21);
    g_sigtt_prv->Draw("P");

    c2->cd(4);
    g_sigtt_fit->SetTitle("Sigma TT Model Fit");
    g_sigtt_fit->Draw("A*");

    TF1* f_sigTT = new TF1("sig_TT", fun_Sig_TT, 0, 0.5, 2); 
    f_sigTT->SetParameters(tt0, tt1);

4) 	g_sigtt_fit->Fit(f_sigTT);

    for(int i =0; i < w_vec.size(); i++) {

            Float_t sigtt_X;
            Float_t q2_term;
            q2_term = tt2 / logq2_vec[i] + tt3 * g_sigtt->GetX()[i] / logq2_vec[i];

            q2_dep = q2_vec[i];

            sigtt_X = (f_sigTT->Eval(g_sigtt->GetX()[i]) / q2_dep + q2_term) * g_vec[i] * sin(th_vec[i] * pi / 180) * sin(th_vec[i] * pi / 180); 
            g_sigtt_fit_tot->SetPoint(i, g_sigtt->GetX()[i], sigtt_X);

    }

5)	fit_status->DrawTextNDC(0.35, 0.8, " Fit Status: " + gMinuit->fCstatu);

    c1->cd(4);

    g_sigtt_fit_tot->SetMarkerStyle(26);
    g_sigtt_fit_tot->SetMarkerColor(2);
    g_sigtt_fit_tot->SetLineColor(2);
    g_sigtt_fit_tot->Draw("LP");

    c1->Print("x_sep/seperated_" + q2_set + ".png");
    c2->Print("x_sep/seperated_" + q2_set + "_fit.png");

    c1->Print("x_sep/seperated_" + q2_set + ".C");
    c2->Print("x_sep/seperated_" + q2_set + "_fit.C");

6)	tt0 = f_sigTT->GetParameter(0);
    tt1 = f_sigTT->GetParameter(1);

    par_vec.push_back(tt0);
    par_vec.push_back(tt1);
    par_vec.push_back(tt2);
    par_vec.push_back(tt3);

    par_err_vec.push_back(f_sigTT->GetParError(0));
    par_err_vec.push_back(f_sigTT->GetParError(1));
    par_err_vec.push_back(0.0);
    par_err_vec.push_back(0.0);

    par_chi2_vec.push_back(f_sigTT->GetChisquare());
    par_chi2_vec.push_back(f_sigTT->GetChisquare());
    par_chi2_vec.push_back(f_sigTT->GetChisquare());
    par_chi2_vec.push_back(f_sigTT->GetChisquare());
'''

'''
1) Output to fit_params/par.pl_{q2_set} as...
          file_out << setw(12) << par_vec[i]  << "   " << setw(15) << setprecision(4) << par_err_vec[i]  << setw(12) << par_chi2_vec[i] << setw(5) << i << endl;

'''

'''

1) L to T ratio


    TGraphErrors* lt_ratio = new TGraphErrors();

    for(int i =0; i < w_vec.size(); i++) {

            l_t_ratio = f_sigL->Eval(g_sigt->GetX()[i]) / f_sigT->Eval(g_sigt->GetX()[i]);

            l_t_ratio_err =  pow( pow(0.03/ g_sigl_fit_tot->GetY()[i],2) + pow(0.03/g_sigt_fit_tot->GetY()[i], 2), 0.5) * l_t_ratio; 

            cout << "ads   "<< g_sigl_fit_tot->GetY()[i] <<  "   " << g_sigt_fit_tot->GetY()[i] << "    " << l_t_ratio << "   " << l_t_ratio_err << endl;

    }

'''

'''

Functions for each sig...

///*--------------------------------------------------*/
/// Function for SigT

Double_t fun_Sig_T(Double_t *x, Double_t *par) {
   Float_t xx = x[0];
//   Double_t f = par[0]*exp(par[1]*xx) + par[2];


///*--------------------------------------------------*/
/// Straight Line
//   Double_t f = par[0]*xx + par[1];

   Double_t f =  par[0] + par[1]*xx;


   return f;
}



///*--------------------------------------------------*/
/// Function for SigL

Double_t fun_Sig_L(Double_t *x, Double_t *par) {
   Float_t xx = x[0];
//   Double_t f = par[0] + par[1]*xx;

//   Double_t f = par[0]*exp(par[1]*xx) + par[2]/xx;
//   Double_t f = par[0]*(xx-par[1])**2 + par[2];
//   Double_t f = par[0];


//   Double_t f = par[0]*exp(par[1]*xx) + par[2];

///*--------------------------------------------------*/
/// Straight Line
///   Double_t f = par[0]*xx + par[1];

//   Double_t f = par[0]*xx + par[1];

   Double_t f =  par[0] + par[1]*xx;

   return f;
}


///*--------------------------------------------------*/
/// Function for SigLT

Double_t fun_Sig_LT(Double_t *x, Double_t *par) {

   Float_t xx = x[0];

//  Double_t f = par[0]*exp(par[1]*xx) + par[2]/xx;
//  Double_t f = par[0]*exp(par[1]*xx) + par[2]/xx;


///*--------------------------------------------------*/
/// Straight Line
//   Double_t f = par[0]*xx + par[1];

//   Double_t f = par[0]*xx + par[1];

   Double_t f =  par[0] + par[1]*xx;

 //  Double_t f = par[0];
   return f;
}


///*--------------------------------------------------*/
/// Function for SigTT

Double_t fun_Sig_TT(Double_t *x, Double_t *par) {

   Float_t xx = x[0];

//  Double_t f = par[0];
//  Double_t f = par[0]*exp(par[1]*xx) + par[2]/xx;

//   Double_t f = par[0]*exp(-par[1]*xx) + par[2];

//   Double_t f = par[0]*exp(par[1]*xx) + par[2]/xx;


///*--------------------------------------------------*/
/// Straight Line
//   Double_t f = par[0]*xx + par[1];
//   Double_t f = par[0]*xx + par[1];

   Double_t f =  par[0] + par[1]*xx;

   return f;
}


'''

# Define constants
pi = 3.14159
mtar_gev = 0.93827231
mpipl=0.139570
mkpl=0.493677

hi_bound =  0.7;
lo_bound = -0.1;

def x_fit_in_t():
    lt_ratio_file = open("x_sep/LT_ratio.txt", "w")

    single_setting("160", lt_ratio_file)
    single_setting("245", lt_ratio_file)

    lt_ratio_file.close()

def single_setting(closest_date, q2_set):
    prv_par_vec = []
    g_vec = []
    w_vec = []
    q2_vec = []
    th_vec = []
    logq2_vec = []
    par_vec = []
    par_err_vec = []
    par_chi2_vec = []

    t0, t1, t2, t3 = 0, 0, 0, 0
    l0, l1, l2, l3 = 0, 0, 0, 0
    lt0, lt1, lt2, lt3 = 0, 0, 0, 0
    tt0, tt1, tt2, tt3 = 0, 0, 0, 0

    fit_status = TText()

    fn_sep = "x_sep.pl_" + q2_set

    nsep = TNtuple("nsep", "nsep", "t:t_e:l:l_e:lt:lt_e:tt:tt_e:chi:u:u_min:t:w:q2")
    nsep.ReadFile(fn_sep)

    file_path = CACHEPATH + closest_date + "par.pl_" + q2_set
    try:
        with open(file_path, 'r') as para_file_in:
            # Continue with your file processing code here
            # You can use para_file_in.readlines() to read lines from the file
            # Make sure to handle file processing based on Python syntax and conventions
    except FileNotFoundError:
        print(f"File {file_path} not found.")

