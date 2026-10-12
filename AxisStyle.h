#ifndef AXISSTYLE_H
#define AXISSTYLE_H

#include "TStyle.h"
#include "TH1.h"

void SetAxis2D(TH1* h) {
    if (!h) return;
    /*
    TObject *oldTitle = h->GetListOfFunctions()->FindObject("title");
    if (oldTitle) {
        h->GetListOfFunctions()->Remove(oldTitle);
        delete oldTitle;
    }*/

    gPad->SetLeftMargin(0.20);
    gPad->SetBottomMargin(0.16);
    gPad->SetTopMargin(0.12);

    gStyle->SetTitleFontSize(0.2);
    gStyle->SetTitleY(1.1);
    h->GetXaxis()->SetTitleFont(43);
    h->GetXaxis()->SetTitleSize(21);
    h->GetXaxis()->SetTitleOffset(0.75);
    h->GetXaxis()->SetLabelSize(0.05);

    h->GetYaxis()->SetTitleFont(43);
    h->GetYaxis()->SetTitleSize(22);
    h->GetYaxis()->SetTitleOffset(1.0);
    h->GetYaxis()->SetLabelSize(0.06);
}

void SetAxis(TH1* h) {
    if (!h) return;

    gPad->SetLeftMargin(0.16);
    gPad->SetBottomMargin(0.16);

    gStyle->SetTitleFontSize(0.06);
    h->GetXaxis()->SetTitleFont(43);
    h->GetXaxis()->SetTitleSize(20);
    h->GetXaxis()->SetTitleOffset(0.76);
    h->GetXaxis()->SetLabelSize(0.055);

    h->GetYaxis()->SetTitleFont(43);
    h->GetYaxis()->SetTitleSize(20);
    h->GetYaxis()->SetTitleOffset(1.5);
    h->GetYaxis()->SetLabelSize(0.05);
}

void SetAxisCF(TH1* h) {
    if (!h) return;

    gPad->SetLeftMargin(0.17);
    gPad->SetBottomMargin(0.18);

    gStyle->SetTitleFontSize(0.07);
    h->GetXaxis()->SetTitleFont(43);
    h->GetXaxis()->SetTitleSize(20);
    h->GetXaxis()->SetTitleOffset(0.9);
    h->GetXaxis()->SetLabelSize(0.045);

    h->GetYaxis()->SetTitleFont(43);
    h->GetYaxis()->SetTitleSize(20);
    h->GetYaxis()->SetTitleOffset(1.0);
    h->GetYaxis()->SetLabelSize(0.045);
}

#endif
