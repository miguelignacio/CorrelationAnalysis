import sys
sys.path.append('/mnt/c/Users/marratia/Linux/PionHadronCorr/CompareHistograms')
from AtlasCommonUtils import *
from Legend import Legend
import ROOT
from ROOT import TLatex
SetAtlasStyle()
latex = TLatex()
latex.SetNDC()
latex.SetTextSize(0.04)


def Format(histo):
   histo.GetYaxis().SetNdivisions(5)
   histo.SetMarkerStyle(1)
   histo.Sumw2()

def DrawLabel():
    latex.DrawLatex(0.38,0.97, "#color[2]{NOT CORRECTED WITH MIXED EVENTS}")
    latex.DrawLatex(0.45,0.90, "8 < p_{T}^{cluster} < 16 GeV")
    latex.DrawLatex(0.45,0.85, "0.55 < DNN <0.85")
    latex.DrawLatex(0.45,0.80,"0.2< z_{T}< 0.4")


file = ROOT.TFile("fout.root")
dPhi_iso = file.Get("dPhi_iso_ztmin2_ztmax4")
dPhi_noniso = file.Get("dPhi_noniso_ztmin2_ztmax4")
dPhi_iso_large = file.Get("dPhiLargedEta_iso_ztmin2_ztmax4")
dPhi_noniso_large = file.Get("dPhiLargedEta_noniso_ztmin2_ztmax4")


def DrawHisto(histo, isolabel,etalabel,name):
    c = ROOT.TCanvas()
    Format(histo)
    histo.Draw("e1x0")
    DrawLabel()
    latex.DrawLatex(0.45,0.75, isolabel)
    latex.DrawLatex(0.45,0.70,etalabel)
    c.SaveAs(name)

DrawHisto(dPhi_noniso,  "Iso > 5 GeV", "|#Delta#eta|< 0.6", "noniso.png")
DrawHisto(dPhi_iso,  "Iso < 2 GeV", "|#Delta#eta|< 0.6", "iso.png")
DrawHisto(dPhi_noniso_large,  "Iso > 5 GeV", "0.8 <|#Delta#eta| <1.4", "noniso_large.png")
DrawHisto(dPhi_iso_large,  "Iso < 2 GeV", "0.8 <|#Delta#eta| <1.4", "iso_large.png")
