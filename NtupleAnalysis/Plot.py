#%matplotlib inline
import ROOT
ROOT.gSystem.Load("/mnt/c/Users/marratia/Linux/RooUnfold/libRooUnfold")
from ROOT import gRandom, TH1, TH1D, cout, THStack
#for unfolding
from ROOT import RooUnfoldResponse
from ROOT import RooUnfold
from ROOT import RooUnfoldBayes
from ROOT import RooUnfoldSvd
from ROOT import RooUnfoldTUnfold
from ROOT import RooUnfoldBinByBin
from ROOT import RooUnfoldIds

#plotting, style
from matplotlib import pyplot as plt
ROOT.gStyle.SetOptStat('')
from AtlasCommonUtils import SetAtlasStyle

from Legend import Legend
from ROOT import TLatex
SetAtlasStyle()
ROOT.gStyle.SetPalette(ROOT.kViridis);

#def myText(x,  y, color, text):#
#
#  ROOT.TLatex l; //l.SetTextAlign(12); l.SetTextSize(tsize);
#  l.SetNDC();
#  l.SetTextColor(color);
#  l.DrawLatex(x,y,text);
#}



def my_callback(option, opt, value, parser):
  setattr(parser.values, option.dest, value.split(','))

import optparse,os
op = optparse.OptionParser(usage="%prog [opts] ")

op.add_option("","--dataname",dest="dataname",
              default=None,
              help="enter data name")

op.add_option("","--ymax",dest="ymax",
              default=0.0,
              help="y max of histogram")

op.add_option("","--purity",dest="purity",
              default=0.0,
              help="purity from 0 to 1.0")

op.add_option("","--rebin",dest="rebin",
              default=1,
              help="argument for TH1::Rebin()")

op.add_option("","--hists",type='string', dest = "hists",
              default ='dPhi',
              action='callback', callback=my_callback,
              help="comma separated list of histograms you would like to Compare, e.g. hists,hists_i,hists_o,hists_ss" )


(p_ops, p_args) = op.parse_args()


dataname = p_ops.dataname


filepath = '/mnt/c/Users/marratia/Linux/PionHadronCorr/NtupleAnalysis/'

labeltitle = ''

Datafile = None 
Datafile = ROOT.TFile('%s%s'%(filepath,dataname),'READ')
if(Datafile==None):
    print 'error, datafile not found'


if '17q' in dataname :
    labeltitle = '5 TeV pp data'
elif '13' in dataname:
    labeltitle = '5 TeV p+Pb data'

if 'JET_R0.3' in dataname:
    jet_R = '0.3'
elif 'JET_R0.4' in dataname:
    jet_R = '0.4'

print 'Trying splitting'
jet_ptmin = float(dataname.split("_PT_",1)[1].split('_Iso')[0]) 
cluster_min = float(dataname.split("_clus",1)[1].split('_max')[0])
cluster_max = float(dataname.split("_max",1)[1].split('_JET')[0])

print 'jetptmin', jet_ptmin
print 'cluster_min', cluster_min
print 'cluster_max', cluster_max

camino = '/mnt/c/Users/marratia/Linux/PionHadronCorr/NtupleAnalysis/'

##if '13' in dataname:
#MCfilename = camino+ 'Out_17g6a1_clus%2.1f_max%2.1f_JET_R%s_PT_%2.1f_Iso1.5.root'%(cluster_min, cluster_max, jet_R,jet_ptmin)
#elif '17q' in dataname:
if '17q' in dataname:
    MCfilename = camino+ 'Out_Skimmed_18b10a_pthat1to6_Lambda_clus%2.1f_max%2.1f_JET_R%s_PT_%2.1f_Iso1.5.root'%(cluster_min, cluster_max, jet_R,jet_ptmin)
else: 
    #MCfilename = camino+  'Out_17g6a1_Lambda_clus%2.1f_max%2.1f_JET_R%s_PT_%2.1f_Iso1.5.root'%(cluster_min, cluster_max, jet_R,jet_ptmin)
    MCfilename = camino+ 'Out_Skimmed_18b10a_pthat1to6_Lambda_clus%2.1f_max%2.1f_JET_R%s_PT_%2.1f_Iso1.5.root'%(cluster_min, cluster_max, jet_R,jet_ptmin)

    
print 'MC file used: ' , MCfilename
MCfile    = ROOT.TFile(MCfilename, 'READ')
MCfile.Print()

responseMatrix = MCfile.Get("h_Xj_Matrix");
#responseMatrix.Rebin2D(nrebin,nrebin)#
#responseMatrix.ClearUnderflowAndOverflow()
mc_truth = responseMatrix.ProjectionX()
mc_reco = responseMatrix.ProjectionY()
mc_truth.SetLineColorAlpha(2,0.8)
mc_reco.SetLineColorAlpha(4,0.8)

##draw mc truth, mc reco, response matrix
c = ROOT.TCanvas('c','c',1600,600)
c.Divide(3)
c.cd(1)
responseMatrix.Draw('colz')
responseMatrix.SetTitle('; x^{true}_{J} ; x^{reco}_{J} ')

ROOT.gPad.SetLogz()
c.cd(2)

label = Legend("")

label.Add(mc_reco,'MC Reco','L')
label.Add(mc_truth,'MC True','L')

hs_mc = ROOT.THStack()
hs_mc.Add(mc_reco)
hs_mc.Add(mc_truth)
hs_mc.Draw('nostack')
hs_mc.SetTitle(';  x^{reco}_{J}; counts')

label.Draw(0.5,.87)
c.cd(3)
ratio_mc = mc_truth.Clone()
ratio_mc.SetMinimum(0)
ratio_mc.SetLineColor(1)
ratio_mc.Divide(mc_reco)
ratio_mc.Draw()
ratio_mc.SetTitle(' ; x^{truth}; Truth/Reco')
c.Draw()
c.SaveAs('UnfoldingMatrix%s.png'%(dataname))

def SetHistoStyle(h,color=2,alpha=0.5):
    h.SetMarkerColorAlpha(color,alpha)
    h.SetLineColorAlpha(color,alpha)
    h.SetMarkerStyle(20)
    h.SetMarkerSize(1)
    #h.SetFillColorAlpha(color,alpha)

def GetHistos(variable='dPhi',rebin=1,purity=0.25):
    
    hweights = Datafile.Get('h_Ntriggers')
    #hweights.Print()
    w_SR = hweights.GetBinContent(1)
    w_BR = hweights.GetBinContent(2)

    print ' SIGNAL REGION ' , w_SR , ' BKG REGION' , w_BR

    SR = Datafile.Get('hSR_'+variable)
    BR = Datafile.Get('hBR_'+variable)
    NonIso = Datafile.Get('hNonIso_'+variable)

    #BR.Print()
    #print 'Integral BR' , BR.Integral()
    
    hMCweights = MCfile.Get('h_weights')
    MC = MCfile.Get('hSR_'+variable)

    #for h in [SR,BR,MC]:
    #    h.Rebin(rebin)



    BR_scaled = BR.Clone('BR_scaled')
    BR_scaled.Scale(1.0-purity)
    
    Sub = SR.Clone()
    Sub.Add(BR_scaled, -1.0)
    Sub.Scale(1.0/purity)
    
    SetHistoStyle(SR, 2, 1.0)
    SetHistoStyle(BR, ROOT.kBlue, 0.70)
    SetHistoStyle(MC, ROOT.kOrange+8, 0.7)
    SetHistoStyle(BR_scaled, 4, 0.7)
    SetHistoStyle(Sub, ROOT.kOrange+8, 0.85)

    SetHistoStyle(NonIso, 1, 0.8)
    histos = [SR, BR, BR_scaled, Sub, MC, hweights, NonIso]
    for h in histos:
        h.Rebin(rebin)
        #h.Scale(1.0/h.GetBinWidth(1))
        
    return histos


#def ObtainSystematicVariation(purity, rebin, variable):
#    histo = {}
#    for key in purity.keys():
#        h = GetHistos(variable,rebin=rebin, purity=purity[key])
#        histo[key] = h[3]
#        hweights = h[5]
#        hweights.Print()
#        print 'key ' , key
#        for ibin in range(1,h[3].GetNbinsX()+1):
#            print key, ' i' , ibin, ' ' , histo[key].GetBinContent(ibin)
#    histo['systematics'] = histo['nominal'].Clone(histo['nominal'].GetName()+'systematic')
#    for ibin in range(1,histo["nominal"].GetNbinsX()+1):
#        maxvariation =  max( abs(histo["nominal"].GetBinContent(ibin) - histo["up"].GetBinContent(ibin)) ,  abs(histo["nominal"].GetBinContent(ibin) - histo["do"].GetBinContent(ibin)))
#        histo['systematics'].SetBinError(ibin, ROOT.TMath.Sqrt(maxvariation*maxvariation + histo['nominal'].GetBinError(ibin)*histo['nominal'].GetBinError(ibin)))
#        print 'i ' , ibin ,' nominal: ' ,histo["nominal"].GetBinContent(ibin) , ' MAXVARIATION: ' , maxvariation
    
#    return histo['systematics']

canvas = ROOT.TCanvas('canvas','canvas')

def Plotting(histos, purity,  variable,rebin=1, title='', ymax = None):
    print ' Entering plotting macro'
    
    hSR, hBR , hBR_scaled, hSub, hMC , hweights, hNonIso = histos
    
    hs = ROOT.THStack()
    hBR.SetTitle(title)
    hSR.SetTitle(title)
    hSub.SetTitle(title)

    hBR.SetMinimum(0)
    hs.Add(hBR)
    #hBR.Draw()
    #hSR_dPhihBR.SetMaximum(hBR.GetMaximum()*1.2)
    hBR.GetYaxis().SetNdivisions(8)
    hBR.GetXaxis().SetNdivisions(8)
    hSR.GetYaxis().SetNdivisions(8)
    hSR.GetXaxis().SetNdivisions(8)

    #hBR_scaled.Draw('samehist')
    hs.Add(hBR_scaled)

    #hSR.Draw('same')
    hs.Add(hSR)
    hs.Draw("nostack")
    label = Legend(labeltitle)
    label.Add(hSR, "Signal Region","P")
    label.Add(hBR, "Bkg Region","P")
    label.Add(hBR_scaled, "Bkg Region, scaled by purity","L")

    label.Draw(0.45,0.90)
    plotname = variable+'_step1_purity%2.3f'%purity['nominal']
    canvas.SaveAs('IntermediateStep_'+plotname+"_DATANAME_%s.pdf"%dataname)
    canvas.Clear()

    ########################
    hMC.SetFillColorAlpha(ROOT.kOrange, 0.70) 
    #hMC.Draw('hist')
    hs = ROOT.THStack()
    hs.Add(hSub)
    hs.Add(hBR)
    hs.Add(hNonIso)
  
    #hs.Add(hBR_scaled,'hist')
    hs.Draw('nostack')
    if(ymax>0):
        hs.SetMaximum(ymax)
    elif(variable=='jetpt'):
        #hs.SetMaximum(0.025)
        hs.SetMinimum(0)

    if(variable=='jetpt'):
        hs.SetMinimum(0.0)
        hs.GetXaxis().SetRangeUser(10,30)
    hs.SetTitle(title)
    hs.GetYaxis().SetNdivisions(8)
    hs.GetXaxis().SetNdivisions(8)
 
    #hSys = ObtainSystematicVariation(purity, rebin, variable)  
    #hs.Add(hSys) 
    #hSys.SetFillColorAlpha(ROOT.kGray,0.5)
    #hSys.SetFillStyle(3144)
   
    #hSys.Draw("E2same")
    #hSub.Draw("same")
    #hSys.Draw("E2same")
    #hMC.Draw("samehist")

    label = Legend(labeltitle)    
    label.Add(hSub, "Iso #gamma-jet","P")
    #label.Add(hMC, "Pythia LO #gamma-jet","L")
    label.Add(hBR, "Iso #pi^{0}-jet ","P")
    label.Add(hNonIso, "#pi^{0}-jet", "P")

    label.Draw(.15,.90)

    latex = TLatex()   
    latex.SetNDC()
    latex.SetTextColor(1);
    xlabel = 0.65
    if(variable=='dPhi'): xlabel = 0.45

    latex.DrawLatex(xlabel,0.88, "%2.0f < p_{T}^{trig} < %2.0f GeV"%(cluster_min,cluster_max))
    latex.DrawLatex(xlabel,0.81, "R=%s, p_{T}^{jet} > %2.0f GeV"%(jet_R,jet_ptmin))
    if(variable !='dPhi'): latex.DrawLatex(xlabel,0.74, "#Delta#phi > #pi/2")
    chi2perdot = hSub.Chi2Test(hBR,"P CHI2/NDF WW")
    pvalue = hSub.Chi2Test(hBR,"P WW ")
    #if(variable !='dPhi'): latex.DrawLatex(xlabel,0.67, "#chi^{2}/dof =%2.1f, p=%2.2f "%(chi2perdot,pvalue))    
    if(variable !='dPhi'): latex.DrawLatex(xlabel,0.67, "p value = %2.2f "%(pvalue))
    canvas.SaveAs('Final_'+plotname+"_DATANAME_%s.pdf"%dataname)
    for i in range(hSub.GetNbinsX()): 
        print  'BinContent1', hSub.GetBinContent(i)
    hSub.Draw()
 
    print 'BinContent1', hSub.GetBinContent(2)
    hSub.GetYaxis().SetNdivisions(8)
    hSub.GetXaxis().SetNdivisions(8)
    hMC.SetFillColorAlpha(ROOT.kGray, 0.50)
    hMC.SetLineColorAlpha(ROOT.kGray, 0.50)
    hMC.SetMarkerColorAlpha(ROOT.kGray, 0.50)

    if(ymax>0):
        hSub.SetMaximum(ymax)
    elif(variable=='jetpt'):
        #hSub.SetMaximum(0.025)
        #hSub.SetMinimum(0)
        print 'oli'
    

    if(variable=='jetpt'): 
        hSub.SetMinimum(0.0)
        hSub.GetXaxis().SetRangeUser(10,30)
    hMC.Draw("histsame")
    hMC.Draw("same")
    hSub.Draw("same")
    latex.DrawLatex(xlabel,0.88, "%2.0f < p_{T}^{trig} < %2.0f GeV"%(cluster_min,cluster_max))
    latex.DrawLatex(xlabel,0.81, "R=%s, p_{T}^{jet} > %2.0f GeV"%(jet_R,jet_ptmin))
    if(variable !='dPhi'): latex.DrawLatex(xlabel,0.74, "#Delta#phi > #pi/2")
    label = Legend(labeltitle)
    label.Add(hSub, "#gamma-jet","P")
    label.Add(hMC, "Pythia LO #gamma-jet","L")
    #label.Add(hBR, "#pi^{0}-jet ","P")
    label.Draw(.2,.90)

    chi2perdot = hSub.Chi2Test(hMC,"P CHI2/NDF WW")
    pvalue = hSub.Chi2Test(hMC,"P WW") 
    #latex.DrawLatex(xlabel,0.67, "#chi^{2}/dof =%2.1f, p=%2.2f "%(chi2perdot,pvalue))   
    if(variable !='dPhi'): latex.DrawLatex(xlabel,0.67, "p value = %2.2f "%(pvalue))
    if(variable=='jetpt_all'):
       ROOT.gPad.SetLogy()
    #fout = ROOT.TFile("finalplots_"+plotname+"%s"%dataname,"RECREATE")
    #hSub.Write("%s_subtracted"%variable)
    #hBR.Write("%s_Bkg"%variable)
    #hMC.Write("%s_MC"%variable)
    #fout.Close()
    canvas.SaveAs('Final_'+plotname+"_DATANAME_%s_WithMC.pdf"%dataname)
    canvas.Clear()   
    SetHistoStyle(hSub,1,1.0)
    SetHistoStyle(hMC,1,1.0)
    SetHistoStyle(hBR,1,1.0)
    
    return hSub,hMC,hBR, hNonIso

##Systematic uncertainty 
purity = {}
purity["nominal"] = float(p_ops.purity)
purity["up"]      = float(p_ops.purity)
purity["do"]      = float(p_ops.purity)

titles = {}
titles['dPhi'] = '; Azimuthal difference #Delta#phi [rad]; 1/N_{#gamma} dN_{pairs}/#Delta#phi'
titles['dEta'] = '; #Delta\eta; 1/N_{#gamma} dN_{pairs}/#Delta#eta'
titles['Multiplicity'] = '; Particle multiplicity in Jet; 1/N_{#gamma} dN_{pairs}/dM'
titles['pTD']  = '; Jet p_{T}D ; 1/N_{#gamma} dN_{pairs}/dp_{T}D'
titles['jetwidth'] = '; Jet #sigma ;  1/N_{#gamma} dN_{pairs}/d#sigma'
titles['Xj']       = '; X_{J} = p_{T}^{jet}/p_{T}^{#gamma}; 1/N_{#gamma} dN_{pairs}/dX_{j}'
titles['XobsPb']  = '; X_{obs}^{Pb}   = (p_{T}^{#gamma}e^{-#eta^{#gamma}} + p_{T}^{jet}e^{-#eta^{jet}})/E_{Pb}; 1/N_{#gamma} dN_{pairs}/dX_{obs}^{Pb}'
titles['AvgEta'] =  '; (#eta^{#gamma} + #eta^{jet})/2.0 ; 1/N_{#gamma} dN_{pairs}/d#Delta#eta'
titles['njet']   =  ' ; number of associated jets ;  1/N_{#gamma} dN_{pairs}'
titles['clusterpt'] = ' ; cluster p_{T} [GeV] ; Normalized units'
titles['clustereta'] = ' ; cluster eta ; Normalized units'
titles['jetpt']      = ' ; Reconstructed jet p_{T} [GeV] ; 1/N_{#gamma} dN_{pairs}/dp_{T}'
titles['jetpt_all']      = ' ; Reconstructed jet p_{T} [GeV] ; 1/N_{#gamma} dN_{pairs}/dp_{T}'
titles['jetphi']      = ' ; jet #phi  [rad] ; 1/N_{#gamma} dN_{pairs}/d#phi'
for i in range(len(p_ops.hists)):
    histname = p_ops.hists[i]
    hSub,hMC,hBR, hNonIso= Plotting(GetHistos(histname,rebin=int(p_ops.rebin), purity=purity['nominal']), purity, histname, rebin=1, title=titles[histname], ymax=float(p_ops.ymax))

    fout = ROOT.TFile("finalplots_"+histname+"%s"%dataname,"RECREATE")
    
    hSub.SetMarkerColor(1)
    hSub.SetLineColor(1)
    hBR.SetMarkerColor(1)
    hBR.SetLineColor(1)
    hMC.SetMarkerColor(1)
    hMC.SetLineColor(1)

    hSub.Write("%s_subtracted"%histname)
    hBR.Write("%s_Bkg"%histname)
    hMC.Write("%s_MC"%histname)
    hNonIso.Write("%s_Pion"%histname)
    fout.Close()

#Plotting(GetHistos('Multiplicity',rebin=2, purity=purity['nominal']), purity, 'Multiplicity', rebin=1, title= '; Particle multiplicity in Jet; 1/N_{#gamma} dN/dM')
#Plotting(GetHistos('dEta',rebin=4, purity=purity['nominal']), purity, 'dEta', rebin=1, title=  '; #Delta\eta_{#gamma-jet}; 1/N_{#gamma} dN/#Delta#eta_{#gamma-jet}')
#Plotting(GetHistos('AvgEta',rebin=4, purity=purity['nominal']), purity, 'AvgEta', rebin=1, title=  '; \eta_{#gamma-jet} = (\eta_{#gamma} + \eta_{jet})/2.0; 1/N_{#gamma} dN/\eta_{#gamma-jet}')
#Plotting(GetHistos('jetpt',rebin=5, purity=purity['nominal']), purity, 'jetpt',  rebin=1, title='; jet p_{T} [GeV] ; 1/N_{#gamma} dN/dp_{T}')
#Plotting(GetHistos('jetphi',rebin=1, purity=purity['nominal']) , purity, 'jetphi', rebin=1, title='; jet #phi [rad] ; 1/N_{#gamma} dN/d#phi')
#Plotting(GetHistos('jeteta',rebin=1, purity=purity['nominal']), purity, 'jeteta', rebin=1, title= '; jet #eta ; 1/N_{#gamma} dN/d#eta')
#Plotting(GetHistos('pTD',rebin=4, purity=purity['nominal']), purity, 'pTD', rebin=1, title='; Jet p_{T}D ; 1/N_{#gamma} dN/dp_{T}D')
#Plotting(GetHistos('Xj',rebin=1, purity=purity['nominal']), purity, 'Xj', rebin=1, title='; X_{J} = p_{T}^{jet}/p_{T}^{#gamma}; 1/N_{#gamma} dN/dX_{j}')
#Plotting(GetHistos('njet',rebin=1, purity=purity['nominal']), purity, 'njet', rebin=1, title= '; Number of associated jets n_{jet}; 1/N_{#gamma} dN/dn_{jet}')
#Plotting(GetHistos('clustereta',rebin=4, purity=purity['nominal']), purity, 'clustereta', rebin=1, title= '; cluster #eta;  1/N_{#gamma} dN/#eta')
#Plotting(GetHistos('clusterphi',rebin=1, purity=purity['nominal']), purity, 'clusterphi', rebin=1, title= '; cluster #phi;  1/N_{#gamma} dN/#phi')


#
#SR, BR, BR_scaled, Sub, MC, hweights, NonIso = GetHistos('clusterpt',rebin=1, purity=purity['nominal'])


#c = ROOT.TCanvas()

#SR.Rebin(2)
#BR.Rebin(2)
#SR.SetTitle('; cluster p_{T} [GeV]; weighted counts')
#SR.DrawNormalized()
#BR.DrawNormalized('same')

#label = Legend(labeltitle)
#label.Add(SR, "Signal region","L")
#label.Add(BR, "Bkg region","L")
#label.Draw(.55,.85)
#c.SaveAs('mcpt.pdf')
