#!/usr/bin/env python
from AtlasCommonUtils import *


def my_callback(option, opt, value, parser):
  setattr(parser.values, option.dest, value.split(','))

import optparse,os
op = optparse.OptionParser(usage="%prog [opts] ")

op.add_option("","--outputDir",dest="outdir",
              default=".",
              help="Output Directory")

op.add_option("-i","--inputFileList",dest="inlist",
              default='input.txt',
              help="Input File list with xsections")

op.add_option("-n","--Normalized",dest="Normalized",
              default=False,
              help="Normalized")

op.add_option("-y","--LogY",dest="LogY",
              default=False,
              help="LogY")


op.add_option("-l","--RatioLimit",dest="RatioLimit",
              default=4,
              help="The limit for the ratio plot (in %)")

op.add_option("-t","--Tag",dest="Tag",
              default="",
              help="tag to add to plots")

op.add_option("","--hists",type='string', dest = "hists",
              action='callback', callback=my_callback,
              help="comma separated list of histograms you would like to Compare, e.g. hists,hists_i,hists_o,hists_ss" )


(p_ops, p_args) = op.parse_args()


output_dir = p_ops.outdir
if not os.path.exists(output_dir):
    os.makedirs(output_dir)
else:
    print "Using output dir ",output_dir,", overwriting existing contents"


if not p_ops.inlist:
    raise RuntimeError("Must specify infile list")


## Array to store a list of names of graphs to be compared
hists = []
if not p_ops.hists:
    hists = ["nhits"]
else:
    hists = p_ops.hists


RatioLimit = float(p_ops.RatioLimit)
Tag = p_ops.Tag

print 'The ratio plot will be plotted at maximum' , RatioLimit
import ROOT 
from Legend import Legend
from ROOT import TLatex
SetAtlasStyle()


color = [1,2,4,417,7,8]
style = [20,21,22,20] 
alphas = [1.0, 0.75, 0.6, 0.5]
def DivideGraphs(den,num): 
    x = ROOT.Double() 
    y1= ROOT.Double()
    y2= ROOT.Double()
    ratio = ROOT.TGraphErrors()
    for i in range(den.GetN()):
        num.GetPoint(i,x,y2)
        den.GetPoint(i,x,y1)
        if(y1!=0 and y2!=0):
            #print "SETTING POINT ", y2/y1
            n = ratio.GetN()
            ratio.SetPoint(n , x,y2/y1*1.0)#-1)*100)
            ratio.SetPointError(n,0,0)
    #num.Print()
    #den.Print()
    #ratio.Print()
    return ratio

def TransformTH1toTGraph(histo):
    graph = ROOT.TGraphErrors(histo)
    graph.GetXaxis().SetTitle(histo.GetXaxis().GetTitle())
    graph.GetYaxis().SetTitle(histo.GetYaxis().GetTitle())
    return graph
            
def NormalizeGraph(graph):

    x = ROOT.Double()
    y= ROOT.Double()
    integral = 0
    for i in range(graph.GetN()):
        graph.GetPoint(i, x, y)
        integral = integral + y
  
    if integral==0: return

    for i in range(graph.GetN()):
        graph.GetPoint(i, x, y)
        
        graph.SetPoint(i, x, y/integral)
        graph.SetPointError(i, graph.GetErrorX(i), graph.GetErrorY(i)/integral)

    return graph


def plotHists(hists,files, xtitle, ytitle, title, pdfname,ly=False):
    c = ROOT.TCanvas("c","c",900,600)
    pad1 = ROOT.TPad("pad1","This is pad1",0.0, 0.0, 0.9, 1.0,         0)
    pad2 = ROOT.TPad("pad2","This is pad2",0.85, 0.0, 1.0, 1.0,        0)
    pad1.Draw()
    pad2.Draw()
    pad1.cd()    
    start = 0
    label = Legend("")
    multi = ROOT.TMultiGraph()

    for i in range(len(hists)): #This loops over the files. nhists(f1,f2,f3) is the format of hists.
        temp = files[i].GetName()
        if(p_ops.Normalized): hists[i] = NormalizeGraph(hists[i])     
        hists[i].SetLineColorAlpha(color[i], alphas[i])
        hists[i].SetMarkerColorAlpha(color[i],alphas[i])
        hists[i].SetMarkerStyle(style[i])
        hists[i].SetMarkerSize(1)
        hists[i].SetLineWidth(2)
       # print 'i ' , i , ' MINIMUM:  ' , hists[i].GetHistogram().GetMinimum() , ' MAXIMUM : ' , hists[i].GetHistogram().GetMaximum(), 
        label.Add(hists[i],tag[i],"P")
        multi.Add(hists[i])

    
    multi.Draw("AP")
    ROOT.gStyle.SetErrorX(0.0001);
    multi.SetMaximum(1.02*multi.GetHistogram().GetMaximum())
    multi.SetTitle(title)
    multi.GetXaxis().SetTitle(xtitle)
    multi.GetYaxis().SetTitle(ytitle)
    multi.GetXaxis().SetNdivisions(8)
    multi.GetYaxis().SetNdivisions(8)
    multi.GetXaxis().CenterTitle()
    multi.GetYaxis().CenterTitle()
    multi.SetMinimum(0)
    #if ly==True:
        #multi.SetMinimum(multi.GetHistogram().GetMinimum()+1e-4) 
    if ly: 
        multi.SetMinimum(multi.GetHistogram().GetMinimum()+1e-4) 
        ROOT.gPad.SetLogy(1)
    pad2.cd()
    label.Draw(0.0,.7)
    latex = TLatex()
    latex.SetNDC()
    c.SaveAs(pdfname)
    c.Clear()
    print "END OF PLOTTING FUNCTION " 

def get_hists(files,hists):
    allhists = []
    for h in hists:
        hlist = []
        for f in files:
            f.Print()
            f.cd()
            histo = f.Get(h) 
            #histo.Print()
            if 'TH1' in str(type(histo)): 
                print 'Found TH1, transforming it into TGraphErrors:'
                histo = TransformTH1toTGraph(histo)
                 
            hlist.append(histo)
        allhists.append(hlist)
    return allhists





files = []
tag = []
xsecs = []
legs = []
weight = []



########### START OF MAIN FUNCTION #########################################
print "Opening file file " 
inlist = open(p_ops.inlist)

if(p_ops.LogY):
    print ' You will be using Log scale in Y axis for all histograms ' 


for line in inlist:
    if line[0] == "#": ##skipping comments
       continue
    print line
    tag.append(line.split(",")[1])
    weight.append(1)
    line = line.split(",")[0] #the line after the ROOT file name is the label 
    print 'Opening file', line
    files.append(ROOT.TFile(line.strip("\n")))
    line = line.split("hists")[0].split("_")[-1]
    tagOutputPlots = line
    
print "Getting all hists"
c = ROOT.TCanvas()

allhists = get_hists(files, hists) #gets all histograms 
temp = get_hists(files,hists)
for i in range(len(hists)):
    print 'Xaxis title = ', allhists[i][0].GetXaxis().GetTitle()
    plotHists(temp[i], files, allhists[i][0].GetXaxis().GetTitle(),allhists[i][0].GetYaxis().GetTitle() ,  allhists[i][0].GetTitle(),  output_dir + "/"+hists[i]+Tag+'_Comparison.pdf',ly= p_ops.LogY)









def FormatGraph(graph, index):
    graph.SetLineColorAlpha(color[index], alphas[index])
    graph.SetMarkerColorAlpha(color[index],alphas[index])
    graph.SetMarkerStyle(style[index])
    graph.SetMarkerSize(1)
    graph.SetLineWidth(2)
    

print 'there are a total of ', len(tag)

cRatio = ROOT.TCanvas()





for i in range(len(hists)):
    multi = ROOT.TMultiGraph()
    multi2 = ROOT.TMultiGraph()
    for j in range(len(tag)):
        ref = NormalizeGraph(allhists[i][0])
        num  = NormalizeGraph(allhists[i][j])

   
        ratio =  DivideGraphs(ref, num)
        FormatGraph(ratio,j )
        FormatGraph(num, j)
        multi2.Add(num)
     
       # print 'i ' , i , ' MINIMUM:  ' , hists[i].GetHistogram().GetMinimum() , ' MAXIMUM : ' , hists[i].GetHistogram().GetMaximum(),
        
        multi.Add(ratio)    

    multi.Draw("APL")
    if float(p_ops.RatioLimit) > 0 :
        minimo = 1 -float(p_ops.RatioLimit)
        maximo = 1+ float(p_ops.RatioLimit)
        multi.GetHistogram().SetMaximum(maximo)
        multi.GetHistogram().SetMinimum(minimo)
#    multi.Draw("APL") 
    multi.GetYaxis().SetTitle("ratio")
    multi.GetXaxis().SetTitle(allhists[i][0].GetXaxis().GetTitle())
    
    #ROOT.gPad.SetLogy() 
    cRatio.SaveAs(output_dir + '/'+hists[i]+Tag+'_Ratio.pdf')



