#!/usr/bin/env python

import sys
import ROOT
from argparse import ArgumentParser
import subprocess as sp
from collections import OrderedDict, namedtuple
from math import sqrt

#-------------------------------------------------------------------------------

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)

plots = OrderedDict()

plots[("etaPhi", "sagitta")]   = OrderedDict()
#plots[("etaPhi", "absBetaMI")] = OrderedDict()
#plots[("etaPhi", "absBetaOI")] = OrderedDict()
#plots[("etaPhi", "absBetaOM")] = OrderedDict()

#plots[("tower", "sagitta")]   = OrderedDict()
#plots[("tower", "absBetaMI")] = OrderedDict()
#plots[("tower", "absBetaOI")] = OrderedDict()
#plots[("tower", "absBetaOM")] = OrderedDict()

COLORS = [ROOT.kRed, ROOT.kBlue, ROOT.kGreen, ROOT.kOrange]
LABELS = ["naive +/- 1.5 sigma", "iterative", "MCP", "double Gaussian"]

rootFile = ROOT.TFile.Open("out.root")

def fetchPlots() :
    for (regionType, deflectionType), plotDict in plots.iteritems() :
        for fitStage in ["deflection", "phiMod", "eta"] :
            #regionType = k[0]
            #deflectionType = k[1]
            #plotName = "h_" + regionType + "_" + deflectionType + "_" + fitStage + "_ptRelDelta"
            plotName = "h_{}_{}_{}_ptRelDelta".format(regionType, deflectionType, fitStage)
            plot = rootFile.Get(plotName)
            #plots[(regionType, deflectionType)].append(plot)
            plotDict[plotName] = plot

#-------------------------------------------------------------------------------

class ResolutionData :
    def __init__(self) :
        self.mean = self.dMean = 0.0
        self.stdDev = self.dStdDev = 0.0
        self.chiSquare = 0.0
        self.meanPlot = ROOT.TGraphErrors()
        self.stdDevPlot = ROOT.TGraphErrors()
        self.chiSquarePlot = ROOT.TGraph()
        self.fit = None

    def updatePlots(self, pt) :
        n = self.meanPlot.GetN()
        self.meanPlot.SetPoint(n, pt, self.mean)
        self.meanPlot.SetPointError(n, 0.25, self.dMean)
        self.stdDevPlot.SetPoint(n, pt, self.stdDev)
        self.stdDevPlot.SetPointError(n, 0.25, self.dStdDev)
        self.chiSquarePlot.SetPoint(n, pt, self.chiSquare)

    def allPlots(self) :
        return [self.meanPlot, self.stdDevPlot, self.chiSquarePlot]


#-------------------------------------------------------------------------------
# pT resolution techniques (naive, iterative, double Gaussian)

def naiveResolution(pt, ptBinHisto, resData, stdDevMultiple=1.5) :
    mean = ptBinHisto.GetMean()
    stdDev = ptBinHisto.GetStdDev()
    #fitMin = mean - args.fitRange * stdDev
    #fitMax = mean + args.fitRange * stdDev
    fitMin = mean - stdDevMultiple * stdDev
    fitMax = mean + stdDevMultiple * stdDev
    fName = "f_{}".format(ptBinHisto.GetName())
    f = ROOT.TF1(fName, "gaus(0)", fitMin, fitMax)
    #f = ROOT.TF1(fName,  "[0]*exp(-0.5*((x-[1])/[2])**2 + [3]*exp(-0.5*((x-[4])/[5])**2", fitMin, fitMax)
    f.SetParameters(1.0, 0.0, 1.0)
    fitOption = "" if abs(args.fitRange) < 1e-6 else "R"
    ptBinHisto.Fit(f, fitOption)
    params = f.GetParameters()
    errors = f.GetParErrors()
    resData.fit = f

    resData.mean = params[1]
    resData.dMean = errors[1]
    resData.stdDev = abs(params[2])
    resData.dStdDev = errors[2]
    resData.chiSquare = f.GetChisquare()/f.GetNDF()
    resData.updatePlots(pt)


def iterativeResolution(pt, ptBinHisto, resData, maxIter=100) :
    # iterate to try to obtain stable RMS upon rescaling
    converged = False
    for i in xrange(maxIter) :
        previousRMS = ptBinHisto.GetRMS()
        ptBinHisto.GetXaxis().SetRangeUser(-2.0*previousRMS, 2.0*previousRMS)
        newRMS = ptBinHisto.GetRMS()

        if abs((newRMS - previousRMS)/previousRMS) < 0.1 :
            converged = True
            break

    if not converged :
        print "Unable to converge on stable RMS after {} iterations...skipping...".format(maxIter)
        return None

    # do Gaussian fit
    fName = "f_{}".format(ptBinHisto.GetName())
    f = ROOT.TF1(fName, "gaus(0)", -2.0*newRMS, 2.0*newRMS)
    f.SetParameters(1.0, 0.0, 1.0)
    ptBinHisto.Fit(f, "R")
    params = f.GetParameters()
    errors = f.GetParErrors()
    resData.fit = f

    resData.mean = params[1]
    resData.dMean = errors[1]
    resData.stdDev = abs(params[2])
    resData.dStdDev = errors[2]
    resData.chiSquare = f.GetChisquare()/f.GetNDF()
    resData.updatePlots(pt)


def mcpResolution(pt, ptBinHisto, resData) :
    rms = ptBinHisto.GetRMS()
    ptBinHisto.GetXaxis().SetRangeUser(-2.0*rms, 2.0*rms)
    distributionMean = ptBinHisto.GetMean()

    # Gaussian fit with mean = distribution mean
    fName = "f_{}".format(ptBinHisto.GetName())
    #f = ROOT.TF1(fName, "[0]*exp(-0.5*((x-{})/[1])**2)".format(distributionMean), -1.5*rms, 1.5*rms)
    f = ROOT.TF1(fName, "gaus".format(distributionMean), -1.5*rms, 1.5*rms)
    f.FixParameter(1, distributionMean)
    f.SetParameters(1.0, 1.0)
    ptBinHisto.Fit(f, "")
    params = f.GetParameters()
    errors = f.GetParErrors()
    resData.fit = f

    resData.mean = distributionMean
    resData.dMean = abs(params[1])/sqrt(ptBinHisto.GetEntries())
    resData.stdDev = abs(params[1])
    resData.dStdDev = errors[1]
    resData.chiSquare = f.GetChisquare()/f.GetNDF()
    resData.updatePlots(pt)


def doubleGaussianResolution(pt, ptBinHisto, resData) :
    mean = 0.0
    rms = ptBinHisto.GetRMS()
    #ptBinHisto.GetXaxis().SetRangeUser(-4.0*rms, 2.0*rms)
    fName = "f_{}".format(ptBinHisto.GetName())
    f = ROOT.TF1(fName, "gaus(0) + gaus(3)", -1.5*rms, 1.5*rms)
    f.SetParameters(1.0, 0.0, 1.0, 1.0, 0.0, 3.0)
    ptBinHisto.Fit(f, "R")
    params = f.GetParameters()
    errors = f.GetParErrors()
    resData.fit = f

    resData.mean = params[1]
    resData.dMean = errors[1]
    resData.stdDev = abs(params[2])
    resData.dStdDev = errors[2]
    resData.chiSquare = f.GetChisquare()/f.GetNDF()
    resData.updatePlots(pt)

#-------------------------------------------------------------------------------

def drawPtSlice(canvas, ptBinHisto, *resData) :
    canvas.cd()

    leg = ROOT.TLegend(0.65, 0.72, 0.93, 0.90)
    leg.SetFillStyle(0)
    leg.SetFillColor(0)
    leg.SetBorderSize(0)
    leg.SetTextFont(42)
    leg.SetNColumns(1)

    ptBinHisto.Draw()
    for i, f in enumerate(resData) :
        f.SetLineColor(COLORS[i])
        leg.AddEntry(f, LABELS[i], "L")
        f.Draw("same")
    leg.Draw()
    canvas.SaveAs("test.png")


#-------------------------------------------------------------------------------

if __name__ == "__main__" :
    parser = ArgumentParser()
    parser.add_argument("-s", "--sampleFit", action="store_true", help="output fits using different fitting strategies")
    parser.add_argument("-r", "--fitRange", type=float, default=0.0, help="number of standard deviations w.r.t. mean to fit or 0 (default) for full range")
    args = parser.parse_args()

    fetchPlots()
    #outFile = ROOT.TFile("resolution_plots.root", "recreate")

    nProjHist = -1
    for (regionType, deflectionType), plotDict in plots.iteritems() :
        print "Making resolution plot for {}, {}".format(regionType, deflectionType)
        nameSuffix = "{}_{}".format(regionType, deflectionType)
        cName = "c_" + nameSuffix
        #canvas = ROOT.TCanvas(cName, cName, 800, 800)
        #canvas = ROOT.TCanvas(cName, cName)
        canvas = ROOT.TCanvas(cName, "multipads", 900, 700)
        canvas.Divide(1, 3, 0.0, 0.0)
        canvas.cd(1)

        naiveData = ResolutionData()
        iterData  = ResolutionData()
        mcpData   = ResolutionData()
        dgData    = ResolutionData()

        histList = []
        for plotName, plot in plotDict.iteritems() :
            #print "{}, {} : {}".format(regionType, deflectionType, plot.GetEntries())
            #ptProfile = plot.ProfileX("profile")
            #ptProfile.Write("test")
            hName = "{}_proj".format(plotName)
            #resolutionHisto = plot.ProjectionX(hName)
            if "_eta_" in plotName :
                print "Adding plot: {}".format(plotName)
                histList.append(plot)

            #print "Drawing `{}`...".format(hName)
            #resolutionHisto.Draw(drawOpt)
            #drawOpt = "same"
            #stack.Add(resolutionHisto)

            if "_eta_" in plotName :  # make graph only for final fit
                # set points in `resolutionPlot` for each bin
                for i in xrange(1, plot.GetNbinsX()) :
                    #if i < plot.GetNbinsX() - 1 : continue

                    ptMin = 10.0 + 0.5*i
                    pt = ptMin + 0.25

                    nProjHist += 1
                    projHistName = "h_{}".format(nProjHist)
                    ptBinHisto = plot.ProjectionY(projHistName, i, i)

                    naiveHisto = ptBinHisto.Clone()
                    naiveResolution(pt, naiveHisto, naiveData)

                    iterHisto = ptBinHisto.Clone()
                    iterativeResolution(pt, iterHisto, iterData)

                    mcpHisto = ptBinHisto.Clone()
                    mcpResolution(pt, mcpHisto, mcpData)

                    dgHisto = ptBinHisto.Clone()
                    doubleGaussianResolution(pt, dgHisto, dgData)

                    print
                    print "chiSquarePerNDF values: {}, {}, {}, {}".format(naiveData.chiSquare, iterData.chiSquare, mcpData.chiSquare, dgData.chiSquare)
                    print

                    if args.sampleFit :
                        # draw fit for one pT slice, then exit
                        drawPtSlice(naiveData, iterData, mcpData)
                        sys.exit()

        canvas.cd(1)
        leg = ROOT.TLegend(0.65, 0.72, 0.93, 0.90)
        leg.SetFillStyle(0)
        leg.SetFillColor(0)
        leg.SetBorderSize(0)
        leg.SetTextFont(42)
        leg.SetNColumns(1)
        #for i, p in enumerate([resolutionPlot_naive, resolutionPlot_iter, resolutionPlot_mcp, resolutionPlot_dg]) :

        # axis labels
        naiveData.chiSquarePlot.GetXaxis().SetTitle("p_T^\\text{off} \\quad \\text{[GeV]}")  # only bottom plot
        naiveData.meanPlot.GetYaxis().SetTitle("\\frac{p_{T}^\\text{on} - p_T^\\text{off}}{p_T^\\text{off}}")
        naiveData.stdDevPlot.GetYaxis().SetTitle("std. dev.")
        naiveData.chiSquarePlot.GetYaxis().SetTitle("\\chi^2/NDF")

        meanPlots = [naiveData.meanPlot, iterData.meanPlot, mcpData.meanPlot, ]
        minMean = min([p.GetHistogram().GetMinimum() for p in meanPlots])
        maxMean = max([p.GetHistogram().GetMaximum() for p in meanPlots])
        for i, p in enumerate(meanPlots) :
            p.SetLineColor(COLORS[i])
            p.SetFillColor(COLORS[i])
            p.SetMarkerColor(COLORS[i])
            leg.AddEntry(p, LABELS[i], "L")
            if i == 0 :
                drawOpt = "ap"
                p.GetXaxis().SetRangeUser(10.0, 40.0)
                meanRange = maxMean - minMean
                p.GetYaxis().SetRangeUser(minMean - 0.1*meanRange, maxMean + 0.1*meanRange)
            else :
                drawOpt = "p"
            p.Draw(drawOpt)
        leg.Draw()

        canvas.cd(2)
        stdDevPlots = [naiveData.stdDevPlot, iterData.stdDevPlot, mcpData.stdDevPlot, ]
        minStdDev = min([p.GetHistogram().GetMinimum() for p in stdDevPlots])
        maxStdDev = max([p.GetHistogram().GetMaximum() for p in stdDevPlots])
        for i, p in enumerate(stdDevPlots) :
            p.SetLineColor(COLORS[i])
            p.SetFillColor(COLORS[i])
            p.SetMarkerColor(COLORS[i])
            if i == 0 :
                drawOpt = "ap"
                p.GetXaxis().SetRangeUser(10.0, 40.0)
                stdDevRange = maxStdDev - minStdDev
                p.GetYaxis().SetRangeUser(minStdDev - 0.1*stdDevRange, maxStdDev + 0.1*stdDevRange)
            else :
                drawOpt = "p"
            p.Draw(drawOpt)

        canvas.cd(3)
        chiSquarePlots = [naiveData.chiSquarePlot, iterData.chiSquarePlot, mcpData.chiSquarePlot, ]
        minChiSquare = min([p.GetHistogram().GetMinimum() for p in chiSquarePlots])
        maxChiSquare = max([p.GetHistogram().GetMaximum() for p in chiSquarePlots])
        for i, p in enumerate(chiSquarePlots) :
            p.SetLineColor(COLORS[i])
            p.SetFillColor(COLORS[i])
            p.SetMarkerColor(COLORS[i])
            p.SetMarkerStyle(4)
            if i == 0 :
                drawOpt = "ap"
                p.GetXaxis().SetRangeUser(10.0, 40.0)
                chiSquareRange = maxChiSquare - minChiSquare
                p.GetYaxis().SetRangeUser(minChiSquare - 0.1*chiSquareRange, maxChiSquare + 0.1*chiSquareRange)
            else :
                drawOpt = "p"
            p.Draw(drawOpt)

        canvas.SaveAs("test.png")
        break


