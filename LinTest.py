from ROOT import TF1, TH1F, gPad, TCanvas, TFile, TGraph, kRed, kBlack, TMultiGraph, TLegend, kBlue, kGreen, kMagenta, kOrange
import numpy as np
import array as arr
import codecs
import math
import scipy 
import json
import os 
import natsort

configfile = 'config.json'
with open(configfile) as datafile:
    config = json.load(datafile)

def plot(datafiles): #Ploeaza histogramele generate din fisierele de tip list pe care le scoate Janus
    
    histsCh0 = []
    histsCh32 =[]
    dataFiles = [os.path.relpath(x) for x in os.scandir(datafiles)]
    dataFiles = natsort.natsorted(dataFiles)
    Ch0 = TCanvas('c1', "Channel 0")
    Ch32 = TCanvas('c2', 'Channel 32')
    for  i in range(len(dataFiles)):
        print(i)
        with codecs.open(dataFiles[i], 'r', 'utf-8', errors='ignore') as file:      
            histCh0 = TH1F(f"h0_{i}", f"Histogram Channel 0_{i}_HoldDelay_" , 2**13, 0, 2**13-1)
            histCh32 = TH1F(f"h32_{i}", f"Histogram Channel 32_{i}_HoldDelay", 2**13, 0, 2**13-1)
            for i, line in enumerate(file):
                if (i>8):
                    line = line.strip()
                    event = line.split()
                    if event[-3] == '00':
                        histCh0.Fill(int(event[-2]))
                    if event[-3] == '32':
                        histCh32.Fill(int(event[-2]))
            histsCh0.append(histCh0)
            histsCh32.append(histCh32)
            # print(f'Max hist {i} = ', histCh0.GetMaximum(), f'Max hist {i} = ', histCh32.GetMaximum())
    print(len(histCh32), " ", len(histCh0))
    

    file = TFile('hists3.root', 'RECREATE')

    for hist in histsCh0:
        if hist.Integral():
            hist.Write()
    for hist in histsCh32:
        if hist.Integral():
            hist.Write()

    for i in range(len(histCh0)):   
        Ch0.cd()
        histsCh0[i].Draw('hist')
        Ch32.cd()
        histsCh32[i].SetLineColor(2)
        histsCh32[i].Draw('hist')
        gPad.WaitPrimitive('ggg')

def getPeak(hist, start, end): #Scoate niste parametrii preliminari pentru a fita gausiene
    mean = 0
    height = 0
    fhwm = 0
    sigma = 0
    baseline = np.array([])
    par = np.array([]), 
    print ('start = ', start, "end = ", end)
    for i in range(start, start + 10):
        baseline = np.append(baseline, hist.GetBinContent(i))

    hist.GetXaxis().SetRange(start, end)
    maxHist = int(hist.GetMaximum())
    for i in range(start, end):
        height = int(hist.GetBinContent(i))
        if (height == maxHist):
            mean = hist.GetBinCenter(i)
            break
    par = np.append(par, height)
    par = np.append(par, mean)
    print("mean = ", mean, 'height = ', height)
    for i in range(hist.FindBin(mean), end):
       
        h = hist.GetBinContent(i)
        c = hist.GetBinCenter(i)
     
        uplim = height/2 + 0.15 * (height/2)
        downlim = height/2 - 0.15 * (height/2)
        if downlim < h and h < uplim:
            fhwm = (c - mean) * 2
            sigma = fhwm/2.355 
            # par = np.append(par, sigma)
            # par = np.append(par, fhwm) 
            break
    print('fhwm = ', fhwm, 'sigma = ', sigma)
    return par

def analyseFiles(datafiles): #Scoate parametrii doriti din fisierele de tip list de la Janus
    voltages = np.array([0.1, 0.11, 0.12, 0.13, 0.14, 0.15, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8,
     0.9, 1.0, 1.1, 1.2, 1.3, 1.35, 1.36, 1.37, 1.38, 1.39, 1.4])
    parCh0 = np.array([])
    parCh32 = np.array([])
 
    histsCh0 = []
    histsCh32 =[]
    dataFiles = [os.path.relpath(x) for x in os.scandir(datafiles)]
    dataFiles = natsort.natsorted(dataFiles)

    c1 = TCanvas("Linearity test", "Linearity Test")
    graphFile = TFile('graphs.root', 'RECREATE')

    graph = TGraph()
    graph.SetName('Linearity Test Ch0')
    graph.GetXaxis().SetTitle("Voltage (V)")
    graph.GetYaxis().SetTitle("Channel")
    graph.SetMarkerSize(1.2)
    graph.SetMarkerStyle(20)
    graph.SetMarkerColor(kBlue)

    graph1 = TGraph()
    graph1.SetName('Linearity Test Ch32')
    graph1.GetXaxis().SetTitle("Voltage (V)")
    graph1.GetYaxis().SetTitle("Channel")
    graph1.SetMarkerSize(1.2)
    graph1.SetMarkerStyle(20)
    graph1.SetMarkerColor(kRed)

    mg = TMultiGraph('graphs', "Linearity Test")
    for  i in range(len(dataFiles)):
        print(i)
        with codecs.open(dataFiles[i], 'r', 'utf-8', errors='ignore') as file:

            histCh0 = TH1F(f"h0_{i}", f"Histogram Channel 0_{i}_HoldDelay" , 2**13, 0, 2**13-1)
            histCh32 = TH1F(f"h32_{i}", f"Histogram Channel 32_{i}_HoldDelay" , 2**13, 0, 2**13-1)
            for j, line in enumerate(file):
                if (j>8):
                    line = line.strip()
                    event = line.split()
                    if event[-3] == '00':
                        histCh0.Fill(int(event[-2]))
                    if event[-3] == '32':
                        if (i > 17):
                            continue
                        histCh32.Fill(int(event[-2]))
        histsCh0.append(histCh0)
        histsCh32.append(histCh32)
        par0 = getPeak(histCh0, histCh0.FindBin(config['start']), histCh0.FindBin(config['end']))
        par32 = getPeak(histCh32, histCh32.FindBin(config['start']), histCh32.FindBin(config['end']))
        print(par32)
        graph.AddPoint(voltages[i], par0[1])
        if (i <= 17):
            graph1.AddPoint(voltages[i], par32[1])
        print(par0)
    mg.Add(graph)
    mg.Add(graph1)
    c1.cd()
    mg.Draw('ap')
    gPad.WaitPrimitive('ggg')
    c1.Write()
    return parCh0



# plot(config['dataFiles'])
analyseFiles(config['dataFiles'])