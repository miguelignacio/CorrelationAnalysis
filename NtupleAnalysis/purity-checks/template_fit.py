import matplotlib.pyplot as plt
import numpy as np

import iminuit

def haveSameLength(*args):
    n = len(args[0])
    return all(len(l) == n for l in args)


def normalize(x):
    return np.array(x, dtype='f') / np.sum(x)


def getBinRange(binEdges, valuemin, valuemax):
    binmin = min([i for i, edge in enumerate(binEdges) if edge >= valuemin])
    binmax = max([i for i, edge in enumerate(binEdges) if edge <= valuemax])
    return binmin, binmax


def getPurity(signal, bkg, binEdges, frac, purityMin, purityMax, returnRange=False):
    # signal and bkg should be normalized to 1
    pmin, pmax = getBinRange(binEdges, purityMin, purityMax)
    purity = np.sum(frac * signal[pmin:pmax]) / (np.sum([(1 - frac) * bkg[pmin:pmax], frac * signal[pmin:pmax]]))

    if returnRange:
        return purity, pmin, pmax
    else:
        return purity


def getErr(dist, disterr):
    if disterr is None:
        return np.sqrt(dist)
    else:
        return np.array(disterr, dtype='f')


def getHistAndErr(df, var, binEdges):
    hist = np.histogram(df[var], bins=binEdges, weights=df['weights'])[0]
    errs = []
    for left, right in zip(binEdges[:-1], binEdges[1:]):
        weights = [w for (x, w) in zip(df[var], df['weights']) if x >= left and x <= right]
        errs.append(np.sqrt(np.sum(np.square(weights))))
    return hist, errs


def getNormHistAndErr(df, fitvar, binEdges):
    hist, err = getHistAndErr(df, fitvar, binEdges)
    nTotal = np.sum(hist)
    return np.divide(hist, nTotal), np.divide(err, nTotal)


class TemplateFit:
    def __init__(self, data, dataerr, signal, signalerr, bkg, bkgerr, binEdges, fitRange=None, verbosity=0):
        self.data = np.array(data, dtype='f')
        self.inputSignal = np.array(signal, dtype='f')
        self.inputBkg = np.array(bkg, dtype='f')
        
        self.N = sum(self.data)

        self.dataerr = getErr(data, dataerr)
        self.inputSignalerr = getErr(signal, signalerr)
        self.inputBkgerr = getErr(bkg, bkgerr)

        self.binEdges = binEdges
        if fitRange:
            self.fitRange = getBinRange(binEdges, *fitRange)
        else:
            self.fitRange = (0, None)

        if not haveSameLength(self.data, self.dataerr, self.inputSignal, self.inputSignalerr, self.inputBkg, self.inputBkgerr, self.binEdges[1:]):
            raise ValueError('Inputs do not have the same length (binEdges should have 1 more than the rest)')

        self.signal = self.inputSignal / np.sum(self.inputSignal)
        self.signalerr = self.inputSignalerr / np.sum(self.inputSignal)
        self.bkg = self.inputBkg / np.sum(self.inputBkg)
        self.bkgerr = self.inputBkgerr / np.sum(self.inputBkg)
        self.binCenters = np.array([(hedge + ledge) / 2.0 for ledge, hedge in zip(binEdges[:-1], binEdges[1:])])
        self.binWidths = np.array([hedge - ledge for ledge, hedge in zip(binEdges[:-1], binEdges[1:])])

        self.signalColor = '#3B7EA1'
        self.bkgColor = '#FDB515'
        self.figureSize = (10, 10)

        self.verbosity = verbosity

        self.doFit()

    def doFit(self):
        def Chi2(f):
            model = self.N * (f * self.signal + (1 - f) * self.bkg)
            totalerror = np.sqrt(np.square(self.dataerr) + np.square(self.N * (1 - f) * self.bkgerr))
            return np.sum(np.power(np.divide(self.data - model, totalerror,
                                             where=np.array(totalerror)!=0), 2.0)[slice(*self.fitRange)])

        mt = iminuit.Minuit(Chi2, f=0.12, error_f=1, errordef=1, print_level=self.verbosity)
        mt.migrad()

        self.fitf = mt.values['f']
        self.fitferr = mt.errors['f']

        self.fitSignal = self.N * self.fitf * self.signal
        self.fitSignalerr = self.N * self.fitf * self.signalerr
        self.fitBkg = self.N * (1 - self.fitf) * self.bkg
        self.fitBkgerr = self.N * (1 - self.fitf) * self.bkgerr

        fitTotal = self.fitSignal + self.fitBkg
        self.residuals = np.divide(fitTotal - self.data, self.dataerr, where=np.array(self.dataerr)!=0)
        self.chi2 = Chi2(self.fitf)
        self.dof = len(self.data[slice(*self.fitRange)]) - 1

    def getChi2DofInRange(self, rangeMin, rangeMax):
        binmin, binmax = getBinRange(self.binEdges, rangeMin, rangeMax)
        dof = len(self.residuals[binmin:binmax]) - 1
        chi2 = np.sum(np.power(self.residuals[binmin:binmax], 2.0))
        
        return chi2/dof

    def getPurity(self, purityMin, purityMax, verbosity=-1):
        if verbosity == -1:
            verbosity = self.verbosity
        
        purity, pmin, pmax = getPurity(self.signal, self.bkg, self.binEdges, self.fitf, purityMin, purityMax, True)
        puritylow = getPurity(self.signal, self.bkg, self.binEdges, self.fitf - self.fitferr, purityMin, purityMax)
        purityhigh = getPurity(self.signal, self.bkg, self.binEdges, self.fitf + self.fitferr, purityMin, purityMax)
        if verbosity == 1:
            print 'Purity = %2.5f, +%2.5f, -%2.5f' % (purity, purityhigh - purity, purity - puritylow)
        return purity, pmin, pmax, puritylow, purityhigh

    def plotFit(self, xlabel='', dataLabel='Data', signalLabel='Signal (MC)', bkgLabel='Bkg (data)', texts=[], legendoptions={}):
        dataplot = plt.errorbar(self.binCenters, self.data, yerr=self.dataerr, label=dataLabel, fmt='ko')
        bkgplot = plt.bar(self.binCenters, self.fitBkg, yerr=self.fitBkgerr, width=self.binWidths,
                          align='center', label=bkgLabel, capsize=0,
                          color=self.bkgColor, ec=self.bkgColor, ecolor=self.bkgColor)
        signalplot = plt.bar(self.binCenters, self.fitSignal, yerr=self.fitSignalerr, bottom=self.fitBkg, width=self.binWidths,
                             align='center', label='{0} + Bkg'.format(signalLabel), capsize=0,
                             color=self.signalColor, ec=self.signalColor, ecolor=self.signalColor)
        chi2text, = plt.plot([], [], ' ', label='Chi2/dof = %2.2f' % (self.chi2 / self.dof))
        textHandles = []
        for text in texts:
            textHandle, = plt.plot([], [], ' ', label=text)
            textHandles.append(textHandle)

        plt.legend(handles=[dataplot, signalplot, bkgplot, chi2text] + textHandles, numpoints=1, loc='best', **legendoptions)
        plt.xlabel(xlabel)
        plt.ylabel('Entries')

    def plotResiduals(self, xlabel, figureFilename=''):
        fig = plt.figure(figsize=self.figureSize)
        plt.plot(self.binCenters, self.residuals, 'o')
        plt.xlabel(xlabel)
        plt.ylabel('(Fit - Data)/Dataerr')

        if figureFilename:
            fig.savefig(figureFilename)

        plt.show()

    def plotTemplates(self, xlabel, figureFilename=''):
        fig = plt.figure(figsize=(self.figureSize[0], self.figureSize[1] / 2.0))
        plt.bar(self.binCenters, self.fitSignal, width=self.binWidths, align='center', label='Signal (MC)',
                color=self.signalColor, ec=self.signalColor, alpha=0.4)
        plt.bar(self.binCenters, self.fitBkg, width=self.binWidths, align='center', label='Background (data)',
                color=self.bkgColor, ec=self.bkgColor, alpha=0.4)
        plt.xlabel(xlabel)
        plt.ylabel('Entries')
        plt.legend(loc='best')

        if figureFilename:
            fig.savefig(figureFilename)

        plt.show()

    def plotNormalizedTemplates(self, xlabel, figureFilename=''):
        fig = plt.figure(figsize=(self.figureSize[0], self.figureSize[1] / 2.0))
        plt.bar(self.binCenters, self.signal, width=self.binWidths, align='center', label='Signal (MC)',
                color=self.signalColor, ec=self.signalColor, alpha=0.4)
        plt.bar(self.binCenters, self.bkg, width=self.binWidths, align='center', label='Background (data)',
                color=self.bkgColor, ec=self.bkgColor, alpha=0.4)
        plt.xlabel(xlabel)
        plt.ylabel('Entries')
        plt.legend(loc='best')

        if figureFilename:
            fig.savefig(figureFilename)

        plt.show()

        
class BackgroundFit:
    def __init__(self, data, dataerr, bkg, bkgerr, binEdges, fitRange, verbosity=0):
        self.data = np.array(data, dtype='f')
        self.bkg = np.array(bkg, dtype='f')

        self.dataerr = getErr(data, dataerr)
        self.bkgerr = getErr(bkg, bkgerr)
        
        self.fitRange = fitRange
        self.fitRangeSlice = slice(*getBinRange(binEdges, *fitRange))
        self.binEdges = binEdges
        self.binCenters = np.array([(hedge+ledge)/2.0 for ledge, hedge in zip(binEdges[:-1], binEdges[1:])])
        self.binWidths = np.array([hedge-ledge for ledge, hedge in zip(binEdges[:-1], binEdges[1:])])
        
        self.bkgColor = '#FDB515'
        self.figureSize = (10, 10)
        
        self.doFit(verbosity)
        
    def doFit(self, verbosity):
        def Chi2(N):
            model = N * self.bkg
            totalerr = np.sqrt(np.square(self.dataerr) + np.square(N*self.bkgerr))
            return np.sum(np.power(np.divide(self.data - model, totalerr,
                                             out=np.zeros_like(totalerr), where=totalerr!=0), 2.0)[self.fitRangeSlice])
        
        mt = iminuit.Minuit(Chi2, N=np.sum(self.data)/np.sum(self.bkg), error_N=1, errordef=1, print_level=verbosity)
        mt.migrad()

        self.fitN = mt.values['N']
        self.fitNerr = mt.errors['N']
        self.fitBkg = self.fitN*self.bkg
        self.fitBkgerr = self.fitN*self.bkgerr
        self.residuals = np.divide(self.fitBkg - self.data, self.dataerr,
                                   out=np.zeros_like(self.dataerr), where=self.dataerr!=0)
        self.chi2 = Chi2(self.fitN)
        self.dof = len(self.data[self.fitRangeSlice])-1
        
    def getChi2DofInRange(self, rangeMin, rangeMax):
        binmin, binmax = getBinRange(self.binEdges, rangeMin, rangeMax)
        dof = len(self.residuals[binmin:binmax]) - 1
        chi2 = np.sum(np.power(self.residuals[binmin:binmax], 2.0))
        
        return chi2/dof
        
    def getPurity(self, purityMin, purityMax, verbosity=0):
        pmin, pmax = getBinRange(self.binEdges, purityMin, purityMax)
        fitSignal = np.subtract(self.data, self.fitBkg)
        purity = np.sum(fitSignal[pmin:pmax])/np.sum(self.data[pmin:pmax])
        lowSignal = np.subtract(self.data, (self.fitN+self.fitNerr)*self.bkg)
        puritylow = np.sum(lowSignal[pmin:pmax])/np.sum(self.data[pmin:pmax])
        highSignal = np.subtract(self.data, (self.fitN-self.fitNerr)*self.bkg)
        purityhigh = np.sum(highSignal[pmin:pmax])/np.sum(self.data[pmin:pmax])
        
        if verbosity == 1:
            print 'Purity = %2.5f, +%2.5f, -%2.5f' % (purity, purityhigh - purity, purity - puritylow)
        
        return purity, puritylow, purityhigh
    
    def plotFit(self, xlabel='var', purityRange=None, legendheader=None, newFigure=True, figureFilename=None):
        if newFigure:
            fig = plt.figure()
        ax = plt.gca()
        
        handles = []
        if legendheader:
            headertext,  = plt.plot([], [], ' ', label=legendheader)
            handles.append(headertext)
        dataplot = plt.errorbar(self.binCenters, self.data, yerr=self.dataerr, label='Data (iso)', fmt='ko')        
        bkgplot = plt.bar(self.binCenters, self.fitBkg, yerr=self.fitBkgerr, width=self.binWidths,
                          align='center', label='Bkg (data anti-iso)', capsize=0,
                          color=self.bkgColor, ec=self.bkgColor, ecolor=self.bkgColor)
        chi2text, = plt.plot([], [], ' ', label='$\chi^2$/dof = %2.2f'%(self.chi2/self.dof))
        handles.extend([dataplot, bkgplot, chi2text])
        bkgminbin, bkgmaxbin = getBinRange(self.binEdges, *self.fitRange)
        ax.axvspan(self.binEdges[bkgminbin], self.binEdges[bkgmaxbin], color='red', alpha=0.2)
        
        if purityRange:
            purityMin, purityMax = purityRange
            pmin, pmax = getBinRange(self.binEdges, purityMin, purityMax)
            purity, puritylow, purityhigh = self.getPurity(purityMin, purityMax)
            ax.axvspan(self.binEdges[pmin], self.binEdges[pmax], color='black', alpha=0.2)
        
        plt.legend(handles=handles, numpoints=1, loc='best', fontsize=12)
        plt.xlabel(xlabel)
        plt.ylabel('Entries')
        
        if newFigure and figureFilename:
            fig.savefig(figureFilename)
            
    def plotResiduals(self, newFigure=True, figureFilename=None):
        if newFigure:
            fig = plt.figure(figsize=self.figureSize)
        plt.plot([self.binEdges[0], self.binEdges[-1]], [0, 0], 'y--')
        plt.plot(self.binCenters[self.fitRange], self.residuals[self.fitRange], 'o')
        plt.ylabel('(Fit - Data)/Dataerr')
        
        if newFigure and figureFilename:
            fig.savefig(figureFilename)
