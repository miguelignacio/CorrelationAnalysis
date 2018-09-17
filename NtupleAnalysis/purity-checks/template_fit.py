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
        return disterr


def getHistAndErr(df, var, binEdges):
    hist = np.histogram(df[var], bins=binEdges, weights=df['weights'])[0]
    errs = []
    for left, right in zip(binEdges[:-1], binEdges[1:]):
        weights = [w for (x, w) in zip(df[var], df['weights']) if x >= left and x <= right]
        errs.append(np.sqrt(np.sum(np.square(weights))))
    return hist, errs


class TemplateFit:
    def __init__(self, data, dataerr, signal, signalerr, bkg, bkgerr, binEdges, fitRange=None, verbosity=0):
        self.data = np.array(data, dtype='f')
        self.inputSignal = np.array(signal, dtype='f')
        self.inputBkg = np.array(bkg, dtype='f')

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
        def Chi2(N, f):
            model = N * (f * self.signal + (1 - f) * self.bkg)
            totalerror = np.sqrt(np.square(self.dataerr) + np.square(N * (1 - f) * self.bkgerr))
            return np.sum(np.power(np.divide(self.data - model, totalerror,
                                             out=np.zeros_like(totalerror), where=totalerror != 0), 2.0)[slice(*self.fitRange)])

        mt = iminuit.Minuit(Chi2, N=np.sum(self.data), f=0.12, error_N=1, error_f=1,
                            errordef=1, print_level=self.verbosity)
        mt.migrad()

        self.fitN = mt.values['N']
        self.fitNerr = mt.errors['N']
        self.fitf = mt.values['f']
        self.fitferr = mt.errors['f']

        self.fitSignal = self.fitN * self.fitf * self.signal
        self.fitSignalerr = self.fitN * self.fitf * self.signalerr
        self.fitBkg = self.fitN * (1 - self.fitf) * self.bkg
        self.fitBkgerr = self.fitN * (1 - self.fitf) * self.bkgerr

        fitTotal = self.fitSignal + self.fitBkg
        self.residuals = np.divide(fitTotal - self.data, self.dataerr,
                                   out=np.zeros_like(self.dataerr), where=self.dataerr != 0)
        self.chi2 = Chi2(self.fitN, self.fitf)
        self.dof = len(self.data[slice(*self.fitRange)]) - 3

    def getPurity(self, purityMin, purityMax):
        purity, pmin, pmax = getPurity(self.signal, self.bkg, self.binEdges, self.fitf, purityMin, purityMax, True)
        puritylow = getPurity(self.signal, self.bkg, self.binEdges, self.fitf - self.fitferr, purityMin, purityMax)
        purityhigh = getPurity(self.signal, self.bkg, self.binEdges, self.fitf + self.fitferr, purityMin, purityMax)
        if self.verbosity == 1:
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
