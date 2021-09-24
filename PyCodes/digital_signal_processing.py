import numpy as np

class DigitalSignalProcessing():

    def __init__(self, scalars, t, dt, sampling_rate, f_window=None):
        super().__init__()
        self.scalars = scalars
        self.t = t
        self.dt = dt
        self.sampling_rate = sampling_rate
        self.inverse_rate = 1 / sampling_rate
        self.bias = np.nanmean(scalars)
        self.signal = scalars - self.bias
        if f_window is None:
            self.windowed_signal = np.copy(self.signal)
        else:
            self.windowed_signal = self.signal * f_window(scalars.size)
        self.period_fmap = {
            'fourier' : self.get_period_by_fourier,
            'auto-correlation' : self.get_period_by_auto_correlation,
            'power spectral density' : self.get_period_by_power_spectral_density,
            'boundary crossings' : self.get_period_by_boundary_crossings}

    def get_period_by_fourier(self):
        ftrans = np.fft.rfft(self.windowed_signal)
        loc = np.argmax(np.abs(ftrans))
        freqs = np.fft.rfftfreq(self.windowed_signal.size, self.inverse_rate)
        lag = 1 / freqs[loc]
        return lag

    def get_period_by_auto_correlation(self):
        acf = np.correlate(self.windowed_signal, self.windowed_signal, 'full')[-self.windowed_signal.size:]
        inflection = np.diff(np.sign(np.diff(acf))) ## second-order difference
        peaks = (inflection < 0).nonzero()[0] + 1 ## second-order difference < 0
        delay = peaks[acf[peaks].argmax()] ## delay at maximum peak
        return delay

    def get_period_by_power_spectral_density(self):
        acf = np.correlate(self.windowed_signal, self.windowed_signal, 'full')[-self.windowed_signal.size:]
        pdg = np.fft.rfft(acf)
        freqs = np.fft.rfftfreq(self.windowed_signal.size, self.dt)
        lag = 1 / freqs[np.argmax(pdg)] ## delay at maximum peak
        return lag

    def get_period_by_boundary_crossings(self):
        ## get consecutive differences
        dr = np.diff(self.windowed_signal)
        ## where product of consecutive differences is negative (ie, up + down OR down + up)
        condition_a = (dr[:-1] * dr[1:] < 0)
        ## where absolute consecutive differences are above tol
        tol = np.sqrt(np.std(self.windowed_signal))
        condition_b = (np.abs(dr) > tol)
        ## indices of condition
        condition = (condition_a & condition_b)
        loc = np.where(condition)[0] + 1
        ## initialize coordinates and classifications
        coordinates, classifications = [], []
        if len(loc) > 1:
            ## get coordinates of hills and valleys
            for r, t in zip(self.windowed_signal[loc], self.t[loc]):
                if r >= 0:
                    extrema_type = 'max'
                else:
                    extrema_type = 'min'
                coordinates.append(t) # .append((r, t))
                classifications.append(extrema_type)
            coordinates, classifications = np.array(coordinates), np.array(classifications)
            ## estimate periods
            periods = []
            for ti, tj, ci, cj in zip(coordinates[:-1], coordinates[1:], classifications[:-1], classifications[1:]):
                dt = tj - ti
                period = 2*dt if ci == cj else dt
                periods.append(period)
            period = np.nanmedian(periods) # np.nanmean(periods)
        else:
            period = np.nan
        return period

    def get_period(self, method):
        if not isinstance(method, str):
            raise ValueError("invalid type(method): {}".format(type(method)))
        if method not in list(self.period_fmap.keys()):
            raise ValueError("invalid method: {}".format(method))
        get_period = self.period_fmap[method]
        return get_period()













##
